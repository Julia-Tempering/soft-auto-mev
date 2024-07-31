using CSV 
using DataFrames
using Dates
using Distributions
using MCMCChains
using BridgeStan
using Random
using SplittableRandoms: SplittableRandom
using Statistics
using StanSample

using Pigeons
using autoHMC
using autoRWMH

###############################################################################
# loading data
###############################################################################

function base_dir()
    base_folder = dirname(dirname(Base.active_project()))
    endswith(base_folder, "soft-auto-mev") || error("please activate the soft-auto-mev julia-environment")
    return base_folder
end

###############################################################################
# utils for processing samples
###############################################################################

get_component_samples(samples::AbstractVector, idx_component::Int) =
    [s[idx_component] for s in samples]
get_component_samples(samples::DataFrame, idx_component) = Array(samples[:,idx_component])

# returns tuple: (margin_idx, mean, std_dev)
special_margin_mean_std(model::String, args...) = 
    if startswith(model,"funnel")
        (1, 0., 3.)
    elseif startswith(model,"banana")
        (1, 0., sqrt(10))
    elseif startswith(model,"eight_schools")    # approx using variational PT (see 8schools_pilot.jl)
        (1, 3.574118538746056, 3.1726880307401455)
    elseif startswith(model, "normal") 
        (1, 0., 1.)
    elseif startswith(model, "two_component_normal")
        (1, 0., first(two_component_normal_stdevs(args...)) ) # the margin with the largest stdev
    else
        throw(KeyError(model))
    end
n_vars(samples::AbstractVector) = length(first(samples))
n_vars(samples::DataFrame) = size(samples,2)

to_chains(samples::AbstractVector) = Chains(samples)
to_chains(samples::DataFrame) = Chains(Array(samples))

function df_to_vec(df::DataFrame) 
    n = size(df,1)
    df_vec = Vector{Vector{Float64}}(undef, n) 
    for i in 1:n 
        df_vec[i] = Vector(df[i, :]) 
    end 
    return df_vec
end 

###############################################################################
# ESS and friends
###############################################################################

include("batch_means.jl")

function margin_ess(samples, model, args...)
    margin_idx, true_mean, true_sd = special_margin_mean_std(model, args...)
    margin_samples = get_component_samples(samples, margin_idx)
    batch_means_ess(margin_samples, true_mean, true_sd) 
end

# minimum ess using batch means
min_ess_batch(samples) = minimum(1:n_vars(samples)) do i
    batch_means_ess(get_component_samples(samples, i))
end

# minimum ess using the MCMCChains approach (which is based on Stan's)
function min_ess_chains(samples)
    chn = to_chains(samples)
    min(
        minimum(ess(chn).nt.ess),            # default is :bulk which computes ess on rank-normalized vars
        minimum(ess(chn,kind=:basic).nt.ess) # :basic is actual ess of the vars.
    )
end

# minimum over dimensions and methods
function min_ess_all_methods(samples, model)
    special_margin_ess = try
        margin_ess(samples, model)
    catch e
        if e isa KeyError
            Inf
        else
            rethrow(e)
        end
    end
    min(special_margin_ess, min(min_ess_chains(samples), min_ess_batch(samples)))
end

# returns ESS, mean and var
function margin_summary(samples, model, args...)
    margin_idx, true_mean, true_sd = special_margin_mean_std(model, args...)
    margin_samples = get_component_samples(samples, margin_idx)
	margin_ess_exact = batch_means_ess(margin_samples, true_mean, true_sd)
	margin_mean = mean(margin_samples) 
	margin_var = var(margin_samples)
    DataFrame(margin_ess_exact = margin_ess_exact, margin_mean = margin_mean, margin_var = margin_var)
end

###############################################################################
# sampling
###############################################################################

function make_explorer(sampler_str, selector_str, int_time_str, jitter_str)
    jitter = jitter_str == "normal" ? Normal(0.0, 0.5) : Dirac(0.0)
	if sampler_str in ("SliceSampler", "NUTS")                        # irrelevant for NUTS since we use cmdstan
		SliceSampler(n_passes=1)
	elseif sampler_str == "SimpleAHMC"                                # we handle autoMALA as special case of SimpleAHMC
		selector = selector_str == "inverted" ?
			autoHMC.AMSelectorInverted() : autoHMC.AMSelectorLegacy() # legacy matches the Pigeons.AutoMALA behavior
		int_time = if int_time_str == "single_step"
            autoHMC.FixedIntegrationTime()                            # Pigeons.AutoMALA is recovered because additionally jitter is Dirac and selector is legacy, see autoHMC tests
        elseif int_time_str == "fixed"
            autoHMC.AdaptiveFixedIntegrationTime()
        else
            autoHMC.AdaptiveRandomIntegrationTime()
        end
		SimpleAHMC(
			n_refresh=1, int_time = int_time, step_size_selector = selector,
			step_jitter_dist = jitter
		)
	elseif sampler_str == "SimpleRWMH"
		selector = selector_str == "inverted" ?
			autoRWMH.MHSelectorInverted() : autoRWMH.MHSelector()
		SimpleRWMH(step_size_selector = selector, step_jitter_dist = jitter)
	else
		throw(ArgumentError("unknown sampler $sampler_str"))
	end
end

#######################################
# pigeons
#######################################

function pt_sample_from_model(model, target, seed, explorer, miness_threshold; max_rounds = 25)
    n_rounds = 1 # NB: cannot start from >1 otherwise we miss the explorer n_steps from all but last round
    pt = PT(Inputs(
        target      = target, 
        seed        = seed,
        n_rounds    = n_rounds,
        n_chains    = 1, 
        record      = [record_default(); Pigeons.explorer_acceptance_pr; Pigeons.traces; Pigeons.reversibility_rate],
        explorer    = explorer, 
        show_report = true
    ))

    # run until minESS threshold is breached
    n_steps = n_samples = 0
    miness = 0.0
    local samples
    while n_rounds < max_rounds # bail after this point
        pt = pigeons(pt)
        n_steps += first(Pigeons.explorer_n_steps(pt))
        samples = get_sample(pt) # only from last round
        n_samples = length(samples)
        miness = n_samples < miness_threshold ? 0.0 : min_ess_all_methods(samples, model) # skip computing ess for low sample sizes (buggy)
        miness > miness_threshold && break
        @info """
            Low ESS after round $n_rounds: miness = $miness < $miness_threshold = miness_threshold
            Running another round.
        """
        pt = Pigeons.increment_n_rounds!(pt, 1)
        n_rounds += 1 
    end
    time = sum(pt.shared.reports.summary.last_round_max_time) # despite name, it is a vector of time elapsed for all rounds
    return time, samples, n_steps, miness
end

#######################################
# cmdstan
#######################################

# StanSample is broken; it doesn't respect the arguments in the SampleModel object
# Therefore, we use its internals but then bypass it by calling cmdstan directly
function nuts_sample_from_model(model, seed, miness_threshold; max_samples = 2^25, kwargs...)
    # make model and data from the arguments
    stan_model = model_string(model; kwargs...)
    sm = SampleModel(model, stan_model) 
    data = stan_data(model; kwargs...)
    StanSample.update_json_files(sm, data, 1, "data")

    # run until minESS threshold is breached
    n_samples = max(1000, ceil(Int, 10*miness_threshold)) # assume ESS/n_samples = 10%. Also run at least 1000 o.w. cmdstan complains
    n_steps = 0
    miness = time = 0.0
    local samples
    while n_samples < max_samples
        cmd = stan_cmd(sm, n_samples, seed)
        time += @elapsed run(cmd)
        info = DataFrame(CSV.File(joinpath(sm.tmpdir, model * "_chain_1.csv"), comment = "#"))
        n_steps += sum(info.n_leapfrog__) # count leapfrogs (including warmup)
        samples = info[(sm.num_warmups+1):end, 8:end] # discard 7 auxiliary variables + warmup
        miness = min_ess_all_methods(samples, model)
        miness > miness_threshold && break
        @info """
            Low ESS: miness = $miness < $miness_threshold = miness_threshold
            Doubling n_samples and re-running.
        """
        n_samples *= 2 
    end
    return time, samples, n_steps, miness
end

function stan_cmd(sm::SampleModel, n_samples, stan_seed; print_every=max(100, round(Int,n_samples/50)))
    cmd = `$(sm.output_base) num_threads=1 sample num_chains=1 num_samples=$n_samples`
    cmd = `$cmd num_warmup=$n_samples save_warmup=true adapt engaged=true`
    cmd = `$cmd id=1 data file=$(first(sm.data_file)) random seed=$stan_seed`
    cmd = `$cmd output file=$(sm.output_base)_chain_1.csv refresh=$print_every`
    return cmd
end

###############################################################################
# model utils
###############################################################################

# scale parametrization for funnel and banana
# same as for autoMALA paper: larger scale -> easier
function get_scale(scale_idx, model_name)
	@assert 1 ≤ scale_idx ≤ 20
	if model_name == "funnel_scale"
		0.2 * scale_idx
	elseif model_name == "banana_scale"
		e = ((-13):6)[scale_idx]
		2.0 ^ e
	else
		throw(ArgumentError("Don't know how to make a scale argument for $model_name"))
	end
end

function model_string(model; dataset=nothing, kwargs...)
    if model == "normal" # dont have the standard normal on Pigeons examples
        return "data {
          int<lower=1> dim;
        }
        parameters {
          vector[dim] x;
        }
        model {
          x ~ std_normal();
        }"
    end
    if startswith(model, "two_component_normal")
        return read(joinpath(
            base_dir(), "stan", "two_component_normal.stan"), String)
    end
    if startswith(model, "horseshoe")
        is_logit = any(Base.Fix1(startswith,dataset), ("prostate", "ionosphere", "sonar"))
        return read(joinpath(
            base_dir(), "stan", "horseshoe_" * (is_logit ? "logit" : "linear") * ".stan"
        ), String)
    end
    if model == "mRNA"
        return read(joinpath(base_dir(), "stan", "mRNA.stan"), String)
    end
    pigeons_stan_dir = joinpath(dirname(dirname(pathof(Pigeons))),"examples","stan")
    if startswith(model, "eight_schools_") 
        return read(joinpath(pigeons_stan_dir,"$model.stan"), String)
    end
    model_class = first(split(model,"_"))
    if model_class in ("banana","funnel") 
        return read(joinpath(pigeons_stan_dir,"$model_class.stan"), String)
    end
    error("model_string: model $model unknown")
end

function stan_data(model::String; dataset=nothing, dim=nothing, scale=nothing) 
    if model in ("funnel", "banana")
        Dict("dim" => dim-1, "scale" => scale)
    elseif model in ("funnel_scale", "banana_scale") 
        Dict("dim" => 1, "scale" => scale) # actual dimension is dim+1
    elseif model == "normal"
        Dict("dim" => dim) 
    elseif model == "two_component_normal_scale"
        s_hi, s_lo = two_component_normal_stdevs(scale)
        Dict("n" => 1, "s_hi" => s_hi, "s_lo" => s_lo)
    elseif model == "horseshoe"
        x,y = isnothing(dim) ? make_HSP_data(dataset) : make_HSP_data(dataset,dim) # interpret dim as n_obs
        Dict("n" => length(y), "d" => size(x,2), "x" => x, "y" => y)
    elseif startswith(model,"eight_schools")
        Dict("J" => 8, "y" => [28, 8, -3, 7, -1, 1, 18, 12],
        "sigma" => [15, 10, 16, 11, 9, 11, 10, 18])
    elseif model == "mRNA"
        dta = DataFrame(CSV.File(joinpath(base_dir(), "data", "transfection.csv")))
        Dict("N" => nrow(dta), "ts" => dta[:,1], "ys" => dta[:,3])
    else
        error("stan_data: unknown model $model") 
    end 
end 

# Two component normal for testing preconditioner
two_component_normal_stdevs(e::Real) = (10. ^(e), 10. ^(-e))
function make_2_comp_norm_target(n, exponent)
    s_hi, s_lo = two_component_normal_stdevs(exponent)
    json_str = Pigeons.json(; n=n, s_hi=s_hi, s_lo=s_lo)
    StanLogPotential(joinpath(
        base_dir(), "stan", "two_component_normal.stan"
    ), json_str)
end

# build the horseshoe prior target with varying number of observations
load_HSP_df(dataset::String) = 
    DataFrame(CSV.File(
        joinpath(base_dir(), "data", dataset * ".csv") ))
make_HSP_data(dataset::String, n_obs::Int=typemax(Int)) = 
    make_HSP_data(dataset,load_HSP_df(dataset),n_obs)
function make_HSP_data(dataset::String, df::DataFrame, n_obs::Int)
    iszero(n_obs) && return (zeros( ( n_obs,size(df,2)-1 ) ), Int64[])
    n = min(n_obs, size(df, 1))
    if startswith(dataset,"prostate")
        x = Matrix(df[1:n,2:end])
        y = df[1:n,1]
    elseif startswith(dataset,"ionosphere")
        x = Matrix(hcat(df[1:n,1], df[1:n,3:(end-1)])) # col 2 is constant
        y = Int.(df[1:n,end] .== "g")
    elseif startswith(dataset,"sonar")
        x = Matrix(df[1:n,1:(end-1)])
        y = Int.(df[1:n,end] .== "Mine")
    end
    x,y
end
function make_HSP_target(dataset::String, n_obs::Int=typemax(Int))
    xmat,y = make_HSP_data(dataset, n_obs)
    d = size(xmat,2)
    json_str = if iszero(n_obs)
        Pigeons.json(; n=n_obs,d=d,x="[[]]",y="[]")
    else
        x = [copy(r) for r in eachrow(xmat)]
        Pigeons.json(; n=length(y), d=d, x=x, y=y)
    end
    is_logit = any(Base.Fix1(startswith,dataset), ("prostate", "ionosphere", "sonar"))
    StanLogPotential(joinpath(
        base_dir(), "stan", "horseshoe_" * (is_logit ? "logit" : "linear") * ".stan"
    ), json_str)
end

