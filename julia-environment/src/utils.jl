using CSV 
using DataFrames
using Dates
using DelimitedFiles
using Distributions
using MCMCChains
using BridgeStan
using Random
using SplittableRandoms: SplittableRandom
using Statistics
using StanSample
using JSON
using Turing
using ReverseDiff

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
    elseif startswith(model,"eight_schools_noncentered")
        (1, 0.316857, 0.988146)
    elseif startswith(model, "garch11")
        (1, 5.05101, 0.12328)
    elseif startswith(model, "gp_pois_regr") 
        (1, 5.68368, 0.708641)
    elseif startswith(model, "lotka_volterra")
        (1, 0.547547, 0.0636274)
    elseif startswith(model, "kilpisjarvi") 
        (1, -58.2412, 30.3941)
    elseif startswith(model, "logearn_logheight_male")
        (1, 3.64354, 2.71194)
    elseif startswith(model, "diamonds")
        (1, 6.72838, 0.22442)
    elseif model == "mRNA"
        (3, -1.8909, 1.0014) # 3 => log(beta), hardest to sample together with delta
    elseif startswith(model, "horseshoe")
        (1, 0.88484, 0.116685)
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
    jitter_dist = if jitter_str == "none"
        Dirac(0.0)
    else
        Normal(0.0, jitter_str in ["adapt", "fixed"] ? 0.5 : parse(Float64, jitter_str))
    end

    if sampler_str in ("SliceSampler", "NUTS")                        # irrelevant for NUTS since we use cmdstan
	SliceSampler(n_passes=1)
    elseif sampler_str == "HitAndRunSlicer"
        HitAndRunSlicer(n_refresh=1)
    elseif sampler_str == "SimpleAHMC"                                # we handle autoMALA as special case of SimpleAHMC
	selector = if selector_str == "inverted"
	        autoHMC.AMSelectorInverted()
        elseif selector_str == "non_adaptive"
            autoHMC.AMNonAdaptiveSelector()
        else 
            autoHMC.AMSelectorLegacy() # legacy matches the Pigeons.AutoMALA behavior
        end
	int_time = if int_time_str == "single_step"
            autoHMC.FixedIntegrationTime()                            # Pigeons.AutoMALA is recovered because additionally jitter is Dirac and selector is legacy, see autoHMC tests
        elseif int_time_str == "fixed"
            autoHMC.AdaptiveFixedIntegrationTime()
        else
            autoHMC.AdaptiveRandomIntegrationTime()
        end
        step_jitter = selector_str == "non_adaptive" ? # do not use step jitter if we want to recover HMC or MALA
            autoHMC.StepJitter(dist = Dirac(0.0), adapt_strategy = autoHMC.FixedStepJitter()) : autoHMC.StepJitter(
            dist = jitter_dist,
            adapt_strategy = jitter_str == "adapt" ? autoHMC.AdaptativeStepJitter() : autoHMC.FixedStepJitter()
        )
		SimpleAHMC(
			n_refresh=1, int_time = int_time, step_size_selector = selector, step_jitter = step_jitter
        )
	elseif sampler_str == "SimpleRWMH"
	    selector = if selector_str == "inverted"
	            autoRWMH.MHSelectorInverted()
            elseif sampler_str == "non_adaptive"
                autoRWMH.MHNonAdaptiveSelector()
            else 
                autoRWMH.MHSelector()
            end
        step_jitter = selector_str == "non_adaptive" ? # do not use step jitter if we want to recover RWMH
            autoRWMH.StepJitter(dist = Dirac(0.0), adapt_strategy = autoRWMH.FixedStepJitter()) : autoRWMH.StepJitter(
            dist = jitter_dist,
            adapt_strategy = jitter_str == "adapt" ? autoRWMH.AdaptativeStepJitter() : autoRWMH.FixedStepJitter()
        )
		SimpleRWMH(n_refresh=1, step_size_selector = selector, step_jitter = step_jitter)
	else
		throw(ArgumentError("unknown sampler $sampler_str"))
	end
end

#######################################
# pigeons
#######################################

# recorder for proposal log difference
# It is a GroupBy(Int, Mean()) so can use as template any other such recorder
explorer_proposal_log_diff() = Pigeons.explorer_acceptance_pr()

function pt_sample_from_model(model, target, seed, explorer, miness_threshold; max_rounds = 25)
    n_rounds = 1 # NB: cannot start from >1 otherwise we miss the explorer n_steps from all but last round
    recorders = [
        record_default(); explorer_proposal_log_diff; Pigeons.explorer_acceptance_pr; Pigeons.traces;
        Pigeons.reversibility_rate; online
    ]
    pt = PT(Inputs(
        target      = target, 
        seed        = seed,
        n_rounds    = n_rounds,
        n_chains    = 1, 
        record      = recorders,
        explorer    = explorer, 
        show_report = true
    ))

    # run until minESS threshold is breached
    n_logprob = n_steps = n_samples = 0
    miness = 0.0
    local samples
    while n_rounds ≤ max_rounds # bail after this point
        pt = pigeons(pt)
        n_steps += first(Pigeons.explorer_n_steps(pt))
        samples = get_sample(pt) # only from last round
        n_samples = length(samples)
        n_logprob += if explorer isa SimpleRWMH || explorer isa SliceSampler || explorer isa HitAndRunSlicer
            n_steps # n_steps record the log potential evaluation for non-gradient based samplers
            else
                first(Pigeons.recorder_values(pt, :explorer_n_logprob))
            end
        miness = n_samples < miness_threshold ? 0.0 : min_ess_all_methods(samples, model)
        miness > miness_threshold && break
        @info """
            Low ESS after round $n_rounds: miness = $miness < $miness_threshold = miness_threshold
            Running another round.
        """
        pt = Pigeons.increment_n_rounds!(pt, 1)
        n_rounds += 1 
    end
    mean_1st_dim = first(mean(pt))
    var_1st_dim = first(var(pt))
    step_size = if explorer isa SliceSampler
        pt.shared.explorer.w
    elseif explorer isa HitAndRunSlicer
        pt.shared.explorer.slicer.w
    else
        pt.shared.explorer.step_size
    end
    energy_jump_distance = if pt.shared.explorer isa SimpleAHMC
        autoHMC.energy_jump_distance
    else
        autoRWMH.energy_jump_distance
    end
    energy_jump_dist = first(Pigeons.recorder_values(pt, :energy_jump_distance))
    time = sum(pt.shared.reports.summary.last_round_max_time) # despite name, it is a vector of time elapsed for all rounds
    acceptance_prob = explorer isa SliceSampler ? zero(miness) : 
        first(Pigeons.recorder_values(pt, :explorer_acceptance_pr))
    jitter_std = isa(pt.shared.explorer, HitAndRunSlicer) ? 0 :
        isa(pt.shared.explorer.step_jitter.dist, Normal) ? std(pt.shared.explorer.step_jitter.dist) : 0
    stats_df = DataFrame(
        mean_1st_dim = mean_1st_dim, var_1st_dim = var_1st_dim, time=time, jitter_std = jitter_std, n_logprob = n_logprob, 
        n_steps=n_steps, miness=miness, acceptance_prob=acceptance_prob, step_size=step_size, n_rounds = n_rounds,
        energy_jump_dist = energy_jump_dist)
    return samples, stats_df
end


# remove miness_threshold, and only run a fixed round
function pt_sample_from_model_fixed(model, target, seed, explorer, num_rounds)
    n_rounds = 1 # NB: cannot start from >1 otherwise we miss the explorer n_steps from all but last round
    recorders = [
        record_default(); explorer_proposal_log_diff; Pigeons.explorer_acceptance_pr; Pigeons.traces;
        Pigeons.reversibility_rate; online
    ]
    pt = PT(Inputs(
        target      = target, 
        seed        = seed,
        n_rounds    = num_rounds,
        n_chains    = 1, 
        record      = recorders,
        explorer    = explorer, 
        show_report = true
    ))

    pt = pigeons(pt)
    n_steps = first(Pigeons.explorer_n_steps(pt))
    samples = get_sample(pt) # only from last round
    n_samples = length(samples)
    n_logprob = if explorer isa SimpleRWMH || explorer isa SliceSampler || explorer isa HitAndRunSlicer
        n_steps # n_steps record the log potential evaluation for non-gradient based samplers
        else
            first(Pigeons.recorder_values(pt, :explorer_n_logprob))
        end
    miness = min_ess_all_methods(samples, model)
    mean_1st_dim = first(mean(pt))
    var_1st_dim = first(var(pt))
    step_size = if explorer isa SliceSampler
        pt.shared.explorer.w
    elseif explorer isa HitAndRunSlicer
        pt.shared.explorer.slicer.w
    else
        pt.shared.explorer.step_size
    end
    energy_jump_distance = if pt.shared.explorer isa SimpleAHMC
        autoHMC.energy_jump_distance
    else
        autoRWMH.energy_jump_distance
    end
    energy_jump_dist = first(Pigeons.recorder_values(pt, :energy_jump_distance))
    time = sum(pt.shared.reports.summary.last_round_max_time) # despite name, it is a vector of time elapsed for all rounds
    acceptance_prob = explorer isa SliceSampler ? zero(miness) : 
        first(Pigeons.recorder_values(pt, :explorer_acceptance_pr))
    jitter_std = isa(pt.shared.explorer, HitAndRunSlicer) ? 0 :
        isa(pt.shared.explorer.step_jitter.dist, Normal) ? std(pt.shared.explorer.step_jitter.dist) : 0
    stats_df = DataFrame(
        mean_1st_dim = mean_1st_dim, var_1st_dim = var_1st_dim, time=time, jitter_std = jitter_std, n_logprob = n_logprob, 
        n_steps=n_steps, miness=miness, acceptance_prob=acceptance_prob, step_size=step_size, n_rounds = n_rounds,
        energy_jump_dist = energy_jump_dist)
    return samples, stats_df
end


# record the jitter distribution for each round for the adaptive_jitter explorer
function pt_sample_from_model_round_by_round(model, target, seed, explorer, miness_threshold; max_rounds = 21)
    n_rounds = 1 # NB: cannot start from >1 otherwise we miss the explorer n_steps from all but last round
    recorders = [
        record_default(); explorer_proposal_log_diff; Pigeons.explorer_acceptance_pr; Pigeons.traces;
        Pigeons.reversibility_rate; online
    ]
    pt = PT(Inputs(
        target      = target, 
        seed        = seed,
        n_rounds    = n_rounds,
        n_chains    = 1, 
        record      = recorders,
        explorer    = explorer, 
        show_report = true
    ))
    stats_df = DataFrame(
        mean_1st_dim = [], var_1st_dim = [], time = [], jitter_std = [], n_steps = [], n_logprob =[], 
        miness = [], acceptance_prob=[], step_size=[], n_rounds = [], energy_jump_dist = [])

    # run until max_rounds
    n_logprob = n_steps = n_samples = 0
    miness = 0.0
    local samples
    while n_rounds ≤ max_rounds # bail after this point
        pt = pigeons(pt)
        n_steps += first(Pigeons.explorer_n_steps(pt))
        energy_jump_dist = first(Pigeons.recorder_values(pt, :energy_jump_distance))
        n_logprob += if explorer isa SimpleRWMH || explorer isa SliceSampler || explorer isa HitAndRunSlicer
            n_steps # n_steps record the log potential evaluation for non-gradient based samplers
            else
                first(Pigeons.recorder_values(pt, :explorer_n_logprob))
            end
        samples = get_sample(pt) # only from last round
        n_samples = length(samples)
        miness = n_samples < miness_threshold ? 0.0 : min_ess_all_methods(samples, model)
        @info """
            $n_rounds: miness = $miness
        """
        pt = Pigeons.increment_n_rounds!(pt, 1)
        n_rounds += 1 
        mean_1st_dim = first(mean(pt))
        var_1st_dim = first(var(pt))
        jitter_std = std(pt.shared.explorer.step_jitter.dist) #record the std of jitter distribution (has to be normal)
        step_size = if explorer isa SliceSampler
            pt.shared.explorer.w
        elseif explorer isa HitAndRunSlicer
            pt.shared.explorer.slicer.w
        else
            pt.shared.explorer.step_size
        end
        energy_jump_distance = if pt.shared.explorer isa SimpleAHMC
            autoHMC.energy_jump_distance
        else
            autoRWMH.energy_jump_distance
        end
        energy_jump_dist = first(Pigeons.recorder_values(pt, :energy_jump_distance))
        time = sum(pt.shared.reports.summary.last_round_max_time) # despite name, it is a vector of time elapsed for all rounds
        acceptance_prob = explorer isa SliceSampler ? zero(miness) : 
            first(Pigeons.recorder_values(pt, :explorer_acceptance_pr))
        push!(stats_df, (mean_1st_dim, var_1st_dim, time, jitter_std, n_logprob, 
                n_steps, miness, acceptance_prob, step_size, n_rounds, energy_jump_dist))
    end

    return samples, stats_df
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
    local info
    while n_samples < max_samples
        cmd = stan_cmd(sm, n_samples, seed)
        time += @elapsed run(cmd)
        info = DataFrame(CSV.File(joinpath(sm.tmpdir, model * "_chain_1.csv"), comment = "#"))
        n_steps += sum(info.n_leapfrog__) # count leapfrogs (including warmup)
        col_idxs = push!(collect(8:size(info, 2)), 1)    # discard 7 aux vars, add model variables, then add back logprob at the end 
        samples = info[(sm.num_warmups+1):end, col_idxs] # discard warmup
        miness = min_ess_all_methods(samples, model)
        miness > miness_threshold && break
        @info """
            Low ESS: miness = $miness < $miness_threshold = miness_threshold
            Doubling n_samples and re-running.
        """
        n_samples *= 2 
    end
    mean_1st_dim = mean(samples[:, 1])
    var_1st_dim = var(samples[:, 1])
    stats_df = DataFrame(
        mean_1st_dim = mean_1st_dim, var_1st_dim = var_1st_dim, time=time, 
        n_steps=n_steps, miness=miness, acceptance_prob=zero(miness), step_size=info[end, 3])
    return samples, stats_df
end

function stan_cmd(sm::SampleModel, n_samples, stan_seed; print_every=max(100, round(Int,n_samples/50)))
    cmd = `$(sm.output_base) num_threads=1 sample num_chains=1 num_samples=$n_samples`
    cmd = `$cmd num_warmup=$n_samples save_warmup=true adapt engaged=true`
    cmd = `$cmd id=1 data file=$(first(sm.data_file)) random seed=$stan_seed`
    cmd = `$cmd output file=$(sm.output_base)_chain_1.csv refresh=$print_every`
    return cmd
end

# using NUTS in Turing.jl
function turing_nuts_sample_from_model(model, seed, miness_threshold; max_samples = 2^25, kwargs...)
    # make model and data from the arguments
    my_data = stan_data(model; kwargs...)
    my_model = turing_nuts_model(model, my_data)
    Random.seed!(seed)

    # run until minESS threshold is reached
    n_samples = ceil(Int, 640*miness_threshold) # start from the 6th round to avoid too many rounds
    n_logprob = n_steps = 0
    miness = my_time = 0.0
    local samples
    local chain
    while n_samples < max_samples
        if startswith(model, "horseshoe") # reversediff complains about bernoullilogit
            my_time += @elapsed chain = sample(my_model, NUTS(max_depth=5), n_samples)
        else
            my_time += @elapsed chain = sample(my_model, NUTS(max_depth=5, adtype = AutoReverseDiff()), n_samples)
        end # reduce max_depth because NUTS is super slow
        n_steps += sum(chain[:n_steps]) # count leapfrogs not including warmup
        n_logprob += 2*n_steps + 2*n_samples # Once at the start, 2n during leapfrog, once more at accept/reject
        samples = [chain[param] for param in names(chain)[1:end-12]] # discard 12 aux vars
        samples = [vec(sample) for sample in samples] # convert to vectors
        samples = [collect(row) for row in eachrow(hcat(samples...))] # convert to format compatible with min_ess_all_methods
        miness = min_ess_all_methods(samples, model)
        miness > miness_threshold && break
        @info """
            Low ESS: miness = $miness < $miness_threshold = miness_threshold
            Doubling n_samples and re-running.
        """
        n_samples *= 2 
    end
    mean_1st_dim = mean(samples[1])
    var_1st_dim = var(samples[1])
    acceptance_prob = mean(chain[:acceptance_rate])
    step_size = mean(chain[:step_size])
    energy_jump_dist = mean(abs.(diff(chain[:hamiltonian_energy], dims=1)))
    stats_df = DataFrame(
        mean_1st_dim = mean_1st_dim, var_1st_dim = var_1st_dim, time=my_time, jitter_std = 0, n_logprob = n_logprob, n_steps=n_steps, 
        miness=miness, acceptance_prob=acceptance_prob, step_size=step_size, n_rounds = log2(n_samples/(10*miness_threshold)), energy_jump_dist = energy_jump_dist)
    return samples, stats_df
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
    elseif model_name == "funnel"
        2.0
    elseif model_name == "banana"
        1.0
    elseif model_name == "normal"
        1.0
    else
		throw(ArgumentError("Don't know how to make a scale argument for $model_name"))
	end
end

# returns a string with the stan code for a model
function model_string(model; dataset=nothing, kwargs...)
    # first: try models in the local stan folder
    file_name = if startswith(model, "horseshoe")
        is_logit = any(Base.Fix1(startswith,dataset), ("prostate", "ionosphere", "sonar"))
        "horseshoe_" * (is_logit ? "logit" : "linear")
    else
        model
    end
    model_str = try
        read(joinpath(base_dir(), "stan", file_name * ".stan"), String)
    catch e
        e isa SystemError || rethrow(e)
        ""
    end
    isempty(model_str) || return model_str
    
    # now see if we find it in Pigeons' examples dir
    pigeons_stan_dir = joinpath(pkgdir(Pigeons), "examples", "stan")
    model_class = first(split(model,"_"))
    file_name = if model_class in ("banana","funnel") 
        model_class
    else
        model
    end
    model_str = read(joinpath(pigeons_stan_dir, file_name * ".stan"), String)
    return model_str
end

function stan_data(model::String; dataset=nothing, dim=nothing, scale=nothing) 
    if (startswith(model, "funnel") || startswith(model, "banana"))
        Dict("dim" => dim-1, "scale" => scale)
    elseif model in ("funnel_scale", "banana_scale") 
        Dict("dim" => 1, "scale" => scale) # actual dimension is dim+1
    elseif model == "normal"
        Dict("dim" => dim) 
    elseif model == "two_component_normal_scale"
        s_hi, s_lo = two_component_normal_stdevs(scale)
        Dict("n" => 1, "s_hi" => s_hi, "s_lo" => s_lo)
    #elseif startswith(model, "horseshoe")
        #x,y = isnothing(dim) ? make_HSP_data(dataset) : make_HSP_data(dataset,dim) # interpret dim as n_obs
        #Dict("n" => length(y), "d" => size(x,2), "x" => x', "y" => y)
    elseif startswith(model,"eight_schools")
        Dict("J" => 8, "y" => [28, 8, -3, 7, -1, 1, 18, 12],
        "sigma" => [15, 10, 16, 11, 9, 11, 10, 18])
    #elseif model == "mRNA"
        #Dict(pairs(mrna_load_data()))
    else
        file_name = if startswith(model, "earn_height")
            "earnings"
        elseif startswith(model,"nes")
            "nes2000"
        elseif startswith(model,"diamonds")
            "diamonds"
        elseif startswith(model,"hmm_example")
            "hmm_example"
        elseif startswith(model, "eight_schools_noncentered")
            "eight_schools"
        elseif startswith(model,"garch11")
            "garch"
        elseif startswith(model,"gp_pois_regr")
            "gp_pois_regr"
        elseif startswith(model,"lotka_volterra")
            "hudson_lynx_hare"
        elseif startswith(model, "kilpisjarvi")
            "kilpisjarvi_mod"
        elseif startswith(model,"logearn_logheight_male")
            "earnings"
        elseif startswith(model, "horseshoe")
            "sonar"
        elseif startswith(model, "mRNA")
            "mRNA"
        else
            error("unknown model $model")
        end
        JSON.parse(read(joinpath(base_dir(), "data", file_name * ".json"), String))
    end 
end

# utility for creating StanLogPotentials for real data models
# (or basically any model that does not have already a StanLogPotential defined in Pigeons)
# IDEA: reuse the machinery for cmdstan, to minimize divergence
# also forces use of temp file; this avoids race conditions when
# multiple nodes are compiling a file in a shared location (like base_dir()/stan)
function stan_logpotential(model; dataset = nothing, dim = nothing, scale = nothing)
    tmpdir = mktempdir()
    isdir(tmpdir) || mkdir(tmpdir)
    
    # write .stan file
    model_str = model_string(model; dataset = dataset)
    stan_fname = joinpath(tmpdir, model * ".stan")
    open(stan_fname, "w") do f
        write(f, model_str)
    end

    # write .json file    
    data_str = json(stan_data(model; dataset = dataset, dim = dim, scale = scale))
    data_fname = joinpath(tmpdir, model * ".json")
    open(data_fname, "w") do f
        write(f, data_str)
    end

    # build the log potential
    StanLogPotential(stan_fname, data_fname)
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

function mrna_load_data()
    dta_path = joinpath(pkgdir(Pigeons), "examples", "data", "Ballnus_et_al_2017_M1a.csv")
    dta = readdlm(dta_path, ',')
    N   = size(dta, 1)
    ts  = dta[:,1]
    ys  = dta[:,2]
    return (;N = N, ts = ts, ys = ys)
end

# utility function to create models for Turing.NUTS()
function turing_nuts_model(model, data)
    if startswith(model, "diamonds")
        @model function diamonds(Y, X, N, K)
            # Data dimensions
            Kc = K - 1
            # Centered version of X without intercept
            Xc = Matrix{Real}(undef, N, Kc)
            means_X = Vector{Real}(undef, Kc)
            for i in 2:K
                means_X[i - 1] = mean(X[:, i])
                Xc[:, i - 1] = X[:, i] .- means_X[i - 1]
            end
            # Parameters
            b = Vector{Real}(undef, Kc)
            Intercept ~ LocationScale(8, 10, TDist(3))  # Intercept
            b .~ Normal(0, 1)  # Population-level effects
            sigma ~ truncated(LocationScale(0, 10, TDist(3)), 0, Inf)  # Residual SD
            # Likelihood
            for i in 1:N
                    Y[i] ~ Normal(sum(Xc[i, :].*b) + Intercept, sigma)
            end
        end
        return diamonds(data["Y"], hcat(data["X"]...)', data["N"], data["K"])
    elseif startswith(model, "horseshoe")
        @model function horseshoe(x, y, n, d)
            # Priors
            tau ~ TDist(1)  # Cauchy(0, 1) equivalent to TDist with 1 degree of freedom
            lambda ~ filldist(TDist(1), d)  # d-dimensional vector with Cauchy(0, 1)
            beta0 ~ LocationScale(0, 1, TDist(3))  # Intercept with StudentT(3, 0, 1)
            beta ~ MvNormal(zeros(d), tau * lambda)  # Multivariate normal with scale tau * lambda
            
            # Likelihood
            for i in 1:n
                y[i] ~ BernoulliLogit(beta0 + sum(x[i, :] .* beta))  # Logistic regression likelihood
            end
        end
        return horseshoe(hcat(data["x"]...)', data["y"], data["n"], data["d"])
    elseif startswith(model, "kilpisjarvi")
        @model function kilpisjarvi(x, y, N, xpred, pmualpha, psalpha, pmubeta, psbeta)
            # Priors
            alpha ~ Normal(pmualpha, psalpha)  # Prior for alpha
            beta ~ Normal(pmubeta, psbeta)     # Prior for beta
            sigma ~ truncated(Normal(0, 1), 0, Inf)  # Prior for sigma, ensuring it's positive
        
            # Likelihood
            for i in 1:N
                y[i] ~ Normal(alpha + beta * x[i], sigma)  # Observations model
            end
        
            # Predicted value at xpred (optional)
            ypred = alpha + beta * xpred
        end
        return kilpisjarvi(data["x"], data["y"], data["N"], data["xpred"], data["pmualpha"], 
                data["psalpha"], data["pmubeta"], data["psbeta"])
    elseif startswith(model, "logearn_logheight_male")
        @model function logearn(N, earn, height, male)
            # Transformed data: log transformations
            log_earn = log.(earn)
            log_height = log.(height)
            # Parameters
            beta ~ MvNormal(zeros(3), ones(3))  # Prior for beta, assuming standard normal priors
            sigma ~ truncated(Normal(0, 1), 0, Inf)  # Positive constraint for sigma
            # Likelihood
            for i in 1:N
                log_earn[i] ~ Normal(beta[1] + beta[2] * log_height[i] + beta[3] * male[i], sigma)
            end
        end
        return logearn(data["N"], data["earn"], data["height"], data["male"])
    elseif startswith(model, "mRNA")
        function exp_a_minus_exp_b(a, b)
            return a > b ? -exp(a) * expm1(b - a) : exp(b) * expm1(a - b)
        end
        function get_mu(tmt0, km0, beta, delta)
            if tmt0 <= 0.0
                return 0.0  # must force mu=0 when t < t0 (reaction hasn't started)
            end
            dmb = delta - beta
            if abs(dmb) < eps()  # `eps()` gives machine precision in Julia
                return km0 * tmt0
            else
                return km0 * exp_a_minus_exp_b(-beta * tmt0, -delta * tmt0) / dmb
            end
        end
        @model function mRNA_model(N, ts, ys)
            # Parameters with transformations (log-transformed priors)
            lt0 ~ Uniform(-2, 1)
            lkm0 ~ Uniform(-5, 5)
            lbeta ~ Uniform(-5, 5)
            ldelta ~ Uniform(-5, 5)
            lsigma ~ Uniform(-2, 2)
            # Transformed parameters
            t0 = 10.0^lt0
            km0 = 10.0^lkm0
            beta = 10.0^lbeta
            delta = 10.0^ldelta
            sigma = 10.0^lsigma
            # Likelihood
            for i in 1:N
                mu_i = get_mu(ts[i] - t0, km0, beta, delta)
                ys[i] ~ Normal(mu_i, sigma)
            end
        end
        return mRNA_model(data["N"], data["ts"], data["ys"])
    elseif startswith(model, "funnel")
        @model function funnel(dim, scale)
            x = Vector{Real}(undef, dim + 1)
            x[1] ~ Normal(0, 3)
            for i in 2:dim
                x[i] ~ Normal(0, exp(x[1] / scale)) 
            end
        end
        return funnel(data["dim"], data["scale"])
    elseif startswith(model, "banana")
        @model function banana(dim, scale)
            x = Vector{Real}(undef, dim + 1)
            x[1] ~ Normal(0, sqrt(10))
            for i in 2:dim
                x[i] ~ Normal(x[1]^2, scale / sqrt(10)) 
            end
        end
        return banana(data["dim"], data["scale"])
    else
        error("unknown model $model")
    end
end


