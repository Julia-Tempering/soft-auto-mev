using Pkg
Pkg.activate("julia-environment") 
include(joinpath("../julia-environment/src/utils.jl")) 
using CairoMakie
using StatsPlots

function make_explorer(explorer)
    if explorer == "HitAndRunSlicer"
        HitAndRunSlicer()
    elseif explorer == "autoRWMH"
        SimpleRWMH(step_jitter = autoRWMH.StepJitter(dist = Normal(0.0, 3.0),adapt_strategy = autoRWMH.FixedStepJitter()))
    elseif explorer == "autoMALA"
        SimpleAHMC(n_refresh=1, int_time = autoHMC.FixedIntegrationTime())
    else
        SimpleAHMC(n_refresh=1, int_time = autoHMC.AdaptiveRandomIntegrationTime())
    end
end

function pairplot(df; title="Pairplot")
    rows, cols = size(df)
    plots = []
	p_mat = Matrix(undef,cols,cols);
    for i = 1:cols, j = 1:cols
		p_mat[i,j] = i>=j ? (label = :°, blank = false) : (label = :_, blank = true)
        if i>j
            subplot = StatsPlots.scatter(df[:,i], df[:,j], legend = false,markersize=1/cols, 
			alpha=0.1, markerstrokecolor=nothing, grid=nothing, markerstrokealpha=0.0, 
			tickfontsize=round(32/cols), tick_direction=:in, 
			xrot=45,showaxis=(i==cols ? (j==1 ? true : :x) : (j==1 ? :y : false)),
			xticks=(i==cols ? :auto : nothing), yticks=(j==1 ? :auto : nothing))
			push!(plots,subplot)
		elseif i==j
			subplot = StatsPlots.histogram(df[:,i], normalize=true, legend=false, alpha=0.4, ticks=nothing, showaxis=(i==cols ? (j==1 ? true : :x) : (j==1 ? :y : false)),xticks=(i==cols ? :auto : nothing), 
				tick_direction=:in, tickfontcolor="white", xrot=45
			)
			StatsPlots.density!(df[:,i], grid=nothing, tickfontsize=round(32/cols), trim=true,
			xlimits=(minimum(df[:,i]),maximum(df[:,i])))
			push!(plots,subplot)
		end
    end
    return StatsPlots.plot(plots..., layout=p_mat, plot_titlevspan=0.05) #plot_title=title, 
end

function get_samples(explorer)
    if explorer == "NUTS"
        my_data = stan_data("mRNA")
        my_model = turing_nuts_model("mRNA", my_data)
        Random.seed!(1)
        chain = sample(my_model, NUTS(max_depth=5, adtype = AutoReverseDiff()), 2^15)
        samples = [chain[param] for param in names(chain)[1:end-12]]
        df = DataFrame(lt0 = vec(samples[1][:]), lkm0 = vec(samples[2][:]), lbeta = vec(samples[3][:]), ldelta = vec(samples[4][:]), lsigma = vec(samples[5][:]))
    else
        pt = pigeons(
        target     = stan_logpotential("mRNA"),
        seed       = 2,
        n_rounds   = 18,
        n_chains   = 1, 
        record     = [record_default(); Pigeons.traces; online],
        explorer   = make_explorer(explorer),
        show_report = true
        )
        samples = vcat(get_sample(pt))
        # Extract x coordinates
        x1 = [sample[1] for sample in samples]
        x2 = [sample[2] for sample in samples]
        x3 = [sample[3] for sample in samples]
        x4 = [sample[4] for sample in samples]
        x5 = [sample[5] for sample in samples]
        df = DataFrame(lt0 = x1, lkm0 = x2, lbeta = x3, ldelta = x4, lsigma = x5)
    end
end

# draw pairplot
for explorer in ["autoRWMH"]#, "autoHMC", "autoMALA", "NUTS", "HitAndRunSlicer"]
    df = get_samples(explorer)
    fig = pairplot(df)
    # fig = StatsPlots.cornerplot(df, compact = true)
    # Display the figure with both trace plots side by side
    # display(fig)
    save("deliverables/pairplot/pairplot_$explorer.png", fig)
end