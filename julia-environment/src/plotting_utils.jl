# run this once when opening a new terminal
using Pkg
Pkg.activate("julia-environment")
Pkg.update()
include("utils.jl")


using AlgebraOfGraphics, CairoMakie, StatsPlots, CategoricalArrays

function get_summary_df(experiment::String)
    base_folder = base_dir()
    csv_path    = joinpath(base_folder, "deliverables", experiment, "aggregated", "summary.csv")
    return DataFrame(CSV.File(csv_path))
end

# ratio of running time of gradient VS log potential
# computed separately, recording the avg of time_gradient/time_log_prob
function log_prob_gradient_ratio(model::AbstractString)
    if startswith(model, "horseshoe")
        35.66540160529861 # another run: 35.45077959104718
    elseif startswith(model, "mRNA")
        5.766575585884267 # 5.793217015357431
    elseif startswith(model, "earning")
        4.129503888073945 # 4.174167584502551
    elseif startswith(model, "kilpisjarvi")
        5.956119530847897 # 6.0015265040396
    elseif startswith(model, "diamonds")
        3.567778885451556 # 3.5331787529624483
    elseif startswith(model, "funnel(2,10)")
        3.571520466162771 # 3.6129133112418566
    elseif startswith(model, "funnel(128,10)")
        62.33694488486751 # 62.743038514474065
    elseif startswith(model, "funnel(128,1)")
        65.43722654676135 # 63.34395522158833
    elseif startswith(model, "funnel")
        4.04648547985939 # 4.0364094087006075
    elseif startswith(model, "banana(2,10)")
        3.431247145051592 # 3.416732917597834
    elseif startswith(model, "banana(128,10)")
        57.05919054113505 # 58.53447763204794
    elseif startswith(model, "banana(128,1)")
        58.94424841180857 # 59.8881857878015
    elseif startswith(model, "banana")
        3.6526912428378147 # 3.634853668305632
    elseif startswith(model, "normal(2,10)")
        5.301028059622657 # 5.313744330782199
    elseif startswith(model, "normal(128,10)")
        54.62782019934652 # 54.60099473608651
    elseif startswith(model, "normal(128,1)")
        54.14752786431609 # 54.29806592940469
    elseif startswith(model, "normal(2,1)")
        5.673805495244712 # 5.785713132240938
    elseif startswith(model, "normal")
        14.19361594425169 # 14.444763963652496  
    else
        throw(KeyError(model))
    end
end

function prepare_df(df::DataFrame)
	hasproperty(df, :miness_per_sec) && return df

	idxs = @. df.sampler_type == "SimpleAHMC" && df.int_time == "single_step"
	df[idxs, :sampler_type] .= "autoMALA"
	idxs = @. df.sampler_type == "SimpleAHMC"
	df[idxs, :sampler_type] .= "autoHMC"
	idxs = @. df.sampler_type == "SimpleRWMH"
	df[idxs, :sampler_type] .= "autoRWMH"
	idxs = @. df.model == "horseshoe_logit"
	df[idxs, :model] .= "horseshoe"
	idxs = @. df.model == "logearn_logheight_male"
	df[idxs, :model] .= "earning"
	idxs = @. df.logstep_jitter == "adapt"
	df[idxs, :logstep_jitter] .= "auto"
    df.energy_jump_dist = ifelse.(
		map(sampler -> sampler in ["NUTS"], df.sampler_type), # NUTS returns total energy jump dist
		df.energy_jump_dist ./ (8000 * 2 .^ df.n_rounds), # first round has 8000 samples
		df.energy_jump_dist
	)
	df.miness_per_sec = df.miness ./ df.time
	df.miness_per_step = df.miness ./ df.n_steps
	df.cost = ifelse.(
		map(sampler -> sampler in ["RWMH", "RWMH0.5", "RWMH2.0", "autoRWMH", "SliceSampler", "HitAndRunSlicer"], df.sampler_type),
		df.n_logprob,  # non gradient-based samplers
		df.n_logprob .+ 2 * df.n_steps .* log_prob_gradient_ratio.(df.model), # 1 leapfrog = 2 gradient eval
	) # gradient based: we use cost = #log_potential_eval + eta * #gradient_eval, where eta is model dependent
	df.miness_per_cost = df.miness ./ df.cost
	df.sampler = map(zip(df.sampler_type, df.selector, df.logstep_jitter)) do (t, s, j)
		t in ("NUTS", "SliceSampler", "HitAndRunSlicer") ? t : t * (s == "inverted" ? "_symmetric" : "") * (j in ["none"] ? "" : "_jitter")
	end
	return df
end

function scatter_miness_cost(experiment::String)
    df = filter(:scale_idx => ==(1), get_summary_df(experiment))
    df = prepare_df(df)
    plots_path = joinpath(base_dir(), "deliverables", experiment)

    foreach((:miness_per_sec, :miness_per_step)) do sym
        fig_plan = data(df) * 
            visual(Scatter) *
            mapping(
                sym => log10 => "minESS / " * (sym == :miness_per_sec ? "second" : "step") * " (log scale)",
                :sampler => identity => "Sampler",
                color=:sampler_type
            )
        fig = draw(fig_plan, axis = (width = 400, height = 400))
        save(joinpath(plots_path,"$sym.png"), fig, px_per_unit = 3) # save high-resolution png
    end
end


#=
High-dimensional scaling plots
=#

function highdim_plots()
    experiment = "highdim"
    df = get_summary_df(experiment)
    likely_miness_threshold = 10^floor(log10(median(df.miness)))
    df = filter(:miness => >=(likely_miness_threshold), df)
    df = prepare_df(df)
    foreach(("banana","funnel","normal")) do model
        highdim_plots_model(filter(:model => ==(model), df))
    end
end

function axis_labeller(sym, trans_str=" (log₁₀)")
    base_str = if sym == :miness_per_sec
        "minESS / second"
    elseif sym == :miness_per_step
        "minESS / step"
    elseif sym == :step_size
        "Step size"
    end
    base_str * trans_str
end

function highdim_plots_model(model_df)
    model = first(model_df.model)
    max_dim = maximum(model_df.dim)
    max_dim_df = filter(:dim => ==(max_dim), model_df)

    plots_path = joinpath(base_dir(), "deliverables", experiment, "$model")
    foreach((:miness_per_sec, :miness_per_step)) do sym
        f = data(max_dim_df) * 
            visual(BoxPlot) *
            mapping(
                :sampler,
                sym => log10 => axis_labeller(sym),
                color=:sampler_type
            )
        p = draw(f, axis = (width = 1400, height = 400, title="Results for $model at highest dimension available ($max_dim)"))
        save(plots_path * "_max_dim_$sym.png", p, px_per_unit = 3) # save high-resolution png
    end

    # scaling of best combinations
    gdf = groupby(max_dim_df, [:sampler,:sampler_type])
    cgdf = combine(gdf, :miness_per_step => median, renamecols=false)
    best_combinations = combine(
        groupby(cgdf, :sampler_type),
        [:sampler, :miness_per_step] => ((s,m) -> s[argmax(m)]) => :sampler
    )
    select!(best_combinations, Not(:sampler_type))
    model_df_small = innerjoin(model_df, best_combinations; on=:sampler)
    foreach((:miness_per_sec, :miness_per_step, :step_size)) do sym
        f = data(model_df_small) * 
            visual(BoxPlot) *
            mapping(
                :dim => log2 => "Dimension (log₂)",
                sym => log10 => axis_labeller(sym),
                color=:sampler,
                dodge=:sampler
            )
        p=draw(f,axis = (width = 800, height = 400, title="Dimensional scaling for $model, best-performing combinations only"))
        save(plots_path * "_scaling_$sym.png", p, px_per_unit = 3) # save high-resolution png
    end
end

#=
jitter-stability plots
=#

function jitter_stability_plots_model(experiment::String)
    base_df = prepare_df(get_summary_df(experiment))
    plots_path = joinpath(base_dir(), "deliverables", experiment)
    base_df = filter(row -> row.logstep_jitter == "auto", base_df) # only consider stability for adaptive jitter
    base_df = filter(row -> !(row.model in ["kilpisjarvi", "earning"]), base_df)
    base_df.jitter_std .= ifelse.(base_df.jitter_std .<= 0.0, 10^(-5), base_df.jitter_std) # prevent zero jitter std

    unique_samplers = unique(base_df[:, :sampler_type])

    # Iterate over each sampler
    foreach(unique_samplers) do sampler
        # Filter the dataframe for the current sampler
        df = filter(:sampler_type => ==(sampler), base_df)
        # Group the data by 'model' and 'n_rounds', and calculate the mean and SEM for 'jitter_std'
        grouped_df = combine(groupby(df, [:model, :n_rounds]), 
                        :jitter_std => mean => :jitter_std_mean, 
                        :jitter_std => (x -> std(x) / sqrt(length(x))) => :jitter_std_sem)

        # Plot the data
        @df grouped_df StatsPlots.plot(:n_rounds, :jitter_std_mean, group=:model, lw=2, legend=:topright,
            ribbon=:jitter_std_sem, xlabel="tuning rounds", ylabel="jitter standard deviation", 
            #title="Convergence of Jitter Standard Deviation for $sampler",
            yaxis=:log10, markershape=:auto, linecolor=:auto) #, xaxis=:log2
    
        savefig(joinpath(plots_path,"$(sampler)_jitter_convergence.png"))
    end
end

#=
within-sampler comparison plots
=#
function within_sampler_comparison_plots_model(experiment::String)
    base_df = prepare_df(get_summary_df(experiment))
    plots_path = joinpath(base_dir(), "deliverables", experiment)

    unique_samplers = unique(base_df[:, :sampler_type])

    # Iterate over each sampler
    foreach(unique_samplers) do my_sampler
        # Filter the dataframe for the current sampler
        df = filter(:sampler_type => ==(my_sampler), base_df)
        sort!(df, :model) # ensure ordering on x-axis
        # Create the grouped boxplot
        yticks = 10.0 .^ LinRange(-6, 6, 13)
        @df df StatsPlots.groupedboxplot(:model, :miness_per_cost, group=:sampler, ylabel="minESS / cost", 
            legend=:bottomright, xrotation = 20, yaxis=:log10, color=:auto, yticks = yticks)
        savefig(joinpath(plots_path,"$(my_sampler)_miness_cost_comparison.png"))

        yticks = 10.0 .^ LinRange(0, 4, 5)
        @df df StatsPlots.groupedboxplot(:model, :miness_per_sec, group=:sampler, ylabel="minESS / sec", 
            legend=:bottomright, xrotation = 20, yaxis=:log10, color=:auto, yticks = yticks)
        savefig(joinpath(plots_path,"$(my_sampler)_miness_sec_comparison.png"))

        yticks = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
        @df df StatsPlots.groupedboxplot(:model, :energy_jump_dist, group=:sampler, ylabel="Average Energy Jump Distance", 
            legend=:topleft, xrotation = 20, color=:auto, yticks = yticks, ylim = (0.05, 0.3))
        savefig(joinpath(plots_path,"$(my_sampler)_energy_jump_comparison.png"))
    end
end

#=
normal jitter comparison plots
=#
function jitter_comparison_plots_model(experiment::String)
    base_df = prepare_df(get_summary_df(experiment))
    base_df = filter(row -> !(row.logstep_jitter == "auto" && row.n_rounds != 21), base_df)
    base_df = filter(row -> !(row.model in ["mRNA", "earning"]), base_df)
    replace!(base_df.logstep_jitter, "none" => "0.0")
    plots_path = joinpath(base_dir(), "deliverables", experiment)

    unique_samplers = unique(base_df[:, :sampler_type])

    # Iterate over each sampler
    foreach(unique_samplers) do my_sampler
        # Filter the dataframe for the current sampler
        df = filter(:sampler_type => ==(my_sampler), base_df)
        sort!(df, :model) # ensure ordering on x-axis
        # Create the grouped boxplot
        yticks = 10.0 .^ LinRange(-6, 6, 13)
        @df df StatsPlots.groupedboxplot(:model, :miness_per_cost, group=:logstep_jitter, xlabel="Model", ylabel="minESS / cost", 
            legend=:bottomleft, yaxis=:log10, color=:auto, yticks = yticks)
        savefig(joinpath(plots_path,"$(my_sampler)_miness_cost_jitter_comparison.png"))

        yticks = 10.0 .^ LinRange(-4, 4, 9)
        @df df StatsPlots.groupedboxplot(:model, :miness_per_sec, group=:logstep_jitter, xlabel="Model", ylabel="minESS / sec", 
            legend=:bottomleft, yaxis=:log10, color=:auto, yticks = yticks)
        savefig(joinpath(plots_path,"$(my_sampler)_miness_sec_jitter_comparison.png"))

        yticks = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
        @df df StatsPlots.groupedboxplot(:model, :energy_jump_dist, group=:logstep_jitter, xlabel="Model", 
        ylabel="Average Energy Jump Distance", legend=:topleft, color=:auto, yticks = yticks, ylim = (0.05, 0.3))
        savefig(joinpath(plots_path,"$(my_sampler)_energy_jump_jitter_comparison.png"))
    end
end

#=
comparison of all autoMCMC samplers and NUTS; experiment = "post_db"
=#
function all_comparison_plots_model(experiment::String)
    df = prepare_df(get_summary_df(experiment))
    plots_path = joinpath(base_dir(), "deliverables", experiment)

    sort!(df, :model) # ensure ordering on x-axis
    # Create the grouped boxplot for minESS/sec
    @df df StatsPlots.groupedboxplot(:model, :miness_per_sec, group=:sampler_type, xlabel="Model", ylabel="minESS / sec", 
        legend=:bottomleft, color=:auto, yaxis=:log10)
    savefig(joinpath(plots_path,"miness_per_sec_comparison.png"))

    # now create the minESS/cost plot
    @df df StatsPlots.groupedboxplot(:model, :miness_per_cost, group=:sampler_type, xlabel="Model", ylabel="minESS / cost", 
        legend=:bottomleft, color=:auto, yaxis=:log10)
    savefig(joinpath(plots_path,"miness_per_cost_comparison.png"))

    # now create the energy jump plot
    @df df StatsPlots.groupedboxplot(:model, :energy_jump_dist, group=:sampler_type, xlabel="Model", 
    ylabel="Average Energy Jump Distance", legend=:topleft, color=:auto, ylim= (-0.03,0.8))
    savefig(joinpath(plots_path,"energy_jump_comparison.png"))

    # now the leapfrog comparison plot
    #df = filter(row -> !(row.sampler_type in ["autoRWMH", "HitAndRunSlicer"]), df)
    #df.leapfrog_per_ess = df.n_steps ./ df.miness
    #@df df StatsPlots.groupedboxplot(:model, :leapfrog_per_ess, group=:sampler_type, xlabel="Model", ylabel="Number of Leapfrogs per minESS", 
    #    legend=:topleft, color=:auto, yaxis=:log10)
    #savefig(joinpath(plots_path,"num_leapfrog_comparison.png"))
end


#=
comparison of autoMCMC with non autoMCMC; experiment = "post_db", "auto_traditional"
=#
function non_auto_plots_model(experiment::String)
	df = prepare_df(get_summary_df(experiment))
    df = filter(row -> !isnan(row.energy_jump_dist), df)
    replace!(df.model, "banana" => "banana(4,1.0)")
    replace!(df.model, "funnel" => "funnel(4,1.0)")
    replace!(df.sampler_type, "RWMH" => "RWMH1.0")
    replace!(df.sampler_type, "MALA" => "MALA1.0")
    replace!(df.sampler_type, "HMC" => "HMC1.0")
	plots_path = joinpath(base_dir(), "deliverables", experiment)

	sort!(df, :model) # ensure ordering on x-axis
	loop_samplers = [
        CategoricalArray(["RWMH", "RWMH0.1", "RWMH0.25", "RWMH1.0", "RWMH4.0", "RWMH10.0", "autoRWMH"], ordered = true),
		CategoricalArray(["MALA", "MALA0.1", "MALA0.25", "MALA1.0", "MALA4.0", "MALA10.0", "autoMALA"], ordered = true),
		CategoricalArray(["HMC", "HMC0.1", "HMC0.25", "HMC1.0", "HMC4.0", "HMC10.0", "autoHMC"], ordered = true)]

	foreach(loop_samplers) do samplers
		df2 = filter(row -> (row.sampler_type in levels(samplers)), df)
        df2.sampler_type = CategoricalArray(df2.sampler_type, ordered=true, levels=samplers)

        yticks = 10.0 .^ LinRange(-6, 6, 13)
        @df df2 StatsPlots.groupedboxplot(:model, :miness_per_sec, group = :sampler_type, #xlabel = "Model", 
            ylabel = "minESS / sec", legend = :bottomleft, color = :auto, yaxis = :log10, yticks = yticks,
            ylim = (0.9 * minimum(df2.miness_per_sec), 1.3 * maximum(df2.miness_per_sec)),
            xrotation = 15)
        # Overlay the number of points as text
        #annotate_model_points!(unique_models, samplers, df2, maximum(df2.miness_per_sec), model_mapping, sampler_offset)
		savefig(joinpath(plots_path, "$(samplers[1])_miness_per_sec.png"))

		# now create the minESS/cost plot
        yticks = 10.0 .^ LinRange(-2, 2, 5)
		@df df2 StatsPlots.groupedboxplot(:model, :miness_per_cost, group = :sampler_type, #xlabel = "Model", 
        ylabel = "minESS / cost", legend = :bottomleft, color = :auto, yaxis = :log10, yticks = yticks,
        ylim = (0.9 * minimum(df2.miness_per_cost), 1.3 * maximum(df2.miness_per_cost)),
        xrotation = 15)
        # Overlay the number of points as text
        #annotate_model_points!(unique_models, samplers, df2, maximum(df2.miness_per_cost), model_mapping, sampler_offset)
		savefig(joinpath(plots_path, "$(samplers[1])_miness_per_cost_.png"))

		# now create the energy jump distance
        yticks = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
		@df df2 StatsPlots.groupedboxplot(:model, :energy_jump_dist, group = :sampler_type, #xlabel = "Model",
			ylabel = "Average Energy Jump Distance", legend = :topright, color = :auto, yticks = yticks, ylim = (-0.03, 0.3),
            xrotation = 15)
        #annotate_model_points!(unique_models, samplers, df2, maximum(df2.energy_jump_dist), model_mapping, sampler_offset)
		savefig(joinpath(plots_path, "$(samplers[1])_energy_jump_dist.png"))
	end
end


function annotate_model_points!(unique_models, samplers, df2, position, model_mapping, sampler_offset)
    for m in unique_models
        for s in samplers
            # Count the number of points for each model and sampler_type combination
            num_points = count((df2.model .== m) .& (df2.sampler_type .== s))
            
            # Find the jittered model value (x-axis position) for this combination
            jittered_x = model_mapping[m] + sampler_offset[s] - 1
            
            # Add text to show the number of points at the appropriate location
            annotate!(jittered_x, 1.2 * position, StatsPlots.text(string(num_points), 6, :center, :black), 
                      halign = :center, valign = :top, color = :black)
        end
    end
end



#=
acceptance prob plot for autostep and fixed-step methods; experiment = "fixed_round"
=#
function acceptance_prob_plots_model(experiment::String)
    df = prepare_df(get_summary_df(experiment))
    replace!(df.sampler_type, "RWMH" => "RWMH1.0")
    replace!(df.sampler_type, "MALA" => "MALA1.0")
    replace!(df.sampler_type, "HMC" => "HMC1.0")
    plots_path = joinpath(base_dir(), "deliverables", experiment)

    sort!(df, :model) # ensure ordering on x-axis
	loop_samplers = [
        CategoricalArray(["RWMH", "RWMH0.1", "RWMH0.25", "RWMH1.0", "RWMH4.0", "RWMH10.0", "autoRWMH"], ordered = true),
		CategoricalArray(["MALA", "MALA0.1", "MALA0.25", "MALA1.0", "MALA4.0", "MALA10.0", "autoMALA"], ordered = true),
		CategoricalArray(["HMC", "HMC0.1", "HMC0.25", "HMC1.0", "HMC4.0", "HMC10.0", "autoHMC"], ordered = true)]

	foreach(loop_samplers) do samplers
		df2 = filter(row -> (row.sampler_type in samplers), df)
        df2.sampler_type = CategoricalArray(df2.sampler_type, ordered=true, levels=samplers)
        yticks = [0,0.2,0.4,0.6,0.8,1.0]

		@df df2 StatsPlots.groupedboxplot(:model, :acceptance_prob, group = :sampler_type, 
			ylabel = "Average acceptance probability", legend = :topright, color = :auto, yticks = yticks, 
            xrotation = 15)
		savefig(joinpath(plots_path, "$(samplers[1])_acceptance_prob.png"))
    end
end