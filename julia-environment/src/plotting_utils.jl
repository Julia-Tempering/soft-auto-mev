using AlgebraOfGraphics, CairoMakie, StatsPlots

#include("utils.jl")

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
    hasproperty(df,:miness_per_sec) && return df

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
    df.miness_per_sec = df.miness ./ df.time
    df.miness_per_step = df.miness ./ df.n_steps
    df.miness_per_cost = ifelse.(
        map(sampler -> sampler in ["autoRWMH", "SliceSampler", "HitAndRunSlicer"], df.sampler_type),
        df.miness ./ df.n_logprob,  # non gradient-based samplers
        df.miness ./ (df.n_logprob .+ 2 * df.n_steps .* log_prob_gradient_ratio.(df.model)) # 1 leapfrog = 2 gradient eval
    ) # gradient based: we use cost = #log_potential_eval + eta * #gradient_eval, where eta is model dependent
    df.sampler = map(zip(df.sampler_type, df.selector, df.logstep_jitter)) do (t,s,j)
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
        @df df StatsPlots.groupedboxplot(:model, :miness_per_cost, group=:sampler, ylabel="minESS / cost", 
            legend=:bottomright, xrotation = 20, yaxis=:log10, color=:auto) #title="minESS per cost Comparison for $(my_sampler)", xlabel="Model", 
        savefig(joinpath(plots_path,"$(my_sampler)_miness_cost_comparison.png"))
        @df df StatsPlots.groupedboxplot(:model, :miness_per_sec, group=:sampler, ylabel="minESS / sec", 
            legend=:bottomright, xrotation = 20, yaxis=:log10, color=:auto) #title="minESS per second Comparison for $(my_sampler)", xlabel="Model", 
        savefig(joinpath(plots_path,"$(my_sampler)_miness_sec_comparison.png"))
    end
end

#=
normal jitter comparison plots
=#
function jitter_comparison_plots_model(experiment::String)
    base_df = prepare_df(get_summary_df(experiment))
    base_df = filter(row -> !(row.logstep_jitter == "auto" && row.n_rounds != 21), base_df)
    base_df = filter(row -> !(row.model in ["mRNA", "earning"]), base_df)
    plots_path = joinpath(base_dir(), "deliverables", experiment)

    unique_samplers = unique(base_df[:, :sampler_type])

    # Iterate over each sampler
    foreach(unique_samplers) do my_sampler
        # Filter the dataframe for the current sampler
        df = filter(:sampler_type => ==(my_sampler), base_df)
        sort!(df, :model) # ensure ordering on x-axis
        # Create the grouped boxplot
        @df df StatsPlots.groupedboxplot(:model, :miness_per_cost, group=:logstep_jitter, xlabel="Model", ylabel="minESS / cost", 
            legend=:bottomleft, yaxis=:log10, color=:auto) #title="minESS per cost by Jitter for $(my_sampler)", 
        savefig(joinpath(plots_path,"$(my_sampler)_miness_cost_jitter_comparison.png"))
        @df df StatsPlots.groupedboxplot(:model, :miness_per_sec, group=:logstep_jitter, xlabel="Model", ylabel="minESS / sec", 
            legend=:bottomleft, yaxis=:log10, color=:auto) #title="minESS per second by Jitter for $(my_sampler)", 
        savefig(joinpath(plots_path,"$(my_sampler)_miness_sec_jitter_comparison.png"))
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
        legend=:bottomright, color=:auto, yaxis=:log10) #title="Comparison of minESS per Second for All Samplers", 
    savefig(joinpath(plots_path,"miness_per_sec_comparison.png"))

    # now create the minESS/cost plot
    @df df StatsPlots.groupedboxplot(:model, :miness_per_cost, group=:sampler_type, xlabel="Model", ylabel="minESS / cost", 
        legend=:bottomright, color=:auto, yaxis=:log10) #title="Comparison of minESS per Cost for All Samplers", 
    savefig(joinpath(plots_path,"miness_per_cost_comparison.png"))

    # now the leapfrog comparison plot
    df = filter(row -> !(row.sampler_type in ["autoRWMH", "HitAndRunSlicer"]), df)
    df.leapfrog_per_ess = df.n_steps ./ df.miness
    @df df StatsPlots.groupedboxplot(:model, :leapfrog_per_ess, group=:sampler_type, xlabel="Model", ylabel="Number of Leapfrogs per minESS", 
        legend=:topleft, color=:auto, yaxis=:log10) #title="Comparison of minESS per Cost for All Samplers", 
    savefig(joinpath(plots_path,"num_leapfrog_comparison.png"))
end