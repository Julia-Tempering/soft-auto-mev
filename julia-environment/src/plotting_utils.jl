using AlgebraOfGraphics, CairoMakie, StatsPlots

include("utils.jl")

function get_summary_df(experiment::String)
    base_folder = base_dir()
    csv_path    = joinpath(base_folder, "deliverables", experiment, "aggregated", "summary.csv")
    return DataFrame(CSV.File(csv_path))
end

# ratio of running time of gradient VS log potential
# computed separately, recording the avg of time_gradient/time_log_prob
function log_prob_gradient_ratio(model::String)
    if startswith(model, "horseshoe")
        35.66540160529861 # another run: 35.45077959104718
    elseif startswith(model, "mRNA")
        5.766575585884267 # another run: 5.793217015357431
    elseif startswith(model, "logearn_logheight_male")
        4.129503888073945 # another run: 4.274167584502551
    elseif startswith(model, "kilpisjarvi")
        5.956119530847897 # another run: 6.0015265040396
    elseif startswith(model, "diamonds")
        3.567778885451556 # another run: 3.5331787529624483
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
    df.miness_per_sec = df.miness ./ df.time
    df.miness_per_step = df.miness ./ df.n_steps
    df.miness_per_cost = ifelse.(
        map(sampler -> sampler in ["autoRWMH", "SliceSampler", "HitAndRunSlicer"], df.sampler_type),
        df.miness ./ df.n_logprob,  # non gradient-based samplers
        df.miness ./ (df.n_logprob .+ 2 * df.n_steps .* log_prob_gradient_ratio.(df.model)) # 1 leapfrog = 2 gradient eval
    ) # gradient based: we use cost = #log_potential_eval + eta * #gradient_eval, where eta is model dependent
    df.sampler = map(zip(df.sampler_type, df.selector, df.logstep_jitter)) do (t,s,j)
        t in ("NUTS", "SliceSampler") ? t : t * (s == "inverted" ? "_inv" : "") * (j in ["adapt", "normal"] ? "_jitter" : "")
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
    base_df = filter(row -> row.logstep_jitter == "adapt", base_df) # only consider stability for adaptive jitter
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
            ribbon=:jitter_std_sem, xlabel="n_steps", ylabel="Mean jitter std (log scale)", title="Mean + SE of jitter_std vs n_rounds",
            yaxis=:log10, markershape=:auto, linecolor=:auto) #, xaxis=:log2
    
        savefig(joinpath(plots_path,"$(sampler)_jitter_evolution.png"))
    end
end

#=
within-sampler comparison plots
=#
function sampler_comparison_plots_model(experiment::String)
    base_df = prepare_df(get_summary_df(experiment))
    plots_path = joinpath(base_dir(), "deliverables", experiment)

    unique_samplers = unique(base_df[:, :sampler_type])

    # Iterate over each sampler
    foreach(unique_samplers) do my_sampler
        # Filter the dataframe for the current sampler
        df = filter(:sampler_type => ==(my_sampler), base_df)
        sort!(df, :model) # ensure ordering on x-axis
        # Create the grouped boxplot
        @df df StatsPlots.groupedboxplot(:model, :miness_per_sec, group=:sampler, xlabel="Model", ylabel="minESS/second (log scale)", 
            legend=:topleft, title="minESS per second Comparison for $(my_sampler)", xrotation = 20, yaxis=:log10, color=:auto)
        savefig(joinpath(plots_path,"$(my_sampler)_miness_comparison.png"))
    end
end

#=
normal jitter comparison plots
=#
function jitter_comparison_plots_model(experiment::String)
    base_df = prepare_df(get_summary_df(experiment))
    base_df = filter(row -> !(row.logstep_jitter == "adapt" && row.n_rounds != 20), base_df)
    plots_path = joinpath(base_dir(), "deliverables", experiment)

    unique_samplers = unique(base_df[:, :sampler_type])

    # Iterate over each sampler
    foreach(unique_samplers) do my_sampler
        # Filter the dataframe for the current sampler
        df = filter(:sampler_type => ==(my_sampler), base_df)
        sort!(df, :model) # ensure ordering on x-axis
        # Create the grouped boxplot
        @df df StatsPlots.groupedboxplot(:model, :miness_per_sec, group=:logstep_jitter, xlabel="Model", ylabel="minESS / second(log scale)", 
            legend=:topleft, title="minESS per second by Jitter for $(my_sampler)", yaxis=:log10, color=:auto)
        savefig(joinpath(plots_path,"$(my_sampler)_miness_comparison.png"))
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
    @df df StatsPlots.groupedboxplot(:model, :miness_per_sec, group=:sampler_type, xlabel="Model", ylabel="minESS / second(log scale)", 
        legend=:outerbottom, title="Comparison of minESS per Second for All Samplers", color=:auto, yaxis=:log10)
    savefig(joinpath(plots_path,"miness_per_sec_comparison.png"))

    # now create the minESS/cost plot
    @df df StatsPlots.groupedboxplot(:model, :miness_per_cost, group=:sampler_type, xlabel="Model", ylabel="minESS / cost(log scale)", 
        legend=:outerbottom, title="Comparison of minESS per Cost for All Samplers", color=:auto, yaxis=:log10)
    savefig(joinpath(plots_path,"miness_per_cost_comparison.png"))
end