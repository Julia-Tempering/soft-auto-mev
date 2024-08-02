using AlgebraOfGraphics, CairoMakie

include("utils.jl")

function get_summary_df(experiment::String)
    base_folder = base_dir()
    csv_path    = joinpath(base_folder, "deliverables", experiment, "aggregated", "summary.csv")
    return DataFrame(CSV.File(csv_path))
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
    df.sampler = map(zip(df.sampler_type, df.selector, df.logstep_jitter)) do (t,s,j)
        t in ("NUTS", "SliceSampler") ? t : t * (s == "inverted" ? "_inv" : "") * (j == "normal" ? "_jitter" : "")
    end
    return df
end

function scatter_miness_cost(experiment::String)
    df = filter(:scale_idx => ==(1), get_summary_df(experiment))
    df = prepare_df(df)
    plots_path = joinpath(base_dir(), "deliverables", experiment)

    foreach((:miness_per_sec, :miness_per_step)) do sym
        fig_steps = data(df) * 
            visual(Scatter) *
            mapping(
                sym => log10 => "minESS / " * (sym == :miness_per_sec ? "second" : "step") * " (log scale)",
                :sampler => identity => "Sampler",
                color=:sampler_type
            ) |> 
            draw(axis = (width = 400, height = 400))
        save(joinpath(plots_path,"$sym.png"), fig_steps, px_per_unit = 3) # save high-resolution png
    end
end