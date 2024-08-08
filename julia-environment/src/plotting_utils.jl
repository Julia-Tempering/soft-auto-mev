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

# TODO: in preparation
function highdim_plots()
    experiment = "highdim"
    df = get_summary_df(experiment)
    df = filter(:miness => >=(100), df) # autoRWMH + inv + Gauss jitter fails bad on funnel
    df = prepare_df(df)

    # sym = :miness_per_step
    # model_df = filter(:model => ==(model), df)
    # model_df = filter(:dim => ==(1024), model_df)

    # fig_steps = data(model_df) * 
    #     visual(BoxPlot,orientation=:horizontal) *
    #     mapping(
    #         :sampler,
    #         sym => log10 => "minESS / " * (sym == :miness_per_sec ? "second" : "step") * " (log₁₀)",
    #         color=:sampler_type
    #     ) |> 
    #     draw(scales(Color = (; palette = :seaborn_colorblind)), axis = (width = 400, height = 400))

    model = "normal"
    dim_df = filter(:model => ==(model), df)
    dim_df = filter(:sampler => in(("autoRWMH_inv","autoHMC_inv","autoMALA_inv","NUTS")), dim_df)
    dim_df.marginess_per_sec = dim_df.margin_ess_exact ./ dim_df.time
    dim_df.marginess_per_step = dim_df.margin_ess_exact ./ dim_df.n_steps

    fig_steps = data(dim_df) * 
        visual(BoxPlot) *
        mapping(
            :dim => log2 => "Dimension (log₂)",
            :marginess_per_sec => log10 => "Exact ESS per second (log₁₀)",
            color=:sampler_type,
            dodge=:sampler_type
        )
    fig_steps=draw(fig_steps,axis = (width = 800, height = 400))
end
