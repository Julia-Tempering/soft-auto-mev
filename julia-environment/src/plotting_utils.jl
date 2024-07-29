using Colors
using DataFrames
using Plots
using Plots.PlotMeasures: px
using StatsPlots

include("utils.jl")

function get_summary_df(experiment::String)
    base_folder = base_dir()
    csv_path    = joinpath(base_folder, "deliverables", experiment, "aggregated", "summary.csv")
    return DataFrame(CSV.File(csv_path))
end

###############################################################################
# plotting utilities
###############################################################################

# Boxplots
const JULIA_AUTO = theme_palette(:auto).colors.colors
const COLORBLIND_IBM = [colorant"#785EF0", colorant"#DC267F", colorant"#FE6100", colorant"#FFB000", colorant"#648FFF"] # IBM palette from https://davidmathlogic.com/colorblind/
const COLORBLIND_WONG = [
    colorant"#000000", colorant"#E69F00", colorant"#56B4E9", colorant"#009E73", 
    colorant"#F0E442", colorant"#0072B2", colorant"#D55E00", colorant"#CC79A7"] # Wong palette from https://davidmathlogic.com/colorblind/
get_palette(n_levels::Int) = 
    ((n_levels <= 2 || n_levels > 8) ? JULIA_AUTO : n_levels <= 5 ? COLORBLIND_IBM : COLORBLIND_WONG)[1:n_levels]
dataset_nickname(d::AbstractString) = d=="sonar" ? "Sonar" : (d=="prostate_small" ? "Prostate" : "Ion.")
function make_boxplots(df::DataFrame; model = first(df.model), fn_end = ".pdf")
    n_samplers = length(unique(df.sampler))
    only_two_samplers = n_samplers == 2
    colors = get_palette(n_samplers)

    # preprocessing
    sort!(df)
    is_funnel = occursin("funnel",model)
    is_highdim = occursin("highdim",model)
    is_banana = !is_highdim && occursin("banana",model) # check if its banana_scale
    is_log2_x = is_banana || is_highdim
    is_log2_x && (df.dim .= log2.(df.dim)) # boxplot not working with xaxis=:log 
    is_hsp = occursin("horseshoe",model) && hasproperty(df, :dataset)
    is_2_comp_norm = occursin("two_component_normal",model)
    xvar = :dim
    if is_hsp
        xvar = :data_nobs
        df[!,xvar] = map(
            t -> (dataset_nickname(t[1]) * "\n" * (t[2] > 100 ? "full" : string(t[2]))),
            zip(df.dataset,df.dim)
        )
        sort!(df, xvar)
    end
    df[!,:min_ess] = min.(df.ess, df.ess_exact)
    df[!,:nleap_to_min_ess] = df.n_leapfrog ./ df.min_ess
    df_means = combine(
        groupby(df, [xvar, :sampler]),
        :margin1_mean => mean,
        :margin1_var => mean,
        :n_leapfrog => mean,
        :min_ess => mean,
        :nleap_to_min_ess => mean, 
        renamecols=false)
    sort!(df_means)

    # common properties
    size   = (650,300)
    xlab   = is_hsp ? "Dataset" : (is_highdim ? "Dimension" : "Inverse scale") * (is_log2_x ? " (log₂)" : "")
    mar    = 15px
    path   = joinpath(base_dir(), "deliverables", model)
    n_samplers = length(unique(df.sampler))

    # plots for known margins
    if !is_hsp
        # margin1 mean
        margin_idx = is_2_comp_norm ? 2 : 1
        p=@df df groupedboxplot(
            :dim, 
            :margin1_mean, 
            group=:sampler,
            bar_position = :dodge,
            size=size,
            xlab=xlab,
            palette=colors,
            ylab="Margin $margin_idx mean",
            left_margin = mar, bottom_margin = mar,
        )
        if !only_two_samplers
            @df df_means plot!(p,
                :dim,
                :margin1_mean,
                group=:sampler,
                palette=colors,
                linewidth=2,
                label=""
            )
        end
        savefig(p,joinpath(path, "boxplots-margin-mean" * fn_end))

        # margin1 var
        p=@df df groupedboxplot(
            :dim, 
            :margin1_var, 
            group=:sampler,
            bar_position = :dodge, 
            size=size,
            xlab=xlab,
            legend = is_2_comp_norm ? :topleft : :best,
            yaxis= is_2_comp_norm ? :log : :identity,
            palette=colors,
            ylab="Margin $margin_idx var",
            left_margin = mar, bottom_margin = mar,
        )
        if !only_two_samplers
            @df df_means plot!(p,
                :dim,
                :margin1_var,
                group=:sampler,
                yaxis= is_2_comp_norm ? :log : :identity,
                palette=colors,
                linewidth=2,
                label=""
            )
        end
        if is_2_comp_norm
            dim_vals = sort(unique(df.dim))
            plot!(dim_vals, 10. .^ (2*dim_vals), label = "true",
            linestyle=:dash, color=colorant"#648FFF")
        end
        savefig(p,joinpath(path, "boxplots-margin-var" * fn_end))
    end

    # n_leapfrog
    p=@df df groupedboxplot(
        (is_hsp ? :data_nobs : :dim), 
        :n_leapfrog, 
        group=:sampler,
        bar_position = :dodge, 
        legend = is_hsp || is_2_comp_norm ? :outerright : :best,
        yaxis= :log,
        palette=colors,
        size=size,
        xlab=xlab,
        ylab="Total number of leapfrog steps",
        left_margin = mar, bottom_margin = mar,
    )
    if !only_two_samplers && is_hsp
        plot!(p,
            repeat(first(first(xticks(p))),inner=n_samplers),
            df_means.n_leapfrog,
            group=df_means.sampler,
            yaxis=:log,
            palette=colors,
            linewidth=2,
            label=""
        )
    elseif !only_two_samplers
        @df df_means plot!(p,
            :dim,
            :n_leapfrog,
            group=:sampler,
            yaxis=:log,
            palette=colors,
            linewidth=2,
            label=""
        )
    end
    savefig(p,joinpath(path, "boxplots-n_leapfrog" * fn_end))

    # min ess
    p=@df df groupedboxplot(
        (is_hsp ? :data_nobs : :dim), 
        :min_ess,
        group=:sampler,
        bar_position = :dodge, 
        legend = is_2_comp_norm ? :outerright : :best,
        yaxis= :log,
        palette=colors,
        size=size,
        xlab=xlab,
        ylab="minESS",
        left_margin = mar, bottom_margin = mar
    )
    if !only_two_samplers && is_hsp
        plot!(p,
            repeat(first(first(xticks(p))),inner=n_samplers),
            df_means.min_ess,
            group=df_means.sampler,
            yaxis=:log,
            palette=colors,
            linewidth=2,
            label=""
        )
    elseif !only_two_samplers
        @df df_means plot!(p,
            :dim,
            :min_ess,
            group=:sampler,
            yaxis=:log,
            palette=colors,
            linewidth=2,
            label=""
        )
    end
    savefig(p,joinpath(path, "boxplots-miness" * fn_end))

    # nleap to miness
    p=@df df groupedboxplot(
        (is_hsp ? :data_nobs : :dim), 
        :nleap_to_min_ess, 
        group=:sampler,
        bar_position = :dodge, 
        legend = (only_two_samplers || is_funnel) ? :bottomright : (is_2_comp_norm ? :topleft : (is_hsp ? :outerright : :best)),
        yaxis= :log,
        palette=colors,
        size=size,
        xlab=xlab,
        ylab="Leapfrogs per minESS",
        left_margin = mar, bottom_margin = mar,
    )
    if !only_two_samplers && is_hsp
        plot!(p,
            repeat(first(first(xticks(p))),inner=n_samplers),
            df_means.nleap_to_min_ess,
            group=df_means.sampler,
            yaxis=:log,
            palette=colors,
            linewidth=2,
            label=""
        )
    elseif !only_two_samplers
        @df df_means plot!(p,
            :dim,
            :nleap_to_min_ess,
            group=:sampler,
            yaxis=:log,
            palette=colors,
            linewidth=2,
            label=""
        )
    end
    savefig(p,joinpath(path, "boxplots-nleap_to_min_ess" * fn_end))
end


