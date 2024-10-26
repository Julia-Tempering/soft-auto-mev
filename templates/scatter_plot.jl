using Pkg
Pkg.activate("julia-environment") 
include(joinpath("../julia-environment/src/utils.jl")) 
using CairoMakie
using LaTeXStrings

function get_selector(selector)
    selector = selector == "autoRWMH_inv" ? AutoStep.ASSelectorInverted() : AutoStep.ASSelector()
    return selector
end

function get_sample_scatter(selector, ax)
    pt = pigeons(
        target     = Pigeons.ScaledPrecisionNormalPath(2.0, 1.0, 2),
        seed       = 1,
        n_rounds   = 13,
        n_chains   = 1, 
        record     = [record_default(); Pigeons.traces; online],
        explorer   = SimpleRWMH(step_size_selector = get_selector(selector),
                                step_jitter = AutoStep.StepJitter(dist = Dirac(0.0),
                                adapt_strategy = AutoStep.FixedStepJitter())),
        show_report = true
    )
    
    samples = vcat(get_sample(pt))
    
    # Extract x coordinates
    x = [sample[1] for sample in samples]
    y = [sample[2] for sample in samples]
    
    # Plot the trace line
    # lines!(ax, x, y, alpha=0.2)

    # Add scatter points (dots) around each x-value
    CairoMakie.scatter!(ax, x, y, alpha=0.2, markersize=6)
    # Set x-axis limits
    CairoMakie.ylims!(ax, -5, 5)
end

# Initialize the figure with a grid layout
fig = Figure(resolution = (800, 400))

# Create two axes side by side
ax1 = Axis(fig[1, 1], xlabel=L"x_1", ylabel=L"x_2")
ax2 = Axis(fig[1, 2], xlabel=L"x_1", ylabel=L"x_2")

# Generate and plot the samples for each selector on separate axes
get_sample_scatter("autoRWMH", ax1)
get_sample_scatter("autoRWMH_inv", ax2)

# Display the figure with both trace plots side by side
# display(fig)
save("deliverables/trace_plot/gaussian_trace_plot.png", fig)