#using Pkg
#Pkg.activate("julia-environment") 
#include(joinpath("../julia-environment/src/utils.jl")) 
using CairoMakie
using LaTeXStrings

function get_explorer(explorer)
    explorer == "RWMH" ? 
    SimpleRWMH(
            step_size = 4 * 0.4179156164274893, 
            estimated_target_std_deviations = [2.7712018780501895, 116.8387973797307, 101.13224298644053, 89.71676443117684],
            step_size_selector = autoRWMH.MHNonAdaptiveSelector(),
            step_jitter = autoRWMH.StepJitter(dist = Dirac(0.0), adapt_strategy = autoRWMH.FixedStepJitter())
    ) : SimpleRWMH()
end

function get_sample_scatter(explorer, ax)
    pt = pigeons(
        target     = Pigeons.stan_funnel(3,1.0),
        seed       = 1,
        n_rounds   = 15,
        n_chains   = 1, 
        record     = [record_default(); Pigeons.traces; online],
        explorer   = get_explorer(explorer),
        show_report = true
    )
    
    samples = vcat(get_sample(pt))
    
    # Extract x coordinates
    x = [sample[1] for sample in samples]
    y = [sample[2] for sample in samples]
    
    # Plot the trace line
    # lines!(ax, x, y, alpha=0.2)

    # Add scatter points (dots) around each x-value
    CairoMakie.scatter!(ax, x, y, alpha=0.07, markersize=4)
    # Set x-axis limits
    # CairoMakie.ylims!(ax, -5, 5)
end

# Initialize the figure with a grid layout
fig = Figure(resolution = (800, 400))

# Create two axes side by side
ax1 = Axis(fig[1, 1], xlabel=L"x_1", ylabel=L"x_2")
ax2 = Axis(fig[1, 2], xlabel=L"x_1", ylabel=L"x_2")

# Generate and plot the samples for each selector on separate axes
get_sample_scatter("RWMH", ax1)
get_sample_scatter("autoRWMH", ax2)

# Display the figure with both trace plots side by side
# display(fig)
save("deliverables/trace_plot/RWMH_scatter.png", fig)