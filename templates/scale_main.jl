#!/usr/bin/env -S julia --heap-size-hint=${task.memory.toGiga()}G
using Pkg
Pkg.activate(joinpath("$baseDir", "$julia_env")) 
include(joinpath("$baseDir", "$julia_env", "src", "utils.jl")) # loads dependencies too

@info """ Running experiments with the following combination of parameters:
	
	${arg}

"""

function main()
	# collect global vars 
	explorer_type = "${arg.sampler_type}"
	explorer = make_explorer(
		"${arg.sampler_type}", "${arg.selector}", "${arg.int_time}", "${arg.logstep_jitter}"
	)
	model = "${arg.model}"
	scale = get_scale(${arg.scale_idx}, model)
	target = ${model_string[arg.model]}
	seed = ${arg.seed}
	min_ess_threshold = ${min_ess_threshold}

	time, samples, n_steps, miness = if explorer_type != "NUTS" # use Pigeons 
	    pt_sample_from_model(target, seed, explorer, min_ess_threshold)
	else # use cmdstan for NUTS
	    nuts_sample_from_model(model, seed, min_ess_threshold; scale=scale)
	end
	margin_idx = 1
	margin_ess_exact = margin_ess(samples, model, margin_idx, scale)
	margin_samples = get_component_samples(samples, margin_idx)
	margin_mean = mean(margin_samples) 
	margin_var = var(margin_samples) 
	df = DataFrame(
	    scale = scale, miness = miness, margin_ess_exact = margin_ess_exact, time = time,
		n_steps = n_steps, margin_mean = margin_mean, margin_var = margin_var, 
	)

	isdir("csvs") || mkdir("csvs")
	CSV.write("csvs/summary.csv", df)
end

main()

