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
		explorer_type, "${arg.selector}", "${arg.int_time}", "${arg.logstep_jitter}"
	)
	model = "${arg.model}"
	scale = get_scale(${arg.scale_idx}, model)
	target = ${model_string[arg.model]}
	seed = ${arg.seed}
	miness_threshold = ${params.dryRun ? 1 : 100}

	time, samples, n_steps, miness, prop_log_diff = if explorer_type != "NUTS" # use Pigeons 
	    pt_sample_from_model(model, target, seed, explorer, miness_threshold)
	else # use cmdstan for NUTS
	    nuts_sample_from_model(model, seed, miness_threshold; scale=scale)
	end

	df = hcat(
		DataFrame(scale = scale, miness = miness, time = time, n_steps = n_steps, prop_log_diff = prop_log_diff),
		margin_summary(samples, model, scale)
	)
	isdir("csvs") || mkdir("csvs")
	CSV.write("csvs/summary.csv", df)
end

main()

