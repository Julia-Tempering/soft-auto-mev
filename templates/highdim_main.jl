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
    dim = ${arg.dim}
	model = "${arg.model}"
	scale = get_scale(1, model)
	target = ${model_string[arg.model]}
	seed = ${arg.seed}
	miness_threshold = ${params.dryRun ? 1 : 100}

	time, samples, n_steps, miness = if explorer_type != "NUTS" # use Pigeons 
	    pt_sample_from_model(model, target, seed, explorer, miness_threshold; max_rounds=21) # autoRWMH breaks down at highdims, need to stop at some point or it keeps going and eating RAM 
	else # use cmdstan for NUTS
	    nuts_sample_from_model(model, seed, miness_threshold; dim=dim, scale=scale)
	end

	df = hcat(
		DataFrame(scale = scale, miness = miness, time = time, n_steps = n_steps),
		margin_summary(samples, model, scale)
	)
	isdir("csvs") || mkdir("csvs")
	CSV.write("csvs/summary.csv", df)
end

main()

