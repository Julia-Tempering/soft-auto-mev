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
	target = ${model_string[arg.model]}
	seed = ${arg.seed}
	miness_threshold = ${params.dryRun ? 1 : 100}

	samples, stats_df = pt_sample_from_model(model, target, seed, explorer, miness_threshold)

	isdir("csvs") || mkdir("csvs")
	CSV.write("csvs/summary.csv", stats_df)
end

main()