#!/usr/bin/env -S julia --heap-size-hint=${task.memory.toGiga()}G
using Pkg
Pkg.activate(joinpath("$baseDir", "$julia_env")) 
include(joinpath("$baseDir", "$julia_env", "src", "utils.jl")) # loads dependencies too
# include(joinpath("$baseDir", "$julia_env", dirname(dirname(pathof(Pigeons))), "test", "supporting", "postdb.jl"))

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
	dim = 4
	scale = model in ["funnel(4,1)", "banana(4,1)"] ? 1 : 0.3
	target = ${model_string[arg.model]}
	seed = ${arg.seed}
	miness_threshold = ${params.dryRun ? 1 : 100}

	samples, stats_df = if explorer_type != "NUTS" # use Pigeons 
	    pt_sample_from_model(model, target, seed, explorer, miness_threshold)
	else
	    turing_nuts_sample_from_model(model, seed, miness_threshold; dim = dim, scale = scale) 
	end

	isdir("csvs") || mkdir("csvs")
	CSV.write("csvs/summary.csv", stats_df)
end

main()

