#!/usr/bin/env -S julia --heap-size-hint=${task.memory.toGiga()}G
using Pkg
Pkg.activate(joinpath("$baseDir", "$julia_env")) 
include(joinpath("$baseDir", "$julia_env", "src", "utils.jl")) # loads dependencies too

function main()
	# collect global vars 
	nleaps = ${arg.nleaps}
	explorer_type = "${arg.sampler}"
	explorer = ${sampler_string[arg.sampler]}
	dim = ${arg.dim} # interpreted as scale
	target = ${model_string[arg.model]}
	model = "${arg.model}"
	seed = ${arg.seed}
	n_rounds = ${n_rounds}
	#model == "funnel_scale" && occursin("autoMALA", explorer_type) && (n_rounds -= 1) # equalize efforts

	time, sample, n_leapfrog = if explorer_type != "NUTS" # use Pigeons 
	    pt_sample_from_model(target, seed, explorer, explorer_type, n_rounds; n_chains = n_chains)
	else # use cmdstan for NUTS
	    nuts_sample_from_model(model, seed, n_rounds; dim=dim)     
	end
	idx_margin = special_margin_idx(model,n_vars(sample))
	ess = margin_ess(sample)
	ess_exact = margin_ess(sample, model, idx_margin, dim) # dim is passed to special_margin_mean_std
	margin1_mean = mean(sample, idx_margin) 
	margin1_var = var(sample, idx_margin) 
	df = DataFrame(
	    ess = ess, ess_exact = ess_exact, n_leapfrog = n_leapfrog,
	    margin1_mean = margin1_mean, margin1_var = margin1_var, 
	)

	!isdir("csvs") ? mkdir("csvs") : nothing
	CSV.write("csvs/summary.csv", df)
end

main()

