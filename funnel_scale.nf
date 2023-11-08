include { crossProduct; collectCSVs; addDeployKey; setupPigeons; head; pow; deliverables; checkGitUpdated; commit } from './utils.nf'
params.dryRun = false

def variables = [
    dim: (1..20).collect{0.2*it},
    seed: (1..20),
    model: ["funnel_scale"],
    nleaps: (0..4).collect{1<<it}, // == 2^it but using bit shift
    sampler: ["AM","AH_simple","NUTS"] // MALA runs alongside autoMALA 
]

model_string = [
    funnel_scale: "Pigeons.stan_funnel(1, 1/dim)", 
]

sampler_string = [ 
    AM: "AutoMALA()",
    AHMC: "SimpleAHMC(n_leaps = nleaps)",
    NUTS: "Pigeons.MALA()", // ignored, just use it to compile
]

n_rounds = params.dryRun ? 4 : 20 
PT_n_chains = 10
def julia_env_dir = file("julia-environment")
def julia_depot_dir = file(".depot")
def deliv = deliverables(workflow)

workflow {
    args = crossProduct(variables, params.dryRun)
        .filter { it.sampler.startsWith("AH") || it.nleaps == 1 } // nleaps only relevant to AHMC
    julia_env_path = addDeployKey(julia_env_dir)
    julia_env = setupPigeons(julia_depot_dir, julia_env_path)
    agg_path = runSimulation(julia_depot_dir, julia_env, args) | collectCSVs 
    //commit(agg_path, params.dryRun) // cannot commit from container, priv keys not available
}

process runSimulation {
    memory {4.GB * task.attempt}
    time { 1.hour * task.attempt }
    errorStrategy 'retry'
    maxRetries '3'
    input:
        env JULIA_DEPOT_PATH
        path julia_env
        val arg
    output:
        tuple val(arg), path('csvs')
  script:
    template 'AM_scale_main.jl'
}

