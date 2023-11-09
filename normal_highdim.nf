include { crossProduct; collectCSVs; setupPigeons; head; pow; deliverables; checkGitUpdated; commit } from './utils.nf'
params.dryRun = false

def variables = [
    dim: (1..10).collect{1<<it}, // == 2^it but using bit shift,
    seed: (1..20),
    model: ["normal"],
    nleaps: (10..30).step(10),
    sampler: ["AM","AH_simple","NUTS"]
]

model_string = [
    normal: "Pigeons.ScaledPrecisionNormalPath(1.0, 1.0, dim)"
]

sampler_string = [ 
    AM: "AutoMALA(base_n_refresh=1)",
    AH_simple: "SimpleAHMC(n_leaps = nleaps, base_n_refresh=1)",
    NUTS: "Pigeons.MALA()", // ignored, just use it to compile
]

scale = 1.0 // keep, otherwise AM_highdim_main.jl will crash
n_rounds = params.dryRun ? 4 : 18
def julia_env_dir = file("julia-environment")
def julia_depot_dir = file(".depot")
def deliv = deliverables(workflow)

workflow {
    args = crossProduct(variables, params.dryRun)
        .filter { it.sampler.startsWith("AH") || it.nleaps == variables.nleaps.first() } // nleaps only relevant to AHMC
    	//.collect()
    	//.view()    
    julia_env = setupPigeons(julia_depot_dir, julia_env_dir)
    agg_path = runSimulation(julia_depot_dir, julia_env, args) | collectCSVs
}

process runSimulation {
    // linearly scale mem with dim*(total number of samples when doing n_rounds = 2^(n_rounds+1)-2 ~ 2^(n_rounds+1))
    // reference is 5G for 1024 dim and 18 rounds for AH_simple
    // NUTS uses roughly 8 times more mem 
    memory { 1.GB * Math.round(
        Math.min(92.0,     // smallest machines on beluga
            Math.max(2.0,   // ~ fixed mem cost
                5.0 * (arg.sampler=="NUTS" ? 8 : 1) * task.attempt * arg.dim * Math.pow(2,n_rounds+1)  / (1024.0 * 524288)
            )
        )
    ) }
    time { 2.hour * task.attempt }
    errorStrategy 'retry'
    maxRetries '1'
    input:
        env JULIA_DEPOT_PATH
        path julia_env
        val arg
    output:
        tuple val(arg), path('csvs')
    script: 
        template 'highdim_main.jl'
}
