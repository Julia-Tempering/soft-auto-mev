include { crossProduct; collectCSVs; setupPigeons; head; pow; deliverables; checkGitUpdated; commit } from './utils.nf'
params.dryRun = false

def variables = [
    dim: (1..10).collect{1<<it}, // == 2^it but using bit shift,
    seed: (1..20),
    model: ["normal"],
    nleaps: [1],// currently not used// (10..30).step(10),
    sampler: ["AM","AH_fix", "AH_unif", "AH_exp","NUTS"]
]

model_string = [
    normal: "Pigeons.ScaledPrecisionNormalPath(1.0, 1.0, dim)"
]

sampler_string = [ 
    AM: "AutoMALA(base_n_refresh=1)",
    AH_fix: "SimpleAHMC()", // base_n_refresh=1 by default on pkg autoHMC. also jitter = Dirac(1.0)
    AH_unif: "SimpleAHMC(jitter_n_leaps=Uniform(0.,2.))",
    AH_exp: "SimpleAHMC(jitter_n_leaps=Exponential(1.))",
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
    // for some reason NUTS uses roughly 10 times more mem 
    memory { 2.GB +            // ~ fixed mem cost
        1.GB * Math.round(
            Math.min(90.0,     // smallest machines on beluga has 92G=90 + 2fixed
                task.attempt * (arg.sampler=="NUTS" ? 10 : 1) * arg.dim * 
                Math.pow(2.0, n_rounds+1.0 + 3.0 - 30.0) // (nrounds+1) + (log_2(bytes per float)) - log_2(bytes per gigabyte)
            )
        ) 
    }
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
