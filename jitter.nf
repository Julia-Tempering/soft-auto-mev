include { crossProduct; collectCSVs; setupPigeons } from './utils.nf'
params.dryRun = false

def variables = [
    seed: (1..30),
    model: ["funnel", "banana", "normal", "kilpisjarvi", "logearn_logheight_male"],
    sampler_type: ["SimpleAHMC", "SimpleRWMH"], //
    selector: ["standard", "inverted"],
    int_time: ["single_step", "rand"], // single_step gives autoMALA
    logstep_jitter: ["none", "0.1", "0.5", "2.0", "4.0", "adapt"]
]

model_string = [
    normal: "Pigeons.ScaledPrecisionNormalPath(1.0, 1.0, 20)",
    funnel: "Pigeons.stan_funnel(3, 3.0)", // NB: funnel and banana have extra parameter
    banana: "Pigeons.stan_banana(3, 3.0)",
    kilpisjarvi: "stan_logpotential(model)",
    logearn_logheight_male: "stan_logpotential(model)"
]

def MAX_RETRIES = params.dryRun ? 0 : 1 // workaround for retry-then-ignore: https://github.com/nextflow-io/nextflow/issues/1090#issuecomment-477964768
def julia_env_dir = file("julia-environment")
def julia_depot_dir = file(".depot")

workflow {
    args = crossProduct(variables, params.dryRun)
        .filter { it.sampler_type.startsWith("Simple") || it.selector == variables.selector.first() } // selector is only relevant for auto types
        .filter { it.sampler_type.startsWith("Simple") || it.logstep_jitter == variables.logstep_jitter.first() } // step jitter only relevant for auto types
        .filter { it.sampler_type == "SimpleAHMC" || it.int_time == variables.int_time.first() } // int_time is only relevant for autoHMC
    	// .view()  
    julia_env = setupPigeons(julia_depot_dir, julia_env_dir)
    agg_path = runSimulation(julia_depot_dir, julia_env, args) | collectCSVs
}

process runSimulation {
    memory { params.dryRun ? 4.GB : ( task.attempt * (8.GB * (arg.sampler_type == "NUTS" ? 2 : 1)) ) } // NUTS needs ~ 65M samples for 200 minESS
    time { 1.hour * task.attempt * (arg.sampler_type == "NUTS" ? 4 : 1) } // same reason as above
    maxRetries { MAX_RETRIES }
    errorStrategy { task.attempt <= MAX_RETRIES ? 'retry' : 'ignore' }
    input:
        env JULIA_DEPOT_PATH
        path julia_env
        val arg
    output:
        tuple val(arg), path('csvs')
  script:
    template 'jitter_main.jl'
}