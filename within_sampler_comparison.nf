include { crossProduct; collectCSVs; setupPigeons } from './utils.nf'
params.dryRun = false

def variables = [
    seed: (1..30),
    model: ["funnel_2d_1c", "funnel_128d_1c", "funnel_128d_10c", "funnel_2d_10c", "banana_2d_1c", "banana_128d_1c", "banana_128d_10c", "banana_2d_10c", "normal_2d_1c", "normal_128d_1c", "normal_128d_10c", "normal_2d_10c"],
    sampler_type: ["SimpleAHMC", "SimpleRWMH"], 
    selector: ["standard", "inverted"],
    int_time: ["single_step", "rand"], // single_step gives autoMALA
    logstep_jitter: ["none", "adapt"]
]

model_string = [
    funnel_2d_1c: "Pigeons.stan_funnel(1, 1.0)",
    funnel_128d_1c: "Pigeons.stan_funnel(127, 1.0)",
    funnel_128d_10c: "Pigeons.stan_funnel(127, 10.0)",
    funnel_2d_10c: "Pigeons.stan_funnel(1, 10.0)",
    banana_2d_1c: "Pigeons.stan_banana(1, 1.0)",
    banana_128d_1c: "Pigeons.stan_banana(127, 1.0)",
    banana_128d_10c: "Pigeons.stan_banana(127, 10.0)",
    banana_2d_10c: "Pigeons.stan_banana(1, 10.0)",
    normal_2d_1c: "Pigeons.ScaledPrecisionNormalPath(1.0, 1.0, 2)",
    normal_128d_1c: "Pigeons.ScaledPrecisionNormalPath(1.0, 1.0, 128)",
    normal_128d_10c: "Pigeons.ScaledPrecisionNormalPath(10.0, 1.0, 128)",
    normal_2d_10c: "Pigeons.ScaledPrecisionNormalPath(10.0, 1.0, 2)"
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
    template 'within_sampler_comparison_main.jl'
}