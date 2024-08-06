include { crossProduct; collectCSVs; setupPigeons } from './utils.nf'
params.dryRun = false

def variables = [
    dim: (1..10).collect{ 1 << it }, // bitshift
    seed: (1..30),
    model: ["banana", "normal", "funnel"],
    sampler_type: ["SimpleAHMC", "SimpleRWMH", "NUTS"],
    selector: ["standard", "inverted"],
    int_time: ["single_step", "rand"], // single_step gives autoMALA
    logstep_jitter: ["none", "normal"]
]

model_string = [
    normal: "Pigeons.ScaledPrecisionNormalPath(1.0, 1.0, dim)",
    funnel: "Pigeons.stan_funnel(dim-1, scale)", // NB: funnel and banana have extra parameter
    banana: "Pigeons.stan_banana(dim-1, scale)"
]

def MAX_RETRIES = params.dryRun ? 0 : 2 // workaround for retry-then-ignore: https://github.com/nextflow-io/nextflow/issues/1090#issuecomment-477964768
def julia_env_dir = file("julia-environment")
def julia_depot_dir = file(".depot")

workflow {
    args = crossProduct(variables, params.dryRun)
        .filter { it.sampler_type.startsWith("Simple") || it.selector == variables.selector.first() } // selector is only relevant for auto types
        .filter { it.sampler_type.startsWith("Simple") || it.logstep_jitter == variables.logstep_jitter.first() } // logstep_jitter is only relevant for auto types
        .filter { it.sampler_type == "SimpleAHMC" || it.int_time == variables.int_time.first() } // int_time is only relevant for autoHMC
    	// .view()  
    julia_env = setupPigeons(julia_depot_dir, julia_env_dir)
    agg_path = runSimulation(julia_depot_dir, julia_env, args) | collectCSVs
}

process runSimulation {
    memory { 1.GB * (6.0 + arg.dim * arg.dim * (90.0 / (1024.0*1024.0))) * task.attempt } // quad dim growth guess
    time { 1.hour * (0.5 + arg.dim * arg.dim * (4.5 / (1024.0*1024.0))) * task.attempt } // quad dim growth guess
    maxRetries { MAX_RETRIES }
    errorStrategy { task.attempt <= MAX_RETRIES ? 'retry' : 'ignore' }
    input:
        env JULIA_DEPOT_PATH
        path julia_env
        val arg
    output:
        tuple val(arg), path('csvs')
  script:
    template 'highdim_main.jl'
}

