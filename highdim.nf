include { crossProduct; collectCSVs; setupPigeons } from './utils.nf'
params.dryRun = false

def variables = [
    dim: (1..20).collect{ 1 << it }, // bitshift
    seed: (1..30),
    model: ["normal", "banana", "funnel"],
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

def MAX_RETRIES = 2 // workaround for retry-then-ignore: https://github.com/nextflow-io/nextflow/issues/1090#issuecomment-477964768
def julia_env_dir = file("julia-environment")
def julia_depot_dir = file(".depot")

workflow {
    args = crossProduct(variables, params.dryRun)
        .filter { it.model == "normal" || it.dim <= 1024 } // only the Gaussian case is feasible in higher dimensions
        .filter { it.sampler_type.startsWith("Simple") || it.selector == variables.selector.first() } // selector is only relevant for auto types
        .filter { it.sampler_type == "SimpleRWMH" || it.logstep_jitter == variables.logstep_jitter.first() } // using logstep_jitter only for RWMH (fails in AHMC outside of Gaussian setting)
        .filter { it.sampler_type == "SimpleAHMC" || it.int_time == variables.int_time.first() } // int_time is only relevant for autoHMC
    	// .view()  
    julia_env = setupPigeons(julia_depot_dir, julia_env_dir)
    agg_path = runSimulation(julia_depot_dir, julia_env, args) | collectCSVs
}

process runSimulation {
    memory { 1.GB * (2.0 + (arg.model == "normal" ? 1 : 1024) * 30.0 *(arg.dim/1048576.0)*(arg.dim/1048576.0)) * task.attempt } // quad dim growth guess (32G @ 2^20 dims) * 1000 if non-Gaussian targets (poor ESS ratio)
    time { 1.hour * (0.5 + (arg.model == "normal" ? 1 : 1024) * 4.5  *(arg.dim/1048576.0)*(arg.dim/1048576.0)) * task.attempt } // similar
    maxRetries { MAX_RETRIES }
    errorStrategy { params.dryRun ? 'finalize' : (task.attempt <= MAX_RETRIES ? 'retry' : 'ignore') }
    input:
        env JULIA_DEPOT_PATH
        path julia_env
        val arg
    output:
        tuple val(arg), path('csvs')
  script:
    template 'highdim_main.jl'
}

