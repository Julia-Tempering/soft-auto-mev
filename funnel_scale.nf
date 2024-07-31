include { crossProduct; collectCSVs; setupPigeons } from './utils.nf'
params.dryRun = false

def variables = [
    scale_idx: (1..20),
    seed: (1..30),
    model: ["funnel_scale"],
    sampler_type: ["SimpleAHMC", "SimpleRWMH", "NUTS", "SliceSampler"],
    selector: ["standard", "inverted"],
    int_time: ["single_step", "rand"], // single_step gives autoMALA
    logstep_jitter: ["none", "normal"]
]

model_string = [
    funnel_scale: "Pigeons.stan_funnel(1, scale)",
]

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
    memory { params.dryRun ? 4.GB : (16.GB * task.attempt) }
    time { 2.hour * task.attempt }
    errorStrategy { params.dryRun ? 'terminate' : ( (task.attempt <= process.maxRetries) ? 'retry' : 'ignore' ) } // retry-then-ignore strategy
    maxRetries 1
    input:
        env JULIA_DEPOT_PATH
        path julia_env
        val arg
    output:
        tuple val(arg), path('csvs')
  script:
    template 'scale_main.jl'
}

