include { crossProduct; collectCSVs; setupPigeons } from './utils.nf'
params.dryRun = false

def variables = [
    seed: (1..30),
    model: ["mRNA", "funnel(4,1)", "funnel(4,0.3)", "banana(4,1)", "banana(4,0.3)", "horseshoe_logit"], 
    sampler_type: ["SimpleAHMC", "SimpleRWMH", "HitAndRunSlicer", "NUTS"], //
    selector: ["inverted"], //no need to include "standard", i.e. non-inverted selector 
    int_time: ["single_step", "rand"], // single_step gives autoMALA
    logstep_jitter: ["adapt"] // no need to include "none" and "fixed"
]

model_string = [
    "funnel(4,1)": "stan_logpotential(\"funnel\"; dim = 4, scale = 1)", 
    "banana(4,1)": "stan_logpotential(\"banana\"; dim = 4, scale = 1)", 
    "funnel(4,0.3)": "stan_logpotential(\"funnel\"; dim = 4, scale = 0.3)",  
    "banana(4,0.3)": "stan_logpotential(\"banana\"; dim = 4, scale = 0.3)", 
    mRNA: "stan_logpotential(model)",
    horseshoe_logit: "stan_logpotential(model; dataset = \"sonar\")"
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
    memory { params.dryRun ? 4.GB : ( task.attempt * (8.GB * (arg.sampler_type == "NUTS" ? 4 : 1)) ) } // NUTS needs ~ 65M samples for 200 minESS
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
    template 'post_db_main.jl'
}

