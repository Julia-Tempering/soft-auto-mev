include { crossProduct; collectCSVs; setupPigeons } from './utils.nf'
params.dryRun = false

def variables = [
    seed: (1..30),
    model: ["earn_height", "nes", "diamonds", "hmm_example"], // # parameters = 3, 10, 27, 111
    sampler_type: ["SimpleAHMC", "SimpleRWMH", "NUTS", "SliceSampler"], //
    selector: ["standard", "inverted"],
    int_time: ["single_step", "rand"], // single_step gives autoMALA
    logstep_jitter: ["none", "normal"]
]

model_string = [
    earn_height: "StanLogPotential(joinpath(\"$baseDir\", \"data/earn_height.stan\"),joinpath(\"$baseDir\", \"data/earnings.json\"))", 
    nes: "StanLogPotential(joinpath(\"$baseDir\", \"data/nes.stan\"),joinpath(\"$baseDir\", \"data/nes2000.json\"))",
    diamonds: "StanLogPotential(joinpath(\"$baseDir\", \"data/diamonds.stan\"),joinpath(\"$baseDir\", \"data/diamonds.json\"))", 
    hmm_example: "StanLogPotential(joinpath(\"$baseDir\", \"data/hmm_example.stan\"),joinpath(\"$baseDir\", \"data/hmm_example.json\"))"
]

def MAX_RETRIES = params.dryRun ? 0 : 1 // workaround for retry-then-ignore: https://github.com/nextflow-io/nextflow/issues/1090#issuecomment-477964768
def julia_env_dir = file("julia-environment")
def julia_depot_dir = file(".depot")

workflow {
    args = crossProduct(variables, params.dryRun)
        .filter { it.sampler_type.startsWith("Simple") || it.selector == variables.selector.first() } // selector is only relevant for auto types
        .filter { it.sampler_type == "SimpleRWMH" || it.logstep_jitter == variables.logstep_jitter.first() } // using logstep_jitter only for RWMH (fails in AHMC outside of Gaussian setting)
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
    template 'post_db_main.jl'
}

