include { crossProduct; collectCSVs; setupPigeons } from './utils.nf'
params.dryRun = false

def variables = [
    seed: (1..30),
    model: ["eight_school_noncentered", "garch11", "gp_pois_regr", "lotka_volterra", "kilpisjarvi", "logearn_logheight_male", "diamonds"], // # parameters = 18, 4, 24, 90, 3, 4, 27
    sampler_type: ["SimpleAHMC", "SimpleRWMH", "NUTS", "HitAndRunSlicer"], //
    selector: ["standard", "inverted"],
    int_time: ["single_step", "rand"], // single_step gives autoMALA
    logstep_jitter: ["none", "fixed", "adapt"]
]

model_string = [
    eight_school_noncentered: "StanLogPotential(joinpath(\"$baseDir\", \"stan/eight_schools_noncentered.stan\"),joinpath(\"$baseDir\", \"data/eight_schools.json\"))", 
    garch11: "StanLogPotential(joinpath(\"$baseDir\", \"stan/garch11.stan\"),joinpath(\"$baseDir\", \"data/garch.json\"))",
    gp_pois_regr: "StanLogPotential(joinpath(\"$baseDir\", \"stan/gp_pois_regr.stan\"),joinpath(\"$baseDir\", \"data/gp_pois_regr.json\"))",
    lotka_volterra: "StanLogPotential(joinpath(\"$baseDir\", \"stan/lotka_volterra.stan\"),joinpath(\"$baseDir\", \"data/hudson_lynx_hare.json\"))",
    kilpisjarvi: "StanLogPotential(joinpath(\"$baseDir\", \"stan/kilpisjarvi.stan\"),joinpath(\"$baseDir\", \"data/kilpisjarvi_mod.json\"))",
    logearn_logheight_male: "StanLogPotential(joinpath(\"$baseDir\", \"stan/logearn_logheight_male.stan\"),joinpath(\"$baseDir\", \"data/earnings.json\"))",
    diamonds: "StanLogPotential(joinpath(\"$baseDir\", \"stan/diamonds.stan\"),joinpath(\"$baseDir\", \"data/diamonds.json\"))"
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
    template 'post_db_main.jl'
}

