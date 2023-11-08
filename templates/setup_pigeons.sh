#!/bin/bash

# if insider the apptainer container, start ssh agent and add key
if [ -v APPTAINER_NAME ] ; then
    echo "initializing ssh agent and adding deploy key"
    chmod 600 $baseDir/keys/id_ed25519
    eval `ssh-agent -s`
    ssh-add $baseDir/keys/id_ed25519
fi

# need to do this here so that nextflow can fill in julia_env and baseDir
# also, cannot run the jl script in another process because then the container is reloaded
# so the ssh setup disappears
cat <<EOF > temp.jl
    using Pkg
    Pkg.activate("$julia_env")
    Pkg.add(name="Pigeons", rev="main")
    Pkg.add(url="git@github.com:Julia-Tempering/autoHMC.git", rev="main")
    Pkg.instantiate()
    Pkg.precompile()

    # force download of BridgeStan artifact because it uses LazyArtifact
    # no other package uses this approach
    # also, need to precompile stan models to avoid race conditions since the
    # .o and .so files are stored in a shared location (pigeons examples folder)
    using BridgeStan
    using Pigeons
    Pigeons.toy_stan_target(1)
    Pigeons.stan_funnel(1)
    Pigeons.stan_banana(1)
    StanLogPotential(
        joinpath("$baseDir", "stan", "horseshoe_logit.stan"), 
        Pigeons.json(; n=0,d=0,x="[[]]",y="[]")
    )
    StanLogPotential(
        joinpath("$baseDir", "stan", "two_component_normal.stan"), 
        Pigeons.json(; n=1, s_lo=0.1, s_hi=10.0)
    )
EOF

julia temp.jl

