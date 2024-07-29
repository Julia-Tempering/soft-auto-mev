#!/bin/bash

# if inside the apptainer container, start ssh agent and add keys
if [ -v APPTAINER_NAME ] ; then
    echo "initializing ssh agent and adding deploy keys"
    chmod 600 $baseDir/keys/id_autoHMC
    chmod 600 $baseDir/keys/id_autoRWMH
    eval `ssh-agent -s`
    ssh-add $baseDir/keys/id_autoHMC
    ssh-add $baseDir/keys/id_autoRWMH
fi

# Create Julia script that activates the julia environment and adds latest dev
# versions of our packages.
# NB: need to do this here so that nextflow can fill in julia_env and baseDir
# also, cannot run the jl script in another process because then the container is reloaded
# so the ssh setup disappears
cat <<EOF > temp.jl
    using Pkg
    Pkg.activate("$julia_env")

    # clone each private repo
    tmp_folder = mktempdir()
    run(
        addenv(
            `git clone git@github.com:Julia-Tempering/autoHMC.git \$tmp_folder/autoHMC`, 
            "GIT_SSH_COMMAND" => "ssh -i keys/id_autoHMC -o IdentitiesOnly=yes"
        )
    )
    run(
        addenv(
            `git clone git@github.com:Julia-Tempering/autoRWMH.git \$tmp_folder/autoRWMH`, 
            "GIT_SSH_COMMAND" => "ssh -i keys/id_autoRWMH -o IdentitiesOnly=yes"
        )
    )

    # need to install them jointly otherwise Julia complains about unregistered pkgs
    Pkg.add([ 
        Pkg.PackageSpec(name="Pigeons", rev="main"),
        Pkg.PackageSpec(path=joinpath(tmp_folder, "autoHMC")),
        Pkg.PackageSpec(path=joinpath(tmp_folder, "autoRWMH"))
    ])
    Pkg.instantiate()
    Pkg.precompile()

    # force download of stan by loading BridgeStan here
    # no other package uses this approach
    # also, need to precompile stan models to avoid race conditions since the
    # .o and .so files are stored in a shared location (pigeons examples folder)
    using BridgeStan, Pigeons
    Pigeons.toy_stan_target(1)
    Pigeons.stan_funnel(1)
    Pigeons.stan_banana(1)
    StanLogPotential(
        joinpath("$baseDir", "stan", "horseshoe_logit.stan"), 
        Pigeons.json(; n=0,d=0,x=zeros((0,0)),y=zeros(0))
    )
    StanLogPotential(
        joinpath("$baseDir", "stan", "two_component_normal.stan"), 
        Pigeons.json(; n=1, s_lo=0.1, s_hi=10.0)
    )
EOF

# Run the Julia script
julia temp.jl

