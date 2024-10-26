using Pkg
Pkg.activate(".")
Pkg.add([ # need to install them jointly otherwise Julia complains AutoStep are not registered
    Pkg.PackageSpec(name="Pigeons", rev="main"),
    Pkg.PackageSpec(url="git@github.com:Julia-Tempering/AutoStep.git", rev="main")
])
Pkg.update()
Pkg.instantiate()
Pkg.precompile()
