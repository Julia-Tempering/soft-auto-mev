include("instantiate_env.jl")
include("utils.jl")

using FiniteDiff
using LogDensityProblems
using LinearAlgebra

function get_reference_samples(target)
    pt = pigeons(
        target = target,
        seed = 1,
        explorer = HitAndRunSlicer(),
        n_chains = 1,
        n_rounds = 15,
        record = [traces]
    )
    sa = sample_array(pt)[:,1:5,1]
    collect(copy(c) for c in eachrow(sa)) # use vector of vectors representation for faster iteration
end
function make_stan_grad(target)
    buff_ad = Pigeons.BufferedAD(target, Pigeons.buffers(), Ref(0.0), Ref{Cstring}())
    (x) -> last(LogDensityProblems.logdensity_and_gradient(buff_ad, x))
end
make_fd_grad(target::StanLogPotential) =
    (x) -> FiniteDiff.finite_difference_gradient(
        Base.Fix1(LogDensityProblems.logdensity, target), x, Val{:central}
    )
function test()
    target = stan_logpotential("mRNA")
    fd_grad = make_fd_grad(target)
    stan_grad = make_stan_grad(target)
    ch = get_reference_samples(target)
    mean(ch) do x
        norm(fd_grad(x)-stan_grad(x))/norm(fd_grad(x))
    end
end

# run
test()
# result = 4.0421319771296056e-10
