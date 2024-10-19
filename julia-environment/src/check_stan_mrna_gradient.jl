include("instantiate_env.jl")
include("utils.jl")

using FiniteDiff
using LogDensityProblems
using LinearAlgebra

# get samples from the prior
function get_prior_samples(;n = 2^15, rng = SplittableRandom(1))
    prior_dist = product_distribution(
        [Uniform(-2,1), Uniform(-5,5), Uniform(-5,5), Uniform(-5,5), Uniform(-2,2)]
    )
    [rand(rng, prior_dist) for _ in 1:n]
end

# get samples from the posterior
function get_posterior_samples(target)
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
function test_grad(target,samples)
    fd_grad = make_fd_grad(target)
    stan_grad = make_stan_grad(target)
    mean(samples) do x
        norm(fd_grad(x)-stan_grad(x))/norm(fd_grad(x))
    end
end

# run
target = stan_logpotential("mrna")
test_grad(target, get_posterior_samples(target)) # 4.0421319771296056e-10
test_grad(target, get_prior_samples()) # 2.1281850912516644e-10
