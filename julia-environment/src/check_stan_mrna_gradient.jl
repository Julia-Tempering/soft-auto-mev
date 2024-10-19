include("instantiate_env.jl")
include("utils.jl")

using FiniteDiff
using LogDensityProblems
using LinearAlgebra
using Bijectors

# get samples from the prior
function make_prior_dist(;only_marginals=false, constrained=false)
    marginals = [Uniform(-2,1), Uniform(-5,5), Uniform(-5,5), Uniform(-5,5), Uniform(-2,2)]
    return if only_marginals
        marginals
    else
        product_distribution(constrained ? marginals : transformed.(marginals))
    end
end
function get_prior_samples_unconstrained(;n = 2^15, rng = SplittableRandom(1))
    prior_dist_unconstrained = make_prior_dist()
    [rand(rng, prior_dist_unconstrained) for _ in 1:n]
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
    b  = Bijectors.bijector(make_prior_dist(constrained=true))
    sa = sample_array(pt)[:,1:5,1]
    [b(c) for c in eachrow(sa)] # return vector of unconstrained samples
end
function make_stan_grad(target)
    buff_ad = Pigeons.BufferedAD(target, Pigeons.buffers(), Ref(0.0), Ref{Cstring}())
    (x) -> last(LogDensityProblems.logdensity_and_gradient(buff_ad, x))
end
make_fd_grad(target::StanLogPotential) =
    (x) -> FiniteDiff.finite_difference_gradient(
        Base.Fix1(LogDensityProblems.logdensity, target), x, Val{:central}
    )
# compute mean relative error
function test_grad(target; prior=true)
    samples = prior ? get_prior_samples_unconstrained() : get_posterior_samples(target)
    fd_grad = make_fd_grad(target)
    stan_grad = make_stan_grad(target)
    mean(samples) do x
        norm(fd_grad(x)-stan_grad(x))/norm(fd_grad(x))
    end
end

# run
target = stan_logpotential("mRNA")
test_grad(target, prior = false) # posterior => err ~ 3.090724925250216e-7
test_grad(target) # prior => err ~ 4.877355846963999e-6
