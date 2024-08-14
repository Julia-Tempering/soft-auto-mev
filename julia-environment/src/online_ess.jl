include("batch_means.jl")
using OnlineStatsBase
using Random

#=
An online exact (i.e., known mean and sd) ESS estimator based on truncated
spectral density at zero (i.e., truncated integrated autocorrelation time).

TODO: results are too dependent on the autocorrelation index cutoff K.
=#
struct ExactESS{TF<:Real,TCB<:CircBuff,TACS} <: OnlineStat{TF}
    μ::TF
    σ::TF
    buffer::TCB
    acs::TACS
    function ExactESS(μ::TF, σ::TF, buffer::CircBuff{TF}, acs::NTuple) where {TF}
        @assert length(buffer.rng.rng) == length(acs)
        new{TF,typeof(buffer),typeof(acs)}(μ, σ, buffer, acs)
    end
end
function ExactESS(μ::Real, σ::Real; K::Integer = 256)
    TF = float(promote_type(typeof(μ), typeof(σ)))
    buffer = CircBuff(TF, K)
    acs = ntuple(i -> Mean(TF), K)
    ExactESS(TF(μ), TF(σ), buffer, acs)
end
OnlineStatsBase.nobs(e::ExactESS) = nobs(e.buffer)
std_spectrum_at_zero(e::ExactESS) = max(zero(e.μ), one(e.σ) + 2sum(value(e.buffer)))
spectrum_at_zero(e::ExactESS) = abs2(e.σ) * std_spectrum_at_zero(e)
relative_ESS(e::ExactESS) = inv(std_spectrum_at_zero(e))
OnlineStatsBase.value(e::ExactESS) = nobs(e) * relative_ESS(e)
function OnlineStatsBase._fit!(e::ExactESS, y)
    z = (y-e.μ)/e.σ                     # standardize observation
    k = length(e.buffer)                # number of available elements in the buffer
    for i in 1:k
        fit!(e.acs[i], e.buffer[i] * z) # update the i-th autocorrelation estimator
    end
    fit!(e.buffer, z)                   # add standardized obs to the circular buffer
end


# test
vals = randn(100000)
# vals = range(0.0,5.0,length=100000)
e = ExactESS(0, 1; K=4)
fit!(e, vals)
value.(e.acs)
std_spectrum_at_zero(e)
batch_means_ess(vals,0.0,1.0)


#=
A standardized (σ=1) AR1 process for testing the estimator
=#

struct StandardAR1Process
    ϕ::Float64
end
spectrum_at_zero(a::StandardAR1Process) = inv(abs2(one(a.ϕ)-a.ϕ))
asymptotic_var(a::StandardAR1Process) = inv(one(a.ϕ)-abs2(a.ϕ))
relative_ESS(a::StandardAR1Process) = (one(a.ϕ)-a.ϕ)/(one(a.ϕ)+a.ϕ) # numerically nicer but equivalent to asymptotic_var(a)/spectrum_at_zero(a)
function simulate(a::StandardAR1Process, rng::AbstractRNG, n::Integer)
    x0 = sqrt(asymptotic_var(a)) * randn(rng)
    tail = Iterators.accumulate((x,ϵ) -> a.ϕ*x + ϵ, (randn(rng) for _ in 2:n); init=x0)
    Iterators.flatten((x0, tail))
end

rng = Xoshiro(1)
a = StandardAR1Process(0.9)
spectrum_at_zero(a)
n = 100_000
res = simulate(a, rng, n)
e = ExactESS(0, sqrt(asymptotic_var(a)); K=63)
nobs(e)
fit!(e,res)
relative_ESS(e)
relative_ESS(a)
value.(e.acs)