include("batch_means.jl")
using OnlineStatsBase
using Random

#=
An online exact (i.e., known mean and sd) ESS estimator based on truncated
spectral density at zero (i.e., truncated integrated autocorrelation time).
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
function ExactESS(μ::Real, σ::Real; K::Integer = 64)
    TF = float(promote_type(typeof(μ), typeof(σ)))
    buffer = CircBuff(TF, K, rev=true) # buffer[1] = most recently added element => used for ac(1)
    acs = ntuple(i -> Mean(TF), K)
    ExactESS(TF(μ), TF(σ), buffer, acs)
end
OnlineStatsBase.nobs(e::ExactESS) = nobs(e.buffer)
std_spectrum_at_zero(e::ExactESS) = max(zero(e.μ), one(e.σ) + 2sum(value, e.acs)) # impose constraint that spectral densities are positive functions
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


#=
Tests
=#

using CairoMakie, AlgebraOfGraphics, DataFrames

n_sims = 100_000
ϕs = 0.99*(2. .^ range(0,1,length=10) .- 1)
seeds = 1:30
Ks = 2 .^ (4:8)
df = mapreduce(vcat, ϕs) do ϕ      # loop ϕ
    a = StandardAR1Process(ϕ)
    σ = sqrt(asymptotic_var(a))
    tru = spectrum_at_zero(a)
    mapreduce(vcat, seeds) do seed # loop seed
        rng = Xoshiro(seed)
        res = simulate(a, rng, n_sims)
        mapreduce(vcat, Ks) do K   # loop K
            e = ExactESS(0, σ; K=K)
            fit!(e, res)
            DataFrame(phi = ϕ, seed=seed, K=K, rerr = abs(spectrum_at_zero(e) - tru)/tru)
        end
    end
end

# plot
s1 = data(df) * 
    visual(BoxPlot) *
    mapping(
        :phi => (ϕ -> string(round(ϕ,digits=2))) => "Autocorrelation parameter",
        :rerr => log10 => "Relative error of IAT estimator (log10)",
        color=:K => nonnumeric,
        dodge=:K => string
    )
# s2 = data(DataFrame(y=0.05)) * mapping(:y => log10) * visual(HLines)
p=draw(s1, axis = (width = 800, height = 400, title="Experimental accuracy of online exact ESS estimator for AR(ϕ) processes"))
save("temp.png", p, px_per_unit = 3) # save high-resolution png
