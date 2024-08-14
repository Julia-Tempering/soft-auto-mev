include("instantiate_env.jl")
include("utils.jl")

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

rng = SplittableRandom(1)
a = StandardAR1Process(0.9)
spectrum_at_zero(a)
n = 9
res = simulate(a, rng, n)
var(res)
asymptotic_var(a)
v = collect(res)
batch_means_ess(v, 0.0, sqrt(asymptotic_var(a)))/n
batch_means_ess(v)/n

k = 3
s = Iterators.Stateful(simulate(a, SplittableRandom(1), n))
collect(Iterators.take(s, k))