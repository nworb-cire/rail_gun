using DifferentialEquations
using Plots
using Unitful, Unitful.DefaultSymbols


const μ = 4e-7π*N/A^2


p = (M=3g, R=1e-5Ω, Eₛ=100V, k=1e-2m)


function eq!(du, u, p, t)
    d, d′, I = u

    d″ = I^2 * μ / (π * p.M)
    I′ = ((π / μ) * (p.Eₛ - p.R*I) - d′*I) / (d + p.k)
    @show du .= [d′, d″, I′]
end

u₀ = [
    0.0m,
    1e-2m/s,
    0.0A,
]

prob = ODEProblem(eq!, u₀, (0.0s, 4e-4s), p)
sol = solve(prob, Tsit5())

rail_length = 0.3m

# callback: stop when distance is equal to rail length
condition(u, t, integrator) = @show ustrip(rail_length - u[1])
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)

import Base.log2
Base.log2(x::Unitful.Quantity) = typeof(x)(Base.log2(ustrip(x)))
import Base.ceil
Base.ceil(x::Unitful.Quantity) = typeof(x)(Base.ceil(ustrip(x)))

sol = solve(prob, Tsit5(), callback=cb)