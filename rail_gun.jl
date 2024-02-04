using DifferentialEquations
using Plots
using Unitful, Unitful.DefaultSymbols
using Distributions

include("solver.jl")


const μ = 4e-7π*N/A^2


p = (
    M=200g, 
    R=1e-5Ω, 
    Eₛ=100V, 
)

const rail_length = 0.3m

function eq!(du, u, p, t)
    d, d′, I = u

    d″ = I^2 * μ / (π * p.M)
    I′ = ((π / μ) * (p.Eₛ - p.R*I) - d′*I) / d
    du .= [d′, d″, I′]
    return nothing
end

u₀ = [
    1cm,
    0.0m/s,
    0.0A,
]

prob = ODEProblem(eq!, u₀, (0.0s, 1ms), p)

# callback: stop when distance is equal to rail length
cb = ContinuousCallback((u, t, i) -> ustrip(rail_length - u[1]), terminate!)
sol = solve(prob, Tsit5() ; callback=cb)


using DataFrames
df = DataFrame()

for M ∈ 3g:10g:200g, R ∈ 1e-5Ω:1e-5Ω:1e-4Ω, Eₛ ∈ 100V:20V:800V
    sol = solve(remake(prob; p=(; M, R, Eₛ)), Tsit5() ; callback=cb)
    println(M, R, Eₛ, sol.t[end], sol.u[end][2])
    push!(df, (M=M, R=R, Eₛ=Eₛ, t=sol.t[end], v=sol.u[end][2]))
end

using CSV
CSV.write("rail_gun.csv", df)
