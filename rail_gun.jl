using DifferentialEquations
using Plots
using Unitful: Quantity, m, cm, s, ms, A, V, Ω, g, N, J
using Distributions

include("solver.jl")


const μ = 4e-7π*N/A^2
const grav = 9.81m/s^2


default_params = (
    M=200g, 
    R=1e-5Ω, 
    Eₛ=100V, 
    c=0.2,  # coefficient of friction, tungsten on tungsten
)

const rail_length = 0.5m
const turns = 2  # the number of turns (including the rails) making the magnetic field
@assert turns >= 1  # there must at least the rails generating a field

function eq!(du, u, p, t)
    d, d′, I, _ = u

    Fₙ = p.M * grav
    Fₛ = p.c * Fₙ * sign(d′)
    d″ = (I^2 * μ * turns / (π * p.M)) - (Fₛ / p.M)
    I′ = ((π / (turns*μ)) * (p.Eₛ - p.R*I) - d′*I) / d
    P = I * p.Eₛ
    du .= [d′, d″, I′, P]
    return nothing
end

u₀ = [
    1cm,
    0.0m/s,
    0.0A,
    0.0J,
]

prob = ODEProblem(eq!, u₀, (0.0s, 50ms), default_params)

# callback: stop when distance is equal to rail length
cb = ContinuousCallback((u, t, i) -> ustrip(rail_length - u[1]), terminate!)
sol = solve(prob, Tsit5() ; callback=cb)


using DataFrames
df = DataFrame()

for M ∈ 4g:1g:6g, R ∈ 10.0.^(-3:0.2:0) .*Ω, Eₛ ∈ 200V:50V:1300V
    local sol = solve(remake(prob; p=(; M, R, Eₛ, default_params.c)), Tsit5() ; callback=cb)
    println(M, R, Eₛ, sol.t[end], sol.u[end][2])
    push!(df, (M=M, R=R, Eₛ=Eₛ, t=sol.t[end], v=sol.u[end][2], energy=sol.u[end][4], retcode=sol.retcode))
end

using CSV
CSV.write("rail_gun.csv", df)
