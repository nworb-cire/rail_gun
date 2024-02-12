using DifferentialEquations
using Plots
using Unitful: Quantity, m, cm, s, ms, A, V, Ω, g, N, J, F, @u_str, ustrip
using Distributions

include("solver.jl")


const μ = 4e-7π*N/A^2
const grav = 9.81m/s^2


default_params = (
    M=50g, 
    R=0.3Ω, 
    C=.0011F,
    c=0.2,  # coefficient of friction, tungsten on tungsten
    n=1, #number of turns
)

const rail_length = 0.5m

function eq!(du, u, p, t)
    d, d′, I, Eₛ, _ = u

    Fₙ = 15N #I don't know, but I suspect this is enough (about 5 newtons per pound)
    Fₛ = p.c * Fₙ * sign(d′)
    d″ = (I^2 * μ * p.n / (π * p.M)) - (Fₛ / p.M)
    I′ = ((π / (p.n*μ)) * (Eₛ - p.R*I) - d′*I) / d
    Eᵦ = (μ*p.n / π)*(d′*I + I′*d)
    E = Eᵦ-Eₛ
    Eₛ′ = E/(p.R*p.C)
    P = I * E
    du .= [d′, d″, I′, Eₛ′, P]
    return nothing
end

u₀ = [
    1cm,
    0.0m/s,
    0.0A,
    800.0V,
    0.0J,
]

# prob = ODEProblem(eq!, u₀, (0.0s, 500ms), default_params)

# callback: stop when distance is equal to rail length
cb = ContinuousCallback((u, t, i) -> ustrip(rail_length - u[1]), terminate!)
# sol = solve(prob, Tsit5() ; callback=cb)


using DataFrames
df = DataFrame()

Eₛ = parse(Float64, ARGS[1])*V
C = parse(Float64, ARGS[2])*F
Rᵢ = parse(Float64, ARGS[3])*Ω
for M ∈ 4g:1g:6g, n ∈ 3:1:10
    R = Rᵢ + .00084Ω*n 
    u₀ = [
        1cm,
        0.0m/s,
        0.0A,
        Eₛ,
        0.0J,
    ]
    prob = ODEProblem(eq!, u₀, (0.0s, 500ms), default_params)
    local sol = solve(remake(prob; p=(; M, R, C, default_params.c, n)), Tsit5() ; callback=cb)
    efficiency = (M*sol.u[end][2])/(C*(Eₛ^2 -sol.u[end][4]^2))
    println(n, ", ",M, ", ", R, ", ", Eₛ, ", ", sol.u[end][2],
            ", ", maximum(sol[3,:]), ", ", efficiency, "%, ",
            sol.u[end][4], ", ", sol.retcode)
    push!(df, (M=M, R=R, Eₛ=Eₛ, C=C, turns=n, t=sol.t[end], v=sol.u[end][2], amps=sol.u[end][3], voltage=sol.u[end][4], power=sol.u[end][5], retcode=sol.retcode))
end

using CSV
CSV.write("rail_gun.csv", df)
