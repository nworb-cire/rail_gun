using DifferentialEquations
using Plots
using Unitful: Quantity, m, s, ms, A, V, Ω, g, N, J, F, @u_str, ustrip, uconvert
using Distributions
using JSON

include("solver.jl")

const μ = 4e-7π*N/A^2

default_params = (
    M=50g, 
    R=0.3Ω, 
    C=.0011F,
    c=0.2,  # coefficient of friction, tungsten on tungsten
    n=1, #number of turns
)

params = JSON.parsefile(ARGS[1])
const rail_length = params["railLength"] * m

function eq!(du, u, p, t)
    d, d′, I, Eₛ, _ = u

    Fₛ = p.c * p.Fₙ * sign(d′)
    d″ = (I^2 * μ * p.n / (π * p.M)) - (Fₛ / p.M)
    I′ = ((π / μ) * (Eₛ - p.R*I) - d′*I) / (p.n*rail_length + d)
    Eᵦ = (μ*p.n / π)*(d′*I + I′*d)
    E = Eᵦ-Eₛ
    Eₛ′ = E/(p.R*p.C)
    P = I * E
    du .= [d′, d″, I′, Eₛ′, P]
    return nothing
end

# callback: stop when distance is equal to rail length
cb = ContinuousCallback((u, t, i) -> ustrip(rail_length - u[1]), terminate!)

using DataFrames
df = DataFrame()
nₚ = params["numCapacitors"]
Eₛ = params["capacitorVoltage"]*V
C = params["capacitance"]*F * nₚ
Rᵢ = params["internalResistance"]*Ω/nₚ
c = params["coefficientOfFriction"]
Fₙ = params["contactPressure"]*N
n = params["numberOfCoils"]
Aₗ = params["wireArea"]m^2
rₚ = params["projectileDiameter"]m/2
lₚ = params["projectileLength"]m
Dₚ = params["projectileDensity"]g/m^3
M = Dₚ*lₚ*π*rₚ^2
ρ = 1.77 * 10^-8*Ω*m # resistivity of copper
Rₘ = ρ*n*rail_length/Aₗ
println(Rₘ, " vs ", Rᵢ)
R = Rᵢ + Rₘ 
u₀ = [
    0m,
    0.0m/s,
    0.0A,
    Eₛ,
    0.0J,
]
prob = ODEProblem(eq!, u₀, (0.0s, 500ms), default_params)
sol = solve(remake(prob; p=(; M, R, C, c, n, Fₙ)), Tsit5() ; callback=cb)
Eₚ = uconvert(J, 0.5 * M*(sol.u[end][2]^2))
Eᵦ = uconvert(J, 0.5 * C * (Eₛ^2 -sol.u[end][4]^2))
println(Eₚ, Eᵦ)
efficiency = Eₚ/Eᵦ
t = argmax(sol[3,:])
du = [0.0m/s,0.0m/(s^2),0.0A/s,0.0V/s,0.0]
eq!(du, sol.u[1], (; M, R, C, c, n, Fₙ), 0)
println(n, ", ",M, ", ", R, ", ", Eₛ, ", ", sol.u[end][2],
        ", ", sol[3,t],", ", du[3], ", ", efficiency, "%, ",
        sol.u[end][4], ", ", sol.retcode)
push!(df, (M=M, R=R, Eₛ=Eₛ, C=C, turns=n, t=sol.t[end], v=sol.u[end][2], amps=sol.u[end][3], voltage=sol.u[end][4], power=sol.u[end][5], retcode=sol.retcode))
#println(sol.u[1], "*******", sol[1], "\n\n******",sol[1,:])
plot(sol.t,sol[1,:])
png("distanceByTime.png")
plot(sol.t,sol[2,:])
png("velocityByTime.png")
plot(sol.t,sol[3,:])
png("ampsByTime.png")
plot(sol.t,sol[4,:])
png("voltageByTime.png")