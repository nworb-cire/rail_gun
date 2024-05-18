using DifferentialEquations
using Plots
using Unitful: Quantity, m, s, ms, A, V, Ω, g, N, J, F, @u_str, ustrip, uconvert
using Distributions
using JSON

include("solver.jl")

const μ₀ = 4e-7π*N/A^2

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
    d, d′, I, Vₛ = u
    r = p.Wᵣ + p.rₚ # r is the radius from the center of the bore to the edge of the rail
    Y = 2*r + 2*p.n*p.Wᵥ # Y is the total width
    X = rail_length
    Alist = [r+n*p.Wᵥ for n=1:p.n]
    M = sum([ ((X-d)/(A^2+(X-d)^2)^0.5 + d/(A^2+d^2)^0.5)/A for A in Alist ])
    M += d / (r^2+d^2)^0.5 / r
    Fₛ = p.c * p.Fₙ * sign(d′)
    d″ = (I^2*p.rₚ*p.μₚ*M / (π * p.M)) - (Fₛ / p.M)

    D = lₚ/(lₚ^2 + r^2)^0.5/r
    E = ((X-d)/((d-X)^2+r^2)^0.5 + d/(d^2+r^2)^0.5)/r
    F = ((r^2+(d+lₚ)^2)^0.5-(r^2+lₚ^2)^0.5+r-(r^2+d^2)^0.5)/r
    G = ((r^2+X^2)^0.5-(r^2+(d-X)^2)^0.5-r+(r^2+d^2)^0.5)/r
    H = sum([ ((A^2+(d+lₚ)^2)^0.5 - (A^2 +(X-d-lₚ)^2)^0.5+(A^2+(X-d)^2)^0.5-(A^2+d^2)^0.5)/A for A in Alist])
    J = sum([ 2*((A^2+X^2)^0.5-A)/A for A in Alist])

    #I′ = ((π / μₚ) * (Vₛ - p.R*I) - d′*I) / (p.n*X + d)
    μ₋ = p.μₚ - μ₀
    I′ = (2*π*Vₛ - 2*π*I*p.R - I*(μ₀*Y*E*d′+μ₋*2*p.rₚ*D*d′))/(μ₀*Y*(J+G)+μ₋*2*p.rₚ*(H+F))
    Vₛ′ = -I/p.C
    du .= [d′, d″, I′, Vₛ′]
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
Wᵥ = params["wireWidth"]m
Wₕ = params["wireHeight"]m
Aₗ = Wᵥ*Wₕ
rₚ = params["projectileDiameter"]m/2
Wᵣ = params["railWidth"]m
lₚ = params["projectileLength"]m
Dₚ = params["projectileDensity"]g/m^3
vᵢ = params["initialVelocity"]m/s
μᵣ = params["relativePermeability"]
μₚ = μᵣ * μ₀
M = Dₚ*lₚ*π*rₚ^2
ρ = 1.77 * 10^-8*Ω*m # resistivity of copper
Rₘ = ρ*n*rail_length/Aₗ
println(Rₘ, " vs ", Rᵢ)
R = Rᵢ + Rₘ 
u₀ = [
    0m,
    vᵢ,
    0.0A,
    Eₛ
]
prob = ODEProblem(eq!, u₀, (0.0s, 500ms), default_params)
sol = solve(remake(prob; p=(; M, R, C, c, n, Fₙ, rₚ, Wᵣ, Wᵥ, μₚ)), Tsit5() ; callback=cb)
Eₚ = uconvert(J, 0.5 * M*(sol.u[end][2]^2 - vᵢ^2))
Eᵦ = uconvert(J, 0.5 * C * (Eₛ^2 -sol.u[end][4]^2))
println(Eₚ, Eᵦ)
efficiency = Eₚ/Eᵦ
t = argmax(sol[3,:])
du = [0.0m/s,0.0m/(s^2),0.0A/s,0.0V/s]
eq!(du, sol.u[1], (; M, R, C, c, n, Fₙ, rₚ, Wᵣ, Wᵥ, μₚ), 0)
println(n, ", ",M, ", ", R, ", ", Eₛ, ", ", sol.u[end][2],
        ", ", sol[3,t],", ", du[3], ", ", efficiency, "%, ",
        sol.u[end][4], ", ", sol.retcode)
push!(df, (M=M, R=R, Eₛ=Eₛ, C=C, turns=n, t=sol.t[end], v=sol.u[end][2], amps=sol.u[end][3], voltage=sol.u[end][4], retcode=sol.retcode))
#println(sol.u[1], "*******", sol[1], "\n\n******",sol[1,:])
plot(sol.t,sol[1,:])
png("distanceByTime.png")
plot(sol.t,sol[2,:])
png("velocityByTime.png")
plot(sol.t,sol[3,:])
png("ampsByTime.png")
plot(sol.t,sol[4,:])
png("voltageByTime.png")