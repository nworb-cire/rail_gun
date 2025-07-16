using DifferentialEquations
using Plots
using Unitful: Quantity, m, s, ms, A, V, Ω, kg, N, J, F, T, @u_str, ustrip, uconvert
using Distributions
using JSON
using DataFrames
using Optim
include("solver.jl")

μ₀ = ustrip(4e-7π*N/A^2)

default_params = (
    M=ustrip(0.05kg), 
    R=ustrip(0.3Ω), 
    C=ustrip(.0011F),
    c=0.2,  # coefficient of friction, tungsten on tungsten
    n=1, #number of turns
)

params = JSON.parsefile(ARGS[1])
beg = ustrip(params["extraCircuitDistance"]m)
rail_length = ustrip(params["railLength"] * m) + beg
df = DataFrame()
nₚ = params["numCapacitors"]
Eₛ = ustrip(params["capacitorVoltage"]*V)
C = ustrip(params["capacitance"]*F * nₚ)
Rᵢ = ustrip(params["internalResistance"]*Ω/nₚ)
c = ustrip(params["coefficientOfFriction"])
Fₙ = ustrip(params["contactPressure"]*N)
n = ustrip(params["numberOfCoils"])
rₚ = ustrip(params["projectileWidth"]m/2)
Hₚ = ustrip(params["projectileHeight"]m)
Wᵣ = ustrip(params["railWidth"]m)
Hᵣ = ustrip(params["railHeight"]m)
lₚ = ustrip(params["projectileLength"]m)
Dₚ = ustrip(params["projectileDensity"]kg/m^3)
vᵢ = ustrip(params["initialVelocity"]m/s)
μH = ustrip(params["externalMagneticField"]T)
M = Dₚ*lₚ*2*rₚ*Hₚ
ρ = ustrip(1.68 * 10^-8*Ω*m) # resistivity of copper
ρ = ustrip(2.82 * 10^-8*Ω*m) # resistivity of aluminum
ρ = ustrip(5.6 * 10^-8*Ω*m) # resistivity of tungsten
ρ = ustrip(5.9 * 10^-8*Ω*m) # resistivity of zink
ρ = ustrip(4.20 * 10^-7*Ω*m) # resistivity of titanium
ρ = ustrip(6.9 * 10^-7*Ω*m) # resistivity of stainless steel
R = Rᵢ 

function eq!(du, u, p, t)
    d, d′, I, Vₛ = u
    Aᵪ = Wᵣ*Hᵣ # cross sectional area of one rail
    R = ρ*p.n*d/Aᵪ # resistance in the rails
    R += p.R # resistance in the whole system
    r = p.Wᵣ/2 + p.rₚ # r is the radius from the center of the bore to the center of the rail
    Y = 2*p.rₚ # Y is the projectile width
    M = p.n*d / (r^2+d^2)^0.5 / r
    #this is a hack to prevent friction from messing things up when velocity is close to zero
    if (d′ < ustrip(0.01m/s))
        Fₛ = ustrip(0N)
    else
        Fₛ = p.c * p.Fₙ * sign(d′)
    end
    d″ = (p.n*I^2*p.rₚ*μ₀*M / (π * p.M)) + (p.n*I*2*p.rₚ*μH/p.M) - (Fₛ / p.M)

    Vᵦ = Vₛ - I*R
    I′ = (Vᵦ - I*((Y*p.n^2*μ₀)/(π))*(d*d′/(r^2+d^2)^0.5/r) - μH*Y*d′)/(μ₀*Y*p.n^2/(π)*((r^2+ d^2)^0.5-r) )
    Vₛ′ = -I/p.C
    du .= [d′, d″, I′, Vₛ′]
    return nothing
end

# callback: stop when distance is equal to rail length
cb = ContinuousCallback((u, t, i) -> ustrip(rail_length - u[1]), terminate!)
function toDifferentiate(vec)
    efficiency, sol = runSim(vec)
    return efficiency
end
function runSim(vec)
    Emf_source = Eₛ
    M = max(ustrip(.001kg),vec[1])
    rail_length = max(ustrip(0.09m), min(ustrip(1m), vec[2]))
    C_new = C
    u₀ = [
        beg,
        vᵢ,
        ustrip(0.0A),
        Emf_source
    ]
    #things that can change M; antecedants to R; r_p, Wᵣ, Wᵥ
    prob = ODEProblem(eq!, u₀, (0.0, 0.5), default_params)
    sol = solve(remake(prob; p=(; M, R=Rᵢ, C=C_new, c, n, Fₙ, rₚ, Wᵣ, rail_length)), Tsit5() ; callback=cb)
    Eₚ = 0.5 * M*(sol.u[end][2]^2 - vᵢ^2)
    Eᵦ = 0.5 * C_new * (Emf_source^2 -sol.u[end][4]^2)
    println(Eₚ," ", Eᵦ)
    println(sol.u[end][2])
    efficiency = Eₚ/Eᵦ
    return -efficiency, sol
end

initialVector = [
    M,
    rail_length
]
print(initialVector);
out, sol = runSim(initialVector);
println(out);
if (ARGS[2] != "skip")
    optimalVector = Optim.minimizer(optimize(toDifferentiate, initialVector, BFGS(); autodiff = :forward))
    effic, sol = runSim(optimalVector);
    println("\n\n And we're done!\n");
    println(sol);
    println(optimalVector)
end
"""
u₀ = [
    0.01m,
    vᵢ,
    0.0A,
    Eₛ
]
prob = ODEProblem(eq!, u₀, (0.0s, 500ms), default_params)
sol = solve(remake(prob; p=(; M, R, C, c, n, Fₙ, rₚ, Wᵣ, Wᵥ)), Tsit5() ; callback=cb)
Eₚ = uconvert(J, 0.5 * M*(sol.u[end][2]^2 - vᵢ^2))
Eᵦ = uconvert(J, 0.5 * C * (Eₛ^2 -sol.u[end][4]^2))
println(Eₚ, Eᵦ)
efficiency = Eₚ/Eᵦ
t = argmax(sol[3,:])
du = [0.0m/s,0.0m/(s^2),0.0A/s,0.0V/s]
eq!(du, sol.u[1], (; M, R, C, c, n, Fₙ, rₚ, Wᵣ, Wᵥ), 0)
println(n, ", ",M, ", ", R, ", ", Eₛ, ", ", sol.u[end][2],
        ", ", sol[3,t],", ", du[3], ", ", efficiency, "%, ",
        sol.u[end][4], ", ", sol.retcode)
push!(df, (M=M, R=R, Eₛ=Eₛ, C=C, turns=n, t=sol.t[end], v=sol.u[end][2], amps=sol.u[end][3], voltage=sol.u[end][4], retcode=sol.retcode))
#println(sol.u[1], "*******", sol[1], "\n\n******",sol[1,:])
"""

plot(sol.t,sol[1,:])
png("distanceByTime.png")
plot(sol.t,sol[2,:])
png("velocityByTime.png")
plot(sol.t,sol[3,:])
png("ampsByTime.png")
plot(sol.t,sol[4,:])
png("voltageByTime.png")