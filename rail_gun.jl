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
rail_length = ustrip(params["railLength"] * m)
df = DataFrame()
nₚ = params["numCapacitors"]
Eₛ = ustrip(params["capacitorVoltage"]*V)
C = ustrip(params["capacitance"]*F * nₚ)
Rᵢ = ustrip(params["internalResistance"]*Ω/nₚ)
c = ustrip(params["coefficientOfFriction"])
Fₙ = ustrip(params["contactPressure"]*N)
n = ustrip(params["numberOfCoils"])
Wᵥ = ustrip(params["wireWidth"]m)
Wₕ = ustrip(params["wireHeight"]m)
Aₗ = Wᵥ*Wₕ
rₚ = ustrip(params["projectileWidth"]m/2)
Hₚ = ustrip(params["projectileHeight"]m)
Wᵣ = ustrip(params["railWidth"]m)
lₚ = ustrip(params["projectileLength"]m)
Dₚ = ustrip(params["projectileDensity"]kg/m^3)
vᵢ = ustrip(params["initialVelocity"]m/s)
μH = ustrip(params["externalMagneticField"]T)
M = Dₚ*lₚ*2*rₚ*Hₚ
ρ = ustrip(1.77 * 10^-8*Ω*m) # resistivity of copper
Rₘ = ρ*(n+1)*rail_length/Aₗ
#println(Rₘ, " vs ", Rᵢ)
R = Rᵢ + Rₘ 

function eq!(du, u, p, t)
    d, d′, I, Vₛ = u
    r = p.Wᵣ + p.rₚ # r is the radius from the center of the bore to the edge of the rail
    Y = 2*r + 2*p.n*p.Wᵥ # Y is the total width
    X = p.rail_length
    Alist = [r+n*p.Wᵥ for n=1:p.n]
    if (p.n == 0)
        M = d / (r^2+d^2)^0.5 / r
        D = ((r^2 + X^2)^0.5 - (r^2 + (d-X)^2)^0.5 -r + (r^2+d^2)^0.5)/r
        F = (2*(r^2+d^2)^0.5 -2*r)/r
        G = 2*d*d′ / (r*(r^2+d^2)^0.5)
    else
        M = sum([ ((X-d)/(A^2+(X-d)^2)^0.5 + d/(A^2+d^2)^0.5)/A for A in Alist ])
        M += d / (r^2+d^2)^0.5 / r
        D = sum([ (2*(A^2 + X^2)^0.5 - 2*A)/A for A in Alist]) + ((r^2 + X^2)^0.5 - (r^2 + (d-X)^2)^0.5 -r + (r^2+d^2)^0.5)/r
        F = sum([ ((A^2+d^2)^0.5 - (A^2 + (X-d)^2)^0.5 - A + (A^2+X^2)^0.5)/A for A in Alist]) + (2*(r^2+d^2)^0.5 -2*r)/r
        G = sum([ (d*d′/(A^2+d^2)^0.5 + (X-d)*d′/(A^2+(X-d)^2)^0.5)/A for A in Alist]) + 2*d*d′ / (r*(r^2+d^2)^0.5)
    end
    E = (d*d′/(r^2+d^2)^0.5 - (d-X)d′/(r^2 + (d-X)^2)^0.5)/r
    #this is a hack to prevent friction from messing things up when velocity is close to zero
    if (d′ < ustrip(0.01m/s))
        Fₛ = ustrip(0N)
    else
        Fₛ = p.c * p.Fₙ * sign(d′)
    end
    d″ = (I^2*p.rₚ*μ₀*M / (π * p.M)) + (I*2*p.rₚ*μH/p.M) - (Fₛ / p.M)


    #I′ = ((π / μₚ) * (Vₛ - p.R*I) - d′*I) / (p.n*X + d)
    Vᵦ = Vₛ - I*p.R
    I′ = (Vᵦ - I*((Y*μ₀)/(2*π))*(p.n*E + G) - μH*Y*d′)/(μ₀*Y/(2*π)*(p.n*D + F) )
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
        ustrip(0.01m),
        vᵢ,
        ustrip(0.0A),
        Emf_source
    ]
    Aₗ = Wᵥ*Wₕ
    Rᵣ = 2*ρ * rail_length/(Wᵣ*ustrip(0.02m))
    Rₘ = 2*ρ*(n)*rail_length/Aₗ
    R = Rᵢ + Rₘ + Rᵣ
    #things that can change M; antecedants to R; r_p, Wᵣ, Wᵥ
    prob = ODEProblem(eq!, u₀, (0.0, 0.5), default_params)
    sol = solve(remake(prob; p=(; M, R, C=C_new, c, n, Fₙ, rₚ, Wᵣ, Wᵥ, rail_length)), Tsit5() ; callback=cb)
    Eₚ = 0.5 * M*(sol.u[end][2]^2 - vᵢ^2)
    Eᵦ = 0.5 * C_new * (Emf_source^2 -sol.u[end][4]^2)
    println(Eₚ," ", Eᵦ)
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