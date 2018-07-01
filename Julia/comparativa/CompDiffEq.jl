
### Requiere crear la función r con Func_r del archivo Solve2body.jl
###    using DiferentialEquations
###
include("Solve2body.jl")

#r,T=Func_r(0.0)
r=1.0
T=2*pi
m=[.5,.5];
s=[-.5,.5];

function Fuerza3c(du,u,p,t)
    m=p[1]
    s=p[2]
    r=p[3]
    R = (r^2*s.^2 +u[1]^2).^(1.5)
    #R = (r(t)^2*s.^2 +u[1]^2).^(1.5)
    du[1] = u[2]
    du[2] = -sum(m./R*u[1])
end
#z₀=1
u0 = [1,0.0]
tspan = (0.0,1000*T)
prob_ode_trescuerpos = ODEProblem(Fuerza3c,u0,tspan,[m,s,r]);
alg =   VCAB3()
sol = solve(prob_ode_trescuerpos,alg,dt=1e-3,progress=true, abstol = 1e-12, reltol = 1e-12);

t=0:.1:1000*T;z=sol.(t);Data=Array{Float64}(length(t),3);Data[:,1]=[a[1] for a in z];Data[:,2]=[a[2] for a in z];Data[:,3]=t;writedlm("ADM54-dt1e-3-reltol-abstol-1e-12.txt",Data,',')

