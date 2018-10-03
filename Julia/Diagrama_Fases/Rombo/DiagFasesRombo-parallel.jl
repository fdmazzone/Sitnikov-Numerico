@everywhere using DifferentialEquations, JLD, Interpolations, ProgressMeter
include("Solve2body.jl")

@everywhere λ=1.3
@everywhere e=.1

v0_lista=[0.0]
z0_lista=linspace(5,5.1,5)

@everywhere CantPer=10;
nombre_arch="Borrame.jld"

size_z0s=length(z0_lista)
size_v0s=length(v0_lista)

@everywhere m₂=4*(1+λ^2)^(3/2)*(8*λ^3-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
@everywhere m₁=4*λ^3*(1+λ^2)^(3/2)*(8-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
@everywhere m=[m₁,m₂]
@everywhere s=[λ 1]

@everywhere r,T=Func_r(e);
@everywhere tspan = (0.0,CantPer*T);
@everywhere alg =  Vern9()
@everywhere t=0:T:CantPer*T;
@everywhere t₀=0.0


@everywhere function Fuerza3c(du,u,p,t)
    m=p[1]
    s=p[2]
    r=p[3]
    R = (r(t)^2*s.^2 +u[1]^2).^(1.5)
    du[1] = u[2]
    du[2] = -2*sum(m./R*u[1])
end



#z=SharedArray{Float64}(CantPer+1,size_z0s,size_v0s)
#v=SharedArray{Float64}(CantPer+1,size_z0s,size_v0s)
#Cont=1

#ProBar= Progress(size_z0s*size_v0s, 1)

@everywhere function DiagFases(u0)
    prob_ode_trescuerpos = ODEProblem(Fuerza3c,u0,tspan,[m,s,r]);
    sol = solve(prob_ode_trescuerpos,alg,dt=1e-3, progress=true,abstol = 1e-12, reltol = 1e-12);
    return sol.(t)
end

u0=[[z,v] for z in z0_lista for v in v0_lista]


c=pmap(DiagFases,u0)




#save(nombre_arch,"z",z,"v",v,"e",e,"z₀",z0_lista,"T",T)
