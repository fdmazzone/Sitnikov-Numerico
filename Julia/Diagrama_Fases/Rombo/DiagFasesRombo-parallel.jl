@everywhere using DifferentialEquations, JLD, Interpolations, ProgressMeter
@everywhere include("Solve2body.jl")

@everywhere λ=1.3
@everywhere e=.1

v0_lista=[0.0]
z0_lista=linspace(5,5.1,5)

@everywhere CantPer=10;
nombre_arch="Borrame.jld"


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




@everywhere function DiagFases(U0,p=nothing)
    ZV=Array{Array{Float64,1},1}
    @parallelprogress p string("Procesando...") for i in 1:length(U0) 
        prob_ode_trescuerpos = ODEProblem(Fuerza3c,U0[i],tspan,[m,s,r]);
        sol = solve(prob_ode_trescuerpos,alg,dt=1e-3, progress=true,abstol = 1e-12, reltol = 1e-12);
        ZV[i]=sol.(t)
    end
    return ZV
end

u0=[[z,v] for z in z0_lista for v in v0_lista]


c=ProgressMeter.pmap(DiagFases,[u0[1:2],u0[3:end]])




#save(nombre_arch,"z",z,"v",v,"e",e,"z₀",z0_lista,"T",T)
