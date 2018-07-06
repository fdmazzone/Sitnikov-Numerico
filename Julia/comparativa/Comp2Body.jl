
### Requiere crear la función r con Func_r del archivo Solve2body.jl
###    using DiferentialEquations
###
#include("Solve2body.jl")

#r,T=Func_r(0.0)
e=.9
T=2*pi
m=1;
dist_peri=1;
a=dist_peri/(1-e);
n=a^(-1.5);
T=2*pi/n;
v_per=sqrt((1+e)/(1-e)/a);

function Force2b(du,u,p,t)
    m=p[1]
    du[1:2] = u[3:4]
    du[3:4]= -m*norm(u[1:2]).^-3*u[1:2]
end

#z₀=1
u0 = [1,0.0,0.0,v_per]
tspan = (0.0,1000*T)
prob_ode_2cuerpos = ODEProblem(Force2b,u0,tspan,[m]);
alg =   DP8()
sol = solve(prob_ode_2cuerpos,alg,dt=1e-3,progress=true, abstol = 1e-16, reltol = 1e-16);
x1=[a[1] for a in sol.u]
x2=[a[2] for a in sol.u]
v1=[a[3] for a in sol.u]
v2=[a[4] for a in sol.u]
E=.5*(v1.^2+v2.^2)-1./sqrt.(x1.^2+x2.^2)
E=abs.(E-E[1])
t=sol.t/T;Data=Array{Float64}(length(t),2);Data[:,1]=E;Data[:,2]=t;writedlm("DP8-1e-16-2body.txt",Data,',')
