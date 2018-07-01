#include("Solve2body.jl")
#r,T=Func_r(0.0)
include("/home/fernando/fer/Investigación/Trabajo en curso/Mecanica Celeste/Julia/metodo_colocacion/collocation.jl")
#using Interpolations, MAT, ProgressMeter
r=1.0
m=[.5,.5];
s=[-.5 .5];

function ForceSit(t,z,r)

    return -(z.^2.+r.^2*s.^2 ).^(-1.5)*m.*z
    #return -(z.^2.+r.(t).^2*s.^2 ).^(-1.5)*m.*z
end

z₀=Array{Float64,2}(1,1)
z₀[1,1] = 1.0
v₀=Array{Float64,2}(1,1)
v₀[1,1] = 0.0
t₀=0.0
t = collect(0.0:.1:1000*T)

z,v=collocation(ForceSit,r,t,z₀,v₀,t₀;dt=.001,order=8,tol=1e-16,dtmin=1e-8,dtmax=1);




Data=Array{Float64}(length(t),3);Data[:,1]=z;Data[:,2]=v;Data[:,3]=t;writedlm("Colloc-dt1e-3-or-8-tol-1e-16.txt",Data,',')
