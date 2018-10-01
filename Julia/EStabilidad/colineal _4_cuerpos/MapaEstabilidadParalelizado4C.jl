#cd("/home/fernando/fer/Investigación/Trabajo en curso/Mecanica Celeste/Julia/EStabilidad/colineal #_4_cuerpos")

#procesos=4
#addprocs(procesos)

@everywhere using DifferentialEquations, QuadGK, JLD, ProgressMeter
@everywhere include("MapaEstabilidad4C.jl")




xini=0.38
xfin=0.41
#precisionx=40
precisionx=5

vini=15
vfin=41
#precisionv=2600
precisionv=10


x=linspace(xini,xfin,procesos*precisionx);
x=reshape(x,(precisionx,procesos));
v=linspace(vini,vfin,precisionv);

M=Array{Any}(procesos)
for i in 1:procesos
	M[i]=[x[:,i],v]
end

c=pmap(MapaEstabilidad4C,vec(M))

#file=string("Est4C-",xini,"_",xfin,"-Vern9.jld")
#save(file,"c",c)
