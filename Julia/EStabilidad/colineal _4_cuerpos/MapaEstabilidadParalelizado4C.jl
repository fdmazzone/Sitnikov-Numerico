#cd("/home/fernando/fer/Investigación/Trabajo en curso/Mecanica #Celeste/Julia/EStabilidad/colineal _4_cuerpos")

#procesos=7
#addprocs(procesos-1)

@everywhere using DifferentialEquations, QuadGK, JLD, ProgressMeter
#@everywhere include("MapaEstabilidad4C.jl")

@everywhere include("/home/fernando/fer/Investigación/Trabajo en curso/Mecanica Celeste/Julia/EStabilidad/Estabilidad-n-colineales.jl")
@everywhere include("/home/fernando/fer/Investigación/Trabajo en curso/Mecanica Celeste/Julia/EStabilidad/HallaMasas.jl")

@everywhere function MapaEstabilidad4C(X,p=nothing)
    x=X[1]
    v=X[2]
    l=length(v)
    k=length(x)
    cc=Array{Bool}(l,k)
    A₁=Array{Complex{Float64}}(l,k)
    A₂=Array{Complex{Float64}}(l,k)
    T=Array{Any}(l,k)
    @parallelprogress p string("Procesando: ",x[1],"-",x[end]," ") for j in 1:k
        m=ColinealInv([-1 -x[j]  x[j] 1])
        (a₁,a₂,t,Δ)=Estabilidad(m[3:4],[x[j] 1], v)
        cc[:,j]=isreal.(Δ).*(abs.(a₁).<=2).*(abs.(a₂).<=2)
        A₁[:,j]=a₁
        A₂[:,j]=a₂
        T[:,j]=t
    end
    xgrid=repmat(x',l,1);
    vgrid = repmat(v,1,k);
    xestable=xgrid[cc];
    zestable=vgrid[cc] ;
    return xestable, zestable,A₁,A₂, T
end







xini=0.38
#xfin=0.41
xfin=.39
vini=15
vfin=16
#Δx=.0001
Δx=.001
X=xini:Δx:xfin
Δv=.1
v=vini:Δv:vfin
N=length(X)
NxProc=Int64(floor(N/procesos))

j_ini=1
j_fin=NxProc





M=Array{Any}(procesos)
for i in 1:procesos
    x=X[j_ini:j_fin]
	M[i]=[x,v]
    j_ini=j_fin+1
    j_fin=j_fin+NxProc
end
c=ProgressMeter.pmap(MapaEstabilidad4C,M)

xestable=c[1][1]
zestable=c[1][2]
a₁=c[1][3]
a₂=c[1][4]
T=c[1][5]

for j in 2:procesos
    xestable=cat(1,xestable,c[j][1])
    zestable=cat(1,zestable,c[j][2])
    a₁=cat(2,a₁,c[j][3])
    a₂=cat(2,a₂,c[j][4])
    T=cat(1,T,c[j][5])
end


file=string("Est4C-Parale-Vern9.jld")
save( file,"x",X,"z₀",v,"xestable",xestable,"zestable",zestable,"a1",a₁,"a2",a₂,"T",T)



#Requiere using PyPlot, DifferentialEquations,  JLD, QuadGK
#ProgressMeter
#y Estabilidad-n-colineales.jl
# Las posiciones de los cuatro primarios son [-1 -x_j  x_j 1]. Se debe #ingresar un vector con los valores x_j
#v=vector de posiciones iniciales de la partículo no grave.






