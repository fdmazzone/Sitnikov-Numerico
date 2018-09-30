
#Requiere using PyPlot, DifferentialEquations,  JLD, QuadGK
#ProgressMeter
#y Estabilidad-n-colineales.jl
# Las posiciones de los cuatro primarios son [-1 -x_j  x_j 1]. Se debe #ingresar un vector con los valores x_j
#v=vector de posiciones iniciales de la partículo no grave.

function MapaEstabilidad4C(X)

x=X[1]
v=X[2]
l=length(v)
k=length(x)
include("/home/fernando/fer/Investigación/Trabajo en curso/Mecanica Celeste/Julia/EStabilidad/Estabilidad-n-colineales.jl")
include("/home/fernando/fer/Investigación/Trabajo en curso/Mecanica Celeste/Julia/EStabilidad/HallaMasas.jl")
cc=Array{Bool}(l,k)
A₁=Array{Complex{Float64}}(l,k)
A₂=Array{Complex{Float64}}(l,k)
T=Array{Any}(l,k)
ProBar= Progress(k, 1)

for j in 1:k
    update!(ProBar,j)
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



#file=string("Est4C-",x[1],"_",x[end],"-Vern9.jld")

#save( file,"x",x,"z₀",v,"xestable",xestable,"zestable",zestable,"a1",A₁,"a2",A₂,"T",T)


return xestable, zestable,A₁,A₂, T
end

#fig, ax=subplots()
#ax[:scatter](xestable,vestable,color="black", marker=:.)
#xlabel(L"x")
#ylabel(L"z(0)")
#title("Stability Map")
