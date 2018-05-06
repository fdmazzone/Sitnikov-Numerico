
#Requiere using PyPlot, y Estabilidad-n-colineales.jl
#El input de esta función es un vector X de posiciones de los primarios intermedios.
#[-1 -x_j  x_j 1]
#v=vector de posiciones iniciales de la partículo no grave.

function MapaEstabilidad4C(x,v)
#x=linspace(0.01,0.41,50)
#v=linspace(50,60,50)
l=size(v)[1]
k=size(x)[1]


cc=Array{Bool}(l,k)
A₁=Array{Complex{Float64}}(l,k)
A₂=Array{Complex{Float64}}(l,k)

for j in 1:k
    print(j/k*100,"% de avance")
    m=ColinealInv([-1 -x[j]  x[j] 1])


    (a₁,a₂,T,Δ)=Estabilidad(m[3:4],[x[j] 1], v)
    cc[:,j]=isreal.(Δ).*(abs.(a₁).<=2).*(abs.(a₂).<=2)
    A₁[:,j]=a₁
    A₂[:,j]=a₂
end

print(size(x))

print(size(v))


xgrid=repmat(x',l,1);
vgrid = repmat(v,1,k);


xestable=xgrid[cc];
vestable=vgrid[cc] ;


#fig, ax=subplots()
#ax[:scatter](xestable,vestable,color="black", marker=:.)
#xlabel(L"x")
#ylabel(L"z(0)")
#title("Stability Map")
return xestable, vestable,A₁,A₂
end
