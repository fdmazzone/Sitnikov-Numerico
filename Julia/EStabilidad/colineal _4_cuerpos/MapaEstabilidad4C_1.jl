
#Requiere using PyPlot, y Estabilidad-n-colineales.jl
# Las posiciones de los cuatro primarios son [-1 -x_j  x_j 1]. Se debe ingresar un vector con los valores x_j
#v=vector de posiciones iniciales de la partículo no grave.

function MapaEstabilidad4C(x,v)
l=length(v)
k=length(x)

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


xgrid=repmat(x',l,1);
vgrid = repmat(v,1,k);


xestable=xgrid[cc];
vestable=vgrid[cc] ;
return xestable, vestable,A₁,A₂
end

#fig, ax=subplots()
#ax[:scatter](xestable,vestable,color="black", marker=:.)
#xlabel(L"x")
#ylabel(L"z(0)")
#title("Stability Map")
