function MapaEstabilidadRombo(Λ,v)
#x=linspace(0.01,0.41,50)
#v=linspace(50,60,50)
#Ejemplo de construcción de gráfico de a₁ y a₂
#fig = figure("Grafico",figsize=(10,10))
#ax = fig[:add_subplot](1,1,1, projection = "3d")
#ax[:plot_surface](xgrid,vgrid,min.(-real.(A₂),30), rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("jet"), alpha=0.8, linewidth=0.25)
#camino="/home/fernando/fer/Investigación/Trabajos GIT/Mecanica Celeste/Estabilidad Soluciones Periodicas Sitnikov Generalizado/Experimentos/4 cuerpos"
#file="/a1.png"
#savefig(string(camino,file))




l=size(v)[1]
m=size(Λ)[1]

cc=Array{Bool}(l,m)
A₁=Array{Complex{Float64}}(l,m)
A₂=Array{Complex{Float64}}(l,m)

for j in 1:m
    (a₁,a₂,T,Δ)=EstabilidadRombo(Λ[j], v)
    cc[:,j]=isreal.(Δ).*(abs.(a₁).<=2).*(abs.(a₂).<=2)
    A₁[:,j]=a₁
    A₂[:,j]=a₂
end


xgrid=repmat(Λ',l,1);
vgrid = repmat(v,1,m);
xestable=xgrid[cc];
vestable=vgrid[cc] ;


fig, ax=subplots()
ax[:scatter](xestable,vestable,color="black", marker=:.)
xlabel(L"\lambda")
ylabel("v")
title("Stability Map")
return xestable, vestable,A₁,A₂
end
