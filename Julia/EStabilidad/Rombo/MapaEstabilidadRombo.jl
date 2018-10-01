function MapaEstabilidadRombo(Λ,z)
#x=linspace(0.01,0.41,50)
#v=linspace(50,60,50)
#Ejemplo de construcción de gráfico de a₁ y a₂
#fig = figure("Grafico",figsize=(10,10))
#ax = fig[:add_subplot](1,1,1, projection = "3d")
#ax[:plot_surface](xgrid,vgrid,min.(-real.(A₂),30), rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("jet"), alpha=0.8, linewidth=0.25)
#camino="/home/fernando/fer/Investigación/Trabajos GIT/Mecanica Celeste/Estabilidad Soluciones Periodicas Sitnikov Generalizado/Experimentos/4 cuerpos"
#file="/a1.png"
#savefig(string(camino,file))



include("EstabilidadRombo.jl")
l=length(z)
m=length(Λ)

cc=Array{Bool}(l,m)
A₁=Array{Complex{Float64}}(l,m)
A₂=Array{Complex{Float64}}(l,m)
T=Array{Float64}(l,m)
ProBar= Progress(m, 1)

for j in 1:m
    (a₁,a₂,TT,Δ)=EstabilidadRombo(Λ[j], z)
    cc[:,j]=isreal.(Δ).*(abs.(a₁).<=2).*(abs.(a₂).<=2)
    A₁[:,j]=a₁
    A₂[:,j]=a₂
    T[:,j]=TT
    update!(ProBar,j)
end


xgrid=repmat(Λ',l,1);
zgrid = repmat(z,1,m);
xestable=xgrid[cc];
zestable=zgrid[cc] ;


file=string("EstRombo","Vern9.jld")

save( file,"Λ",Λ,"z₀",z,"xestable",xestable,"zestable",zestable,"a1",A₁,"a2",A₂,"T",T)


#fig, ax=subplots()
#ax[:scatter](xestable,zestable,color="black", marker=:.)
#xlabel(L"\lambda")
#ylabel("z")
#title("Stability Map")
return xestable, zestable,A₁,A₂
end
