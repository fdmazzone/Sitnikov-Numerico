#using PyPlot, JLD
#camino="/media/fernando/Seagate Expansion Drive2/fer/Astronomia/Simulaciones/Estabilidad Sitnikov/MapaEstabilidadDatos6CC-7.jld";
#Data=load(camino)


camino="/home/fernando/fer/Investigación/Trabajos GIT/Mecanica Celeste/Estabilidad Soluciones Periodicas Sitnikov Generalizado/Experimentos/6 cuerpos/"

x=Data["x"]
y=Data["y"]
z=Data["z(0)"]

I=1:10:size(z)[1];
for i in I
    ind_zi=find(x -> x == z[i],z[I])
    fig,ax=subplots()
    ind_zi=I[ind_zi]
    ax[:scatter](x[ind_zi],y[ind_zi])
    xlim(0,1)
    ylim(0,1)
    savefig(string(camino,"graf-",i,".jpg"))
    close(fig)
end

