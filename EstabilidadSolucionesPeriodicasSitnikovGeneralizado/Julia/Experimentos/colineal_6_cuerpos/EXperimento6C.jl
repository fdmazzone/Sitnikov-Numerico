x=linspace(.1,0.99,5);
y=linspace(.1,0.99,5);
camino="/home/fernando/fer/Investigación/Trabajos GIT/Mecanica Celeste/Estabilidad Soluciones Periodicas Sitnikov Generalizado/Experimentos/6 cuerpos/"

for j in 1:8
    v=linspace(20+(j-1)*10,20+j*10,7)
    Ind, Xestable, Yestable, Vestable,A₁,A₂,T=MapaEstabilidad6C(x,y,v)
    file=string(camino, "MapaEstabilidadDatos6CC-",j,".jld")
    save(file,"Indice",Ind,"x",Xestable,"y",Yestable, "z(0)",Vestable,"a1",A₁,"a2",A₂,"T",T)
end
