x=linspace(.1,0.417220824736089,100);
v=linspace(10,100,7000)

xestable, vestable,A₁,A₂=MapaEstabilidad4C(x,v)

camino="/home/fernando/fer/Investigación/Trabajos GIT/Mecanica Celeste/Estabilidad Soluciones Periodicas Sitnikov Generalizado/Experimentos/4 cuerpos/MapaEstabilidadDatos.jld"
save(camino,"x",xestable, "z(0)",vestable,"a1",A₁,"a2",A₂)
