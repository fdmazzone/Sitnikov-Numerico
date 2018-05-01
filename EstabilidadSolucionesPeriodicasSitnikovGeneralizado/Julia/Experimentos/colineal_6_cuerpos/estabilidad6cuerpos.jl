x=0.4
y=0.4
v=linspace(5,10,500)
    include("/home/fernando/fer/Investigación/Trabajos GIT/Mecanica Celeste/Estabilidad Soluciones Periodicas Sitnikov Generalizado/Experimentos/HallaMasas.jl")
    include("/home/fernando/fer/Investigación/Trabajos GIT/Mecanica Celeste/Estabilidad Soluciones Periodicas Sitnikov Generalizado/Experimentos/Estabilidad-n-colineales.jl")
    s=[-1 -x -x*y x*y x 1]
    m=ColinealInv(s)
    a=Estabilidad(m,s,v)

    plot(v,real(-a[1]))
    plot!(v,real(-a[2]))
