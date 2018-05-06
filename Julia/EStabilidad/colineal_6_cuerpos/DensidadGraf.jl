#using PyPlot, JLD
camino="/home/fernando/fer/Astronomia/Simulaciones/Estabilidad Sitnikov/MapaEstabilidadDatos6CC-";

camino2="/home/fernando/fer/Investigación/Trabajos GIT/Mecanica Celeste/Estabilidad Soluciones Periodicas Sitnikov Generalizado/Experimentos/6 cuerpos/DensidadEstabilidad-6CC-";

#Data=load(camino)
#using StatsBase
#include hallaMAsas.jl

for h in 1:7
    file=string(camino,h,".jld")
    Data=load(file)
    xest=Data["x"]
    yest=Data["y"]
    zest=Data["z(0)"]


    XY=[[xest[i],yest[i]] for i in 1:size(xest)[1]];
    ContXY = countmap(XY);
    b = hcat([[key, val] for (key, val) in ContXY]...)'
    x=linspace(.1,.99,50)
    y=linspace(.1,.99,50)

    z=zeros(50,50)
    for i in 1:50
        for j in 1:50
            if [x[j] y[i]] in b[:,1]
                z[i,j]=ContXY[[x[j], y[i]]]
            end
        end
    end
              
    k=size(x)[1]
    n=size(y)[1]
    M = [ ColinealInv([-1 -x1 -x1*y1 x1*y1 x1 1]) for x1 in x, y1 in y ]
    EsPosit=[any(k->k<0,l) for l in M];
    Ind=collect(Iterators.product(1:k, 1:n))
    Ind=Ind[EsPosit]
    zz=zeros(50,50)
    for i in Ind
          z[i[2],i[1]]=z[i[2],i[1]]-100
    end
      
    fig1,ax1=subplots()
    cf1=ax1[:pcolormesh](x,y,z/1000, cmap=ColorMap("Greys"))
    fig1[:colorbar](cf1,ax=ax1)
    ylim(.1,.51)
    xlabel(L"x")
    ylabel(L"y")
    title(string("Stability Density for ",(h-1)*10+20,L"\leq z(0)\leq",h*10+20))
    savefig(string(camino2,h,".png"))
end
