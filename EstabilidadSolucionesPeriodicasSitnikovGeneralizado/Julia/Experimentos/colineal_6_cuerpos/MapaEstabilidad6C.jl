#Requiere using PyPlot, y Estabilidad-n-colineales.jl ProgressMeter.jl
#El input de esta función es un vector X de posiciones de los primarios intermedios.
#[-1 -x_j  x_j 1]
#v=vector de posiciones iniciales de la partículo no grave.

function MapaEstabilidad6C(x,y,v)
#x=linspace(0.01,0.41,50)
#v=linspace(50,60,50)
l=size(v)[1]
k=size(x)[1]
n=size(y)[1]

#MasaPosi=Array{Bool}(k,n,l)
#A₁=Array{Complex{Float64}}(k,n,l)
#A₂=Array{Complex{Float64}}(k,n,l)

M = [ ColinealInv([-1 -x1 -x1*y1 x1*y1 x1 1]) for x1 in x, y1 in y ]
EsPosit=[~any(k->k<0,l) for l in M];

Ind=collect(Iterators.product(1:k, 1:n))
Ind=Ind[EsPosit]
#Estable=[Estabilidad(M[j[1],j[2]][4:6],[x[j[1]]*y[j[2]] x[j[1]] 1], v) for j in Ind ]
EsEstable=[]
A₁=[]
A₂=[]
Xestable=[]
Yestable=[]
Vestable=[]
TT=[]

###
 @showprogress 1 "Calculando ..." for j in Ind
        m=M[j[1],j[2]]

        (a₁,a₂,Δ,T)=Estabilidad(m[4:6],[x[j[1]]*y[j[2]] x[j[1]] 1], v)

        valorFV=isreal.(Δ).*(abs.(a₁).<=2).*(abs.(a₂).<=2)
        vest=v[valorFV]
        xest=x[j[1]]*ones(size(vest))
        yest=y[j[2]]*ones(size(vest))

        append!(Yestable,yest)
        append!(Xestable,xest)
        append!(Vestable,vest)
        append!(A₁,a₁)
        append!(A₂,a₂)
        append!(TT,T)
    end



#
# xgrid=repmat(x,l,1);
# ygrid=repmat()
# vgrid = repmat(v,1,k);
#
#
# xestable=xgrid[cc];
# vestable=vgrid[cc] ;
###


#fig, ax=subplots()
#ax[:scatter](xestable,vestable,color="black", marker=:.)
#xlabel(L"x")
#ylabel(L"z(0)")
#title("Stability Map")
#return xestable, vestable,A₁,A₂
return Xestable, Yestable, Vestable,A₁,A₂,TT
end
