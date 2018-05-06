#Resuelve la ecuación de kepler E=M-e*sen(E)
function solkepler(M::Array{Float64},e::Float64,tol::Float64)
@assert((0<=e && e<1),"Range excentricity is [0,1)")




filtro=1:length(M);
delta=tol+1
E=deepcopy(M);

while delta>tol 
    E1=M[filtro]+e.*sin.(E[filtro]);
    Error=abs.(E[filtro]-E1);
    delta=maximum(Error);
    E[filtro]=E1;
    filtro=filtro[Error.>tol];
end
return E

end
