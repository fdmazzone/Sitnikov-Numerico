# Solve kepler equation  E=M-e*sen(E) by fixed point iteration
# solkepler(M::Array{Float64},e::Float64,tol::Float64)
# M: mean anomaly Array{Float64}
# e: excentricity Float64
# tol: tolerance, a difference less that tol finished iteration


#include("SolKepler.jl")
#using Interpolations

function Func_r(e::Float64)
	@assert((0<=e && e<1),"Range excentricity is [0,1)")
	dist_peri=1;
	a=dist_peri/(1-e);
	n=a^(-1.5);
	T=2*pi/n;
	t=Array(0:.01:T);
	M=n*t
	E=solkepler(M,e,.0001);
	x=a*(cos.(E)-e);
	y=a*sqrt(1-e^2)*sin.(E);
	t=0:.01:T;
	sol_x_a = interpolate(x, BSpline(Cubic(Line())), OnGrid())
	sol_x=scale(sol_x_a,t)
	sol_y_a = interpolate(y, BSpline(Cubic(Line())), OnGrid())
	sol_y=scale(sol_y_a,t)
	r(t)=sqrt(sol_x[t-floor(t/T)*T]^2+sol_y[t-floor(t/T)*T]^2)
	return r,T,sol_x,sol_y
end

function solkepler(M::Array{Float64},e::Float64,tol::Float64)

"""
Solve kepler equation  E=M-e*sen(E) by fixed point iteration
solkepler(M::Array{Float64},e::Float64,tol::Float64)
M: mean anomaly Array{Float64}
e: excentricity Float64
tol: tolerance, a difference less that tol finished iteration

"""

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

