#include("SolKepler.jl")
#using Interpolations

function Func_r(e::Float64)
	r=1;
	a=r/(1-e);
	n=a^(-1.5);
	T=2*pi/n;
	t=Array(0:.01:T);
	M=n*t
	E=solkepler(M,e,.0001);
	dist_peri=1;
	a=dist_peri/(1-e);
	x=a*(cos.(E)-e);
	y=a*sqrt(1-e^2)*sin.(E);
	t=0:.01:T;
	sol_x_a = interpolate(x, BSpline(Cubic(Line())), OnGrid())
	sol_x=scale(sol_x_a,t)
	sol_y_a = interpolate(y, BSpline(Cubic(Line())), OnGrid())
	sol_y=scale(sol_y_a,t)
	r(t)=sqrt(sol_x[t-floor(t/T)*T]^2+sol_y[t-floor(t/T)*T]^2)
	return r,T
end
