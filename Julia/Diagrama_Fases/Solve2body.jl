include("SolKepler.jl")
e=.1
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

using Interpolations
t=0:.01:T;
sol_x_a = interpolate(x, BSpline(Cubic(Line())), OnGrid())
sol_x=scale(sol_x_a,t)

