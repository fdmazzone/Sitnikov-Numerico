m=[.5,.5];
s=[-.5,.5];
function Fuerza3c(du,u,p,t)
    m=p[1]
    s=p[2]
    r=p[3]
    R = (r(t)^2*s.^2 +u[1]^2).^(1.5)
    du[1] = u[2]
    du[2] = -sum(m./R*u[1])
end
z₀=1
u0 = [z₀,0.0]
tspan = (0.0,100*T)
prob_ode_trescuerpos = ODEProblem(Fuerza3c,u0,tspan,[m,s,r]);
alg =  ARKODE(order=8)
sol = solve(prob_ode_trescuerpos,alg);
print("Fin")

