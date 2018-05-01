#Problema tres cuerpos
#using OrdinaryDiffEq, Plots, QuadGK

# Three Body
n=3
@assert n>0 "Cantidad de cuerpos negativo"
θ=pi/n
Rₙ=1/(2*sin(θ))
m=(2*(sin(θ))^3*(sum(1/(sin(i*θ)) for i in 1:n-1)))^(-1)

ν=1
M=[1 0 0 0;0 -1 0 0; 0 0 -1 0;0 0 0 1]
u1=eye(4)
Fuerza3c = (t,u,du) -> begin
  du[1] = u[2]
  du[2] = -n*m*u[1]/(Rₙ^2+u[1]^2)^(1.5)
end
#Altura máxima de la partícula

#function CoefEst(z₀)
function CoefEst(v₀)
#Periodo de la partícula

  E=v₀^2/2-n*m/Rₙ
  #E=(z₀^2+1/4)^(-.5)
  #z₀=  v₀*sqrt(8-v₀^2 )/(4-v₀^2 )/2

  f(z)=1/sqrt(E+n*m*(Rₙ^2+z^2)^(-0.5))

  zₑ=sqrt((n*m/E)^2-Rₙ^2)

  T₀=quadgk(f,0,zₑ-1e-12)

  T= T₀[1]/2^.5
  # Sistema no lineal

  #u0 = [z₀,0.0]
  u0 = [0.0,v₀]
  tspan = (0.0,T)

  prob_ode_trescuerpos = ODEProblem(Fuerza3c,u0,tspan)


  sol = solve(prob_ode_trescuerpos)

  #Matríz Variacional
  Φ(t)=0.25*csc(θ)^2+sol(t)[1]^2
  F(t)=1-n*m/Φ(t)^(1.5)+3*n*m*csc(θ)^2/(8*Φ(t)^(2.5))

  #F₂(t)=1-1/Φ(t)^(3/2)
  A(t)=[0 0 1 0;0 0 0 1; F(t) 0 0 2*ν; 0 F(t) -2*ν 0]

  #Sistema Ecuaciones Variacionales
  variacional(t,u)=A(t)*u


  prob_ode_variacional = ODEProblem(variacional,u1,tspan)
  so_var = solve(prob_ode_variacional)
  #Matríz monodromía

  G=(M*inv(so_var(T))*M*so_var(T))^2
  #Coeficientes de estabilidad
  p=-trace(G)
  q=.5*(p^2-trace(G^2))
  Δ=sqrtm(p^2-4(q-2))
  a₁=.5*(p+Δ)
  a₂=.5*(p-Δ)
  return (a₁,a₂,T)
end
#
v=linspace(1.3,	1.8,	1000)
k=length(v)
a₁=Array{Complex{Float64}}(k)
a₂ =Array{Complex{Float64}}(k)
T =Array{Float64}(k)
for i=1:k
  a₁[i],a₂[i],T[i]=CoefEst(v[i])
end
