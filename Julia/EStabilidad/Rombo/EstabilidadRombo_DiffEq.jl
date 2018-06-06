function EstabilidadRombo(λ,v)


s=[λ,1]
m₂=4*(1+λ^2)^(3/2)*(8*λ^3-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
m₁=4*λ^3*(1+λ^2)^(3/2)*(8-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
m=[m₁,m₂]
k=length(v)
a₁=Array{Complex{Float64}}(k)
a₂ =Array{Complex{Float64}}(k)
T = Array{Float64}(k)
Δ =Array{Complex{Float64}}(k)
for i=1:k
  a₁[i],a₂[i],T[i],Δ[i]=CoefEstRombo(m,s,v[i])
end
    return (a₁,a₂,T,Δ)
end



function Fuerza3c(du,u,p,t)

    m=p[1]
    s=p[2]

  r = (s.^2 +u[1]^2).^(1.5)
  du[1] = u[2]
  du[2] = sum(m./r)*(-2*u[1])

  end
#Altura máxima de la partícula

function CoefEstRombo(m,s,z₀)
    M=[1 0 0 0;0 -1 0 0; 0 0 -1 0;0 0 0 1]
u1=eye(4)
#function CoefEst(v₀)
#Periodo de la partícula
  #E=v₀^2/2-2
r = (s.^2 +z₀^2).^(.5)
          E=-2*sum(m./r)

  #z₀=  v₀*sqrt(8-v₀^2 )/(4-v₀^2 )/2
f(z)=1/sqrt(E+2sum(m./sqrt.(s.^2+z^2)))
  T₀=quadgk(f,0,z₀-1e-10)

  T= T₀[1]/2^.5
  # Sistema no lineal


  u0 = [z₀,0.0]
  #u0 = [0,v₀]
  tspan = (0.0,T)

  prob_ode_trescuerpos = ODEProblem(Fuerza3c,u0,tspan,[m,s])

  sol = solve(prob_ode_trescuerpos)
  #Matríz Variacional
  # a=s[1]
  # b=s[2]
  # m₁=m[1]
  # m₂=m[2]
  #
  # Φ₁(t)=a^2+(sol(t)[1])^2
  # Φ₂(t)=b^2+(sol(t)[1])^2
  #
  #
  # F₁(t)=1-2m₁*(-3a^2/Φ₁(t)^(5/2)+1/Φ₁(t)^(3/2))-2m₂/Φ₂(t)^(3/2)
  # F₂(t)= 1-2m₂*(-3b^2/Φ₂(t)^(5/2)+1/Φ₂(t)^(3/2))-2m₁/Φ₁(t)^(3/2)


  Φ(t)=s.^2+(sol(t)[1])^2

# REVISAR



F(t)=1-2*sum(m.*(Φ(t).^(-1.5)))

F₂(t)=F(t)+6*m[2]*s[2]^2*Φ(t)[2]^(-2.5)
F₁(t)=F(t)+6*m[1]*s[1]^2*Φ(t)[1]^(-2.5)





  A(t)=[0 0 1 0;0 0 0 1; F₁(t) 0 0 2; 0 F₂(t) -2 0 ]

  #Sistema Ecuaciones Variacionales
variacional(u,p,t)=A(t)*u


  prob_ode_variacional = ODEProblem(variacional,u1,tspan)
  so_var = solve(prob_ode_variacional)
  #Matríz monodromía

  G=(M*inv(so_var(T))*M*so_var(T))^2
  #Coeficientes de estabilidad
  p=-trace(G)
  q=.5*(p^2-trace(G^2))

  Δ=sqrtm(p^2-4(q-2))
# pensar definición a_1 y a_2 para que sean funciones suaves
  a₁=.5*(p+Δ)
  a₂=.5*(p-Δ)
  return (a₁,a₂,T,Δ)
end
#
