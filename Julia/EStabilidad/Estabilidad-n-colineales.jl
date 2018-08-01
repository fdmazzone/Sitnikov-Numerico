#El input de esta función  son un vector de masas y posiciones de una CC. Como la configuracion debe ser admisible se supone que hay una cantidad par y solo hay que sumuinistrar un vector de masas correspondiente a posiciones positivas. El centro de masa se supone que es el origen.
#m= vector de masas correspondiente a las posiciones positivas.
#s= vector de posiciones positivas
#v=vector de posiciones iniciales de la partícula no grave.La velocidad inicial se supone cero.


# Requiere using DifferentialEquations, QuadGK, y HallaMasas.jl



function ColinealInv(x)
  #x vector fila de posiciones
  l=size(x)[2]
  D=ones(l,1)*x
  C=D-D'
  A=((abs.(C)).^(-3)).*C
  for i=1:l
    A[i,i]=0
  end
  return A\(-x')
end

function Estabilidad(m,s,v) #CARGAR SOLO LAS POSICIONES POSITIVAS

# m=ColinealInv(s)
# proveemos m=vector de masas, s=vector de posiciones
#CARGAR SOLO LAS POSICIONES POSITIVAS
#v=linspace(10.3,	15.13,	1000), vector de posiciones iniciaeles

k=length(v)
a₁=Array{Complex{Float64}}(k)
a₂ =Array{Complex{Float64}}(k)
T =Array{Float64}(k)
Δ =Array{Complex{Float64}}(k)
for i=1:k
  a₁[i],a₂[i],T[i],Δ[i]=CoefEst(m,s,v[i])
end

    return (a₁, a₂, Δ, T)
end


#Asi se programa las funciones para el solver F(du, variable independet , param, var. indepent.)
function Fuerza3c(du,u,p,t)

     m=p[1]
    s=p[2]

  r = (s.^2 +u[1]^2).^(1.5)
  du[1] = u[2]
  du[2] = sum(m./r)*(-2*u[1])
  end


  function CoefEst(m,s,z₀)
    M=[1 0 0 0;0 -1 0 0; 0 0 -1 0;0 0 0 1]
    u1=eye(4)

 r = (s.^2 +z₀^2).^(.5)

  E=-2*sum(m./r)

f(z)=1/sqrt(E+2sum(m./sqrt.(s.^2+z^2)))


  T₀=quadgk(f,0,z₀-1e-10)
  T= T₀[1]/2^.5

  # Sistema no lineal


  u0 = [z₀,0.0]
  #u0 = [0,v₀]
  tspan = (0.0,T)

  prob_ode_trescuerpos = ODEProblem(Fuerza3c,u0,tspan,[m,s])
    alg= Vern9();
  sol = solve(prob_ode_trescuerpos,alg)
  #Matríz Variacional

  Φ(t)=s.^2+(sol(t)[1])^2


F₂(t)=1-2*sum(m.*(Φ(t).^(-1.5)))
F₁(t)=F₂(t)+6*sum(m.*(s.^2.*Φ(t).^(-2.5)))



  A(t)=[0 0 1 0;0 0 0 1; F₁(t) 0 0 2; 0 F₂(t) -2 0 ]

#-------
  #Sistema Ecuaciones Variacionales
  variacional(u,p,t)=A(t)*u


  prob_ode_variacional = ODEProblem(variacional,u1,tspan)
  so_var = solve(prob_ode_variacional,alg)
  #Matríz monodromía


  G=(M*inv(so_var(T))*M*so_var(T))^2
  #Coeficientes de estabilidad
  p=-trace(G)
  q=.5*(p^2-trace(G^2))

  Δ=sqrtm(p^2-4(q-2))

  a₁=.5*(p+Δ)
  a₂=.5*(p-Δ)
  return (a₁,a₂,T,Δ)
end
