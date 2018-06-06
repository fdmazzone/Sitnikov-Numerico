function EstabilidadRombo(λ,z₀)
    include("/home/fernando/fer/Investigación/Trabajo en curso/Mecanica Celeste/Julia/metodo_colocacion/collocation.jl")
    global s,m,T,M,u1, v₀, t₀, t, tp

    t₀=0.0
    v₀=Array{Float64}(1,1)
    v₀[1,1]=0.0
    M=[1 0 0 0;0 -1 0 0; 0 0 -1 0;0 0 0 1]
    u1=eye(4)
    s=[λ 1]
    m₂=4*(1+λ^2)^(3/2)*(8*λ^3-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
    m₁=4*λ^3*(1+λ^2)^(3/2)*(8-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
    m=[m₁ m₂]
    r₀ = (s.^2 .+z₀.^2).^(.5)
    E=-2*sum(m./r₀)
    f(z)=1/sqrt(E+2sum(m./sqrt.(s.^2.+z.^2)))

    T₀=quadgk(f,0,z₀[1,1]-1e-20)
    T= T₀[1]/2^.5
    tp=linspace(0,T,1000)
    t=collect(tp);

    #r=[1.41421, 1.41421]
    #m=[1.04482, 1.04482]
    #E-2.9551845001450343
    #T=4.619566061527369

    a₁,a₂,T,Δ=CoefEstRombo(z₀)

    # k=length(v)
    # a₁=Array{Complex{Float64}}(k)
    # a₂ =Array{Complex{Float64}}(k)
    # T = Array{Float64}(k)
    # Δ =Array{Complex{Float64}}(k)
    # for i=1:k
    #   a₁[i],a₂[i],T[i],Δ[i]=CoefEstRombo(m,s,v[i])
    # end
    # return (a₁,a₂,T,Δ)
    return γ
end


function Fuerza3c(t,z,p)
    return -2(z.^2.+s.^2 ).^(-1.5)*m'.*z
end
#Altura máxima de la partícula

function VariationalSystem(z₀)

    z,v=collocation(Fuerza3c,[],t,z₀,v₀,t₀;dt=.001,order=10,tol=1e-22,dtmin=.0001,dtmax=10)
    z_inter_ad = interpolate(z[:], BSpline(Cubic(Line())), OnGrid())
    γ=scale(z_inter_ad,tp)
    Φ(t)=s.^2.+γ[t]^2
    F(t)=1-2*sum(m.*(Φ(t).^(-1.5)))
    F₂(t)=F(t)+6*m[2]*s[2]^2*Φ(t)[2]^(-2.5)
    F₁(t)=F(t)+6*m[1]*s[1]^2*Φ(t)[1]^(-2.5)

    function VarFor(t,u,p)
        l=Int64(length(u)/4)
        v=reshape(u,(l,4))
        y=v*[0 0 F₁(t) 0; 0 0 0 F₂(t); 1 0 0 -2; 0 1 2 0 ]
        return y[:]'
    end

   return VarFor
end



#
#
#
#
 function CoefEstRombo(z₀)


      VarFor=VariationalSystem(z₀)                                                                                                                                                                                                                                                         

#
#
#     #
#     #   prob_ode_variacional = ODEProblem(variacional,u1,tspan)
#     #   so_var = solve(prob_ode_variacional)
#     #   #Matríz monodromía
#     #
#     #   G=(M*inv(so_var(T))*M*so_var(T))^2
#     #   #Coeficientes de estabilidad
#     #   p=-trace(G)
#     #   q=.5*(p^2-trace(G^2))
#     #
#     #   Δ=sqrtm(p^2-4(q-2))
#     # # pensar definición a_1 y a_2 para que sean funciones suaves
#     #   a₁=.5*(p+Δ)
#     #   a₂=.5*(p-Δ)
#     #   return (a₁,a₂,T,Δ)
    return A
end
# #
