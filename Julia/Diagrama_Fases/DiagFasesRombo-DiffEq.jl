
λ=1.3
e=.1
#e_lista=.0:.05:.5
#size_e=length(e_lista)
v0_lista=collect(linspace(-.1,.1,10));
size_v0s=length(v0_lista)
z0_lista=collect(linspace(5.80,6.1,10));
size_z0s=length(z0_lista)

m₂=4*(1+λ^2)^(3/2)*(8*λ^3-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
m₁=4*λ^3*(1+λ^2)^(3/2)*(8-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
m=[m₁,m₂]
s=[λ 1]



include("Solve2body.jl")

function Fuerza3c(du,u,p,t)
    m=p[1]
    s=p[2]
    r=p[3]
    R = (r(t)^2*s.^2 +u[1]^2).^(1.5)
    du[1] = u[2]
    du[2] = -2*sum(m./R*u[1])
end


t₀=0.0
CantPer=2000
z=Array{Float64}(CantPer+1,size_z0s,size_v0s)
v=Array{Float64}(CantPer+1,size_z0s,size_v0s)
#z=Array{Float64}(CantPer+1,size_z0s)
#v=Array{Float64}(CantPer+1,size_z0s)
control=-1

Cont=1
#for j_e in 1:size_e
	r,T=Func_r(e);
	tspan = (0.0,CantPer*T);
    alg =  Vern9()
    t=0:T:CantPer*T;
    #ProBar= Progress(size_z0s, 1)
    ProBar= Progress(size_z0s*size_v0s, 1)
    for i ∈ 1:size_v0s
    	for j ∈ 1:size_z0s
    		control=-control
    		#u0=[control*z0_lista[j],0.0]
            u0=[z0_lista[j],v0_lista[i]]
    		Cont=Cont+1
    		#msj=string("e=",e_lista[j_e],", z₀ =",z₀[1,1], ", Total= ",total, "%, Integrador= ")
    		prob_ode_trescuerpos = ODEProblem(Fuerza3c,u0,tspan,[m,s,r]);

    		sol = solve(prob_ode_trescuerpos,alg,dt=1e-3, progress=true,abstol = 1e-12, reltol = 1e-12);

            #zz=sol.(t);z[:,j]=[a[1] for a in zz];v[:,j]=[a[2] for a in zz];
            #update!(ProBar,j)
    		zz=sol.(t);z[:,j,i]=[a[1] for a in zz];v[:,j,i]=[a[2] for a in zz];
            update!(ProBar,j+(i-1)*10)
    	end
    end
	file=string("FasesRombo-e-.1_res-es-1:10-II.jld")
	save(file,"z",z,"v",v,"e",e,"z₀",z0_lista,"T",T)

#end





#fig, ax=subplots()
#xlabel(L"z")
#ylabel(L"V")
#title(L"\lambda="*string(λ)*","*L"e="*string(e))
#xlim(-6,6)
#ylim(-4,4)


#etiqueta=string("DiagFases_e_",e,"_lambda_",λ ,".svg")
#savefig(etiqueta)
