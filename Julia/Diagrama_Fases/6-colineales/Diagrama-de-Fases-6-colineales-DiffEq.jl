#using DifferentialEquations,  JLD, Interpolations, ProgressMeter, MAT
#posiciones=[(0.05040816326530612, 0.05040816326530612), (0.04, 0.04) , (0.1,0.4), (0.3, 0.3), (0.4,0.4), (0.6,0.4), (0.8,0.45) ]
posiciones=[(0.1,0.3)]
exentricidades=[0.1]
#exentricidades=[0.1]
cant_posiciones_iniciales=[.001]#0.1:0.02:5
#cant_posiciones_iniciales=4:0.2:5
n=100
A=Array{Float64}(size(cant_posiciones_iniciales)[1]*(n+1),2)
size_z0s=length(cant_posiciones_iniciales)
function Fuerza3c(du,u,p,t)
	m=p[1]
	s=p[2]
	r=p[3]
	R = (r(t)^2*s.^2 +u[1]^2).^(1.5)
	du[1] = u[2]
	du[2] = -2*sum(m./R*u[1])
end

for x ∈ posiciones


#for x ∈ []

	for e ∈ exentricidades
        x₁=x[1]
        y₁=x[2]

		#m=ColinealInv([-1 -x₁ -x₁*y₁ x₁*y₁ x₁ 1])
		#s=[-1 -x₁ -x₁*y₁	x₁*y₁ x₁ 1]

		m=ColinealInv([-1 -x₁ -x₁*y₁ x₁*y₁ x₁ 1])[4:6]
		s=[x₁*y₁ x₁ 1]
		#include("Solve2body.jl")
		alg=Vern7()
		r,T=Func_r(e);
		tspan=(0.0,n*T);

		z₀=Array{Float64}(1,1)
		v₀=Array{Float64}(1,1)
		z₀[1,1]=1
		v₀[1,1]=0.0

		control=-1
		z=Array{Float64}(n+1,size_z0s)
		v=Array{Float64}(n+1,size_z0s)

		z_len=length(cant_posiciones_iniciales)
		#k=1
		ProBar = Progress(z_len, 1)




		t=0:T:n*T;
		#tic()
		#ii=1
		for j ∈  1:size_z0s
			control=-control
			u0=[cant_posiciones_iniciales[j]*control,0.0]

			#Cont=Cont+1
			#msj=string("e=",e_lista[j_e],", z₀ =",z₀[1,1], ", Total= ",total, "%, Integrador= ")
			prob_ode_trescuerpos = ODEProblem(Fuerza3c,u0,tspan,[m,s,r]);

			sol = solve(prob_ode_trescuerpos,alg,dt=1e-3, progress=true,abstol = 1e-12, reltol = 1e-12);


			zz=sol.(t);z[:,j]=[a[1] for a in zz];v[:,j]=[a[2] for a in zz];
			#print(ii)
	        update!(ProBar,j)
		end

		#update!(ProBar , k)
		#k=k+1
		#toc()
        #A=hcat(z,v)


		#etiqueta=string("DiagFases_e_",e,"_lambda_",λ ,".svg")
		#savefig(etiqueta)
		#camino="C:/Users/Caro/Desktop/Diagrama de Fases 6 colineales/"
		file=string("DiffEq_Datos_e_",e,"alg",alg,"pos_ini",x,".jld")
		save(file,"z",z,"v",v,"e",e,"posini",x,"T",T)

	end

end
print("Listo!!!!")
