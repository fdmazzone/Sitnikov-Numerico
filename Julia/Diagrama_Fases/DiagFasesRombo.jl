λ=1.3
#e_lista=.1:.05:.6
e_lista=.0:.05:.5
size_e=length(e_lista)
z0_lista=collect(1:.02:5);
size_z0s=length(z0_lista)



m₂=4*(1+λ^2)^(3/2)*(8*λ^3-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
m₁=4*λ^3*(1+λ^2)^(3/2)*(8-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
m=[m₁,m₂]
s=[λ 1]



include("Solve2body.jl")


function ForceSitRombo(t,z,r)
	return -2(z.^2.+r.(t).^2*s.^2 ).^(-1.5)*m.*z
end

t₀=0.0
#CantPer=1000
CantPer=1000

z₀=Array{Float64}(1,1)
v₀=Array{Float64}(1,1)
v₀[1,1]=0.0
control=-1


z=Array{Float64}(CantPer+1,size_z0s)
v=Array{Float64}(CantPer+1,size_z0s)

Cont=1
for j_e in 1:size_e
	r,T=Func_r(e_lista[j_e]);
	t=collect(linspace(0,CantPer*T,CantPer+1));



	for j ∈ 1:size_z0s
		control=-control
		z₀[1,1]=control*z0_lista[j]
		total=Float16(Cont/size_z0s/size_e)
		Cont=Cont+1
		msj=string("e=",e_lista[j_e],", z₀ =",z₀[1,1], ", Total= ",total, "%, Integrador= ")
		zz,vv=collocation(ForceSitRombo,r,t,z₀,v₀,t₀;dt=.001,order=10,tol=1e-19,dtmin=.0001,dtmax=10,message=msj);
		#ax[:scatter](z,v,s=1)
		z[:,j]=zz
		v[:,j]=vv

	end
	file=string("Datos_e_",e_lista[j_e],".jld")
	save(file,"z",z,"v",v,"e",e_lista[j_e],"z₀",z0_lista,"T",T)
end





#fig, ax=subplots()
#xlabel(L"z")
#ylabel(L"V")
#title(L"\lambda="*string(λ)*","*L"e="*string(e))
#xlim(-6,6)
#ylim(-4,4)


#etiqueta=string("DiagFases_e_",e,"_lambda_",λ ,".svg")
#savefig(etiqueta)
