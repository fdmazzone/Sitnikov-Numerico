λ=1.3



m₂=4*(1+λ^2)^(3/2)*(8*λ^3-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
m₁=4*λ^3*(1+λ^2)^(3/2)*(8-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
m=[m₁,m₂]
s=[λ 1]



include("Solve2body.jl")


function ForceSitRombo(t,z,r)
	return -2(z.^2.+r.(t).^2*s.^2 ).^(-1.5)*m.*z
end


r,T=Func_r(.2);

cont=1
z₀=Array{Float64}(1,1)
v₀=Array{Float64}(1,1)
v₀[1,1]=0
CantPer=4000
t=collect(linspace(0,CantPer*T,CantPer+1));
k=length(t)



t₀=0.0

z0_lista=collect(3.373:.0001:3.381);
size_z0s=length(z0_lista)

z=Array{Float64}(k,size_z0s)
v=Array{Float64}(k,size_z0s)

for j ∈ 1:size_z0s
	z₀[1,1]=z0_lista[j]
	zz,vv=collocation(ForceSitRombo,r,t,z₀,v₀,t₀;dt=.001,order=10,tol=1e-19,dtmin=.0001,dtmax=10);
	z[:,cont]=zz
	v[:,cont]=vv
	cont=cont+1

end
file=string("CanToriZoom2.jld")
save(file,"z",z,"v",v)
