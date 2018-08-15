using DifferentialEquations, JLD, Interpolations, ProgressMeter
include("DiagFasesRombo-DiffEq.jl")
include("Solve2body.jl")
λ=1.3
e=.1

#Experimento 1
v0_lista=[0.0]
z0_lista=linspace(5.73,5.8,200)
CantPer=3000;
nombre_arch="CantToriIsla1:10.jld"
FasesRombo(λ,e,v0_lista,z0_lista,CantPer,nombre_arch)

#Experimento 2

z0_lista=linspace(3.22,3.31,200)
nombre_arch="CantToriIslaPrin.jld"
FasesRombo(λ,e,v0_lista,z0_lista,CantPer,nombre_arch)

#Experimento 3
v0_lista=linspace(3.62345,3.6237,200)
z0_lista=[0,0]
CantPer=3000;
nombre_arch="CantToriIsla1Res1:1_zoomx4_cerca_ines.jld"
FasesRombo(λ,e,v0_lista,z0_lista,CantPer,nombre_arch)
