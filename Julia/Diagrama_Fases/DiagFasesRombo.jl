for λ ∈ [1.7 1.3 1]
	for e ∈ [.1 .3 .6]

		m₂=4*(1+λ^2)^(3/2)*(8*λ^3-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
		m₁=4*λ^3*(1+λ^2)^(3/2)*(8-(1+λ^2)^(3/2))/(64*λ^3-(1+λ^2)^3)
		m=[m₁,m₂]
		s=[λ 1]
		include("Solve2body.jl")

		r,T=Func_r(e);



		function ForceSitRombo(t,z,r)
		    return -2(z.^2.+r.(t).^2*s.^2 ).^(-1.5)*m.*z
		end

		t₀=0.0
		t=collect(0:T:500*T);
		z₀=Array{Float64}(1,1)
		v₀=Array{Float64}(1,1)
		z₀[1,1]=1
		v₀[1,1]=0.0


		fig, ax=subplots()
		xlabel(L"z")
		ylabel(L"V")
		title(L"\lambda="*string(λ)*","*L"e="*string(e))
   		xlim(-6,6)
		ylim(-4,4)
		control=-1
		for j ∈  .3:.01:4
			control=-control
		    z₀[1,1]=control*j
		    z,v=collocation2(ForceSitRombo,r,t,z₀,v₀,t₀;dt=.001,order=10,tol=1e-19,dtmin=.0001,dtmax=10);
		    ax[:scatter](z,v,s=1)
		end
		etiqueta=string("DiagFases_e_",e,"_lambda_",λ ,".svg")
		savefig(etiqueta)
	end
end
