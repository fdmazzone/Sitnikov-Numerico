using PyPlot

using LaTeXStrings

include("/home/fernando/fer/Investigación/Trabajos GIT/Mecanica Celeste/Estabilidad Soluciones Periodicas Sitnikov Generalizado/Experimentos/HallaMasas.jl")


x=linspace(.01,.99,100)

m=[ColinealInv([-1 -y y 1]) for y in x]

m₁=[k[1] for k in m];

m₂=[k[2] for k in m];

fig, ax = subplots()
ax[:plot](x, m₁, "r-", linewidth=2, label=L"m_1", alpha=0.6)
ax[:legend](loc="upper center")

ax[:plot](x, m₂, "b-", linewidth=2, label=L"m_2", alpha=0.6)
ax[:legend](loc="upper center")
