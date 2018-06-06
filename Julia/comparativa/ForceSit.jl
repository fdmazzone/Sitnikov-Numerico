
""" Requiere crear la función r con Func_r del archivo Solve2body.jl
    using DiferentialEquations
"""

m=[.5,.5]; #debe ser vector 1d
s=[-.5 .5]; # asi es un array 2d 1x2
function ForceSit(t,z,r)
    return -(z.^2.+r.(t).^2*s.^2 ).^(-1.5)*m.*z
end
