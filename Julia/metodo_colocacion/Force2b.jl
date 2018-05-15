function Force2b(t,x,p)
    m=1
    D=diagm(vec(sum(x.^2,2).^-1.5) )
    return -m*D*x
end
