function f=Force2b(x,funcion_datos)
    m=1;
    D=diag(sum(x.^2,2).^-1.5) ;
    f=-m*D*x;
end