"""
Implement the method collocation adaptativa for solve system of second order
differential equations  x''(t)=f(t,x(t)) (poner cita).
Usage
======
collocation2(f,p, t,z₀,v₀,t₀; dt=ddd,order=n,iter=3,tol=eps,dtmin=xxx,dtmax=yyy)

f function defined the system arguments t, x. t is a vector (column) and x a 2d array
(Matrix). The j-row correpond to the value of x in t[j].
p: parameter for function
t: vector of times where we want to find the solutions must to be t₀<=t[1],
z₀: initial z(t₀)
v₀: z'(t₀)
t₀: initial time

optional argeuments
===================
dt: size initial step
order: the orden method (6-10)
tol=eps: tolerance (order=8 tol=1e-19)
dtmin= minimum step
dtmax=maximun step
message= message in progress bar

"""

function collocation(f::Function,p, t::Array{Float64},z₀::Array{Float64},v₀::Array{Float64},t₀::Float64; dt=1e-3,order=10,iter=3,tol=1e-19,dtmin=0,dtmax=Inf,message="")
    @assert(t₀<=t[1],"t₀ must does no to be less than t[1] ")


    t_len=length(t);
    z_len=size(z₀,2)
    z=zeros(t_len,z_len)
    v=zeros(t_len,z_len);
    ind_out=1; #Number of t solved

    ProBar = Progress(t_len, 1,message)


    # Precomputing certain matrix
    file = matopen("matrices.mat")
        M=read(file, "Matriz_iteracion")
        D=read(file, "Matriz_diferencias")
    close(file)
    M=M[1:order-1,1:order-1];
    D=D[1:order-1,1:order-1];
    MD=M*D;







    n=convert(Float64, order)
    l=collect( 2:order)'
    exponentes=repmat(l,order-1,1)
    j=collect(0:order-2)
    matriz_epoca=(n-2)^(-2)./repmat(l.*(l-1),order-1,1)
    matriz_epocas_escala=matriz_epoca.*repmat(j,1,order-1).^exponentes;
    matriz_velocidad=((n-2).^j )./(j+1)
    matriz_epoca_vel=(n-2)^(-1)./(l-1);

    #Coefficient adaptation
    V=(l+order-3)./(l-1)
    coef_adap=factorial(order-2)*abs(V*M[:,order-1])

    t₁=t₀+dt;
    tₚ=collect(t₀:dt/(n-2):t₁);
    zₚ1=  repmat(z₀,order-1,1)+(tₚ-t₀)*v₀
    zₚ= zₚ1 +.5*(tₚ-t₀).^2*f(t₀,z₀,p);

    matriz_fuerza=f(tₚ,zₚ,p);


    estado_epocas=similar(zₚ1)

    global C
    for i ∈ 2:(order-1)

        C=MD*matriz_fuerza;
        estado_epocas=zₚ1+dt^2*matriz_epocas_escala*C;
        matriz_fuerza=f(tₚ,estado_epocas,p);
    end
    Ind=(t.<t₁).*(t.>=t₀)
    t_aux=t[Ind]
    k_aux=length(t_aux)
    if k_aux>0
        Z1=repmat(z₀,k_aux,1)+(t_aux-t₀)*v₀
        Z2=dt^2*taylor((order-2)*dt^-1*(t_aux-t₀),order,2)
        Z3=diagm(vec(matriz_epoca[1,:]))*C
        z[ind_out:ind_out+k_aux-1,:]=Z1+Z2*Z3

        V1=repmat(v₀,k_aux,1)
        V2=dt*taylor((order-2)*dt^-1*(t_aux-t₀),order,1);
        V3=diagm(vec(matriz_epoca_vel))*C
        v[ind_out:ind_out+k_aux-1,:]=V1+V2*V3
    end

    ind_out=ind_out+k_aux
    if ind_out>t_len
        return z,v
    end

    iteraciones=1

    ##################################################################
    ##COMIENZO DEL CICLO PRINCIPAL
    ##################################################################5
    while ind_out<=t_len
        dt_pas=dt
        dt=min(dtmax,max(dtmin,(n-2)*(tol*dt_pas^(n-2)/(coef_adap*maximum(abs.(C[order-1,:]))))^(1/(n-1))));
        r=dt/dt_pas
        t₀_pas=t₀
        t₀=t₁
        t₁=t₁+dt

        tₚ=collect(t₀:dt/(n-2):t₁)

        MV1=(n-2)+r*collect(0:(n-2))

        MV=repmat(MV1,1,order-1).^exponentes;

        zₚ=repmat(z₀,order-1,1)+(tₚ-t₀_pas)*v₀+dt_pas^2*matriz_epoca.*MV*C;

        z₀=zₚ[1,:]'
        v₀=v₀+dt_pas*matriz_velocidad'*C
        zₚ1=  repmat(z₀,order-1,1)+(tₚ-t₀)*v₀



        for k ∈ 1:iter
            matriz_fuerza=f(tₚ,zₚ,p);
            C=MD*matriz_fuerza;
            estado_epocas=zₚ1+dt^2*matriz_epocas_escala*C;
            matriz_fuerza=f(tₚ,estado_epocas,p);
        end


        Ind=(t.<t₁).*(t.>=t₀)
        t_aux=t[Ind]
        k_aux=length(t_aux)

        if k_aux>0

            Z1=repmat(z₀,k_aux,1)+(t_aux-t₀)*v₀
            Z2=dt^2*taylor((order-2)*dt^-1*(t_aux-t₀),order,2)
            Z3=diagm(vec(matriz_epoca[1,:]))*C
            z[ind_out:ind_out+k_aux-1,:]=Z1+Z2*Z3

            V1=repmat(v₀,k_aux,1)
            V2=dt*taylor((order-2)*dt^-1*(t_aux-t₀),order,1);
            V3=diagm(vec(matriz_epoca_vel))*C
	        v[ind_out:ind_out+k_aux-1,:]=V1+V2*V3
        end

        ind_out=ind_out+k_aux
        update!(ProBar , ind_out)
    end




    return z,v

end


function taylor(t,order,exp_ini)
    return repmat(t,1,order-1).^repmat((exp_ini:order-2+exp_ini)',length(t),1);
end
