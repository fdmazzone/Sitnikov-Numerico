%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Esta  funci�n implementa un m�todo de colocaci�n adaptativo para      %%
%% resolver problemas aut�nomos de segundo orden. En particular se aplica%%
%% al problema de los n-cuerpos                                                                %%
%% efemerides_epocas= epocas en la que se quiere almacenar la soluci�n   %%
%%                                                               %%
%% iteraciones_pasos=cantidad de iteraciones para resolver la ecuac�n no %%
%% lineal en cada paso del m�todo                                        %%
%% orden=orden del m�todo                                                %%
%% delta_t=longitud del paso inicial                                     %%
%% tol=tolerancia para la adaptac�n                                      %%
%% posicion_out y velocidad_out solucion              %%
%% en efemerides_epocas                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [posicion_out,velocidad_out]=colocacion_nb_adap(funcion,funcion_datos,efemerides_epocas,posicion,velocidad,Jini,parametros)

%%%Extraccion parametros

iteraciones_pasos=parametros.iter;
orden=parametros.orden;
tol=parametros.tol;
delta_t=parametros.paso;
if not(efemerides_epocas(end)-Jini==0)
    delta_t=delta_t*sign(efemerides_epocas(end)-Jini);
end


mensaje_control=logical(0);
if length(fields(parametros))>4
    mensaje=parametros.mensaje;
    mensaje_control=logical(1);
end
    


%%Organizacion del tiempo

Jfin=efemerides_epocas(end);


%%Pre-alocando
cant=length(efemerides_epocas);
posicion_out=zeros(cant,size(posicion,2));
velocidad_out=zeros(cant,size(posicion,2));


indice_out=1; %contador que lleva registro de cuantas �pocas de las 
%efemerides fueron resueltas


     
    



%Matrices del metodo numerico
load('matrices');
Matriz_iteracion=Matriz_iteracion(1:orden-1,1:orden-1);
Matriz_diferencias=Matriz_diferencias(1:orden-1,1:orden-1);
Matriz=Matriz_iteracion*Matriz_diferencias;


l=[2:orden];
exponentes=repmat(l,[orden-1,1]);
j=(0:orden-2)';

matriz_epoca=(orden-2)^(-2)./repmat(l.*(l-1),[orden-1,1]);
matriz_epocas_escala=matriz_epoca.*(repmat(j,[1,orden-1]).^exponentes);
matriz_velocidad=[(orden-2).^(0:orden-2)]./(1:orden-1);
matriz_epoca_vel=(orden-2)^(-1)./((2:orden)-1);






%% Coeficientes para adaptacion

V=(orden-1:1:2*orden-3)./(1:orden-1);
coef_adap=factorial(orden-2)*abs(V*Matriz_iteracion(:,orden-1));

%Procedimiento de inicio, calculando la solucion entre las epocas de
%colocacion entre t(1) y t(2)






colocacion_epocas_ini=Jini;
colocacion_epocas_fin=Jini+delta_t;
colocacion_epocas=(colocacion_epocas_ini:delta_t/(orden-2):colocacion_epocas_fin)';

estado_epocas=ones(orden-1,1)*posicion+(colocacion_epocas-colocacion_epocas_ini)*velocidad...
    +.5*(colocacion_epocas-colocacion_epocas_ini).^2*funcion(posicion,funcion_datos);

matriz_fuerza=funcion(estado_epocas,funcion_datos);

taylor_2_term=ones(orden-1,1)*posicion+(colocacion_epocas-Jini)*velocidad;

for indi=2:1:(orden-1)
    C=Matriz*matriz_fuerza;
    estado_epocas=taylor_2_term+delta_t^2*matriz_epocas_escala*C;
    matriz_fuerza=funcion(estado_epocas,funcion_datos);
end
%%Si hay epocas de las efemerides dentro la almacenamos



Indices=find(sign(delta_t)*efemerides_epocas>=sign(delta_t)*colocacion_epocas_ini & sign(delta_t)*efemerides_epocas<sign(delta_t)*colocacion_epocas_fin);
if ~isempty(Indices)

    efemerides_epocas_aux=efemerides_epocas(Indices);
    
    posicion_out(indice_out:indice_out+length(Indices)-1,:)=ones(length(Indices),1)*posicion...
        +(efemerides_epocas(Indices)-colocacion_epocas_ini)*velocidad...
         +delta_t^2*monomios_taylor((orden-2)*delta_t^(-1)*(efemerides_epocas(Indices)-colocacion_epocas_ini),orden,2)*diag(matriz_epoca(1,:))*C;

     velocidad_out(indice_out:indice_out+length(Indices)-1,:)=ones(length(Indices),1)*velocidad...
         +delta_t*monomios_taylor((orden-2)*delta_t^(-1)*(efemerides_epocas(Indices)-colocacion_epocas_ini),orden,1)*diag(matriz_epoca_vel)*C;

    
    indice_out=length(Indices)+indice_out;
    if indice_out>cant
        return
    end
end
        





iteraciones=1; %contador para las iteraciones de la soluci�n del sistema 
%% no lineal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%COMIENZO DEL CICLO PRINCIPAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

if mensaje_control
    barra_progreso = waitbar(0,mensaje);
end

while indice_out<=cant
  
   %%nuevo delta
   delta_pas=delta_t;
    delta_t=sign(delta_t)*(orden-2)*(tol*delta_pas^(orden-2)/(coef_adap*max(abs(C(orden-1,:)))))^(1/(orden-1));
    r=delta_t/delta_pas; %razones entre epocas
  
%%nuevas epocas 
   colocacion_epocas_pas=colocacion_epocas_ini;
   colocacion_epocas_ini=colocacion_epocas_fin;
   colocacion_epocas_fin=colocacion_epocas_ini+delta_t;
   colocacion_epocas=(colocacion_epocas_ini:delta_t/(orden-2):colocacion_epocas_fin)';
   
   
   matriz_variable=repmat((orden-2:r:(orden-2)*(1+r))',[1,orden-1]).^exponentes;
   
%%Estimacion de nuevos estados a partir del paso anterior
    estado_epocas=ones(orden-1,1)*posicion...
        +(colocacion_epocas-colocacion_epocas_pas)*velocidad...    
        +delta_pas^2*matriz_epoca.*matriz_variable*C;
   
%%nueva posicion y velocidad
   posicion=estado_epocas(1,:);
   velocidad=velocidad+delta_pas*matriz_velocidad*C;
    
    
    
   while iteraciones<=iteraciones_pasos
        matriz_fuerza=funcion(estado_epocas,funcion_datos);
        C=Matriz*matriz_fuerza;
        estado_epocas=ones(orden-1,1)*posicion...
            +(colocacion_epocas-colocacion_epocas_ini)*velocidad...
            +delta_t^2*matriz_epocas_escala*C;
        iteraciones=iteraciones+1;
   end
   
   %Si hay epocas de la efemerides en el rango las calculamos
    Indices=find(sign(delta_t)*efemerides_epocas>=sign(delta_t)*colocacion_epocas_ini & sign(delta_t)*efemerides_epocas<sign(delta_t)*colocacion_epocas_fin);
    

    
    if ~isempty(Indices)

        %efmerides_epocas_aux=efemerides_epocas(Indices);
        posicion_out(indice_out:indice_out+length(Indices)-1,:)=ones(length(Indices),1)*posicion...
            +(efemerides_epocas(Indices)-colocacion_epocas_ini)*velocidad...
            +delta_t^2*monomios_taylor((orden-2)*delta_t^(-1)*(efemerides_epocas(Indices)-colocacion_epocas_ini),orden,2)*diag(matriz_epoca(1,:))*C;

        velocidad_out(indice_out:indice_out+length(Indices)-1,:)=ones(length(Indices),1)*velocidad...
            +delta_t*monomios_taylor((orden-2)*delta_t^(-1)*(efemerides_epocas(Indices)-colocacion_epocas_ini),orden,1)*diag(matriz_epoca_vel)*C;


        indice_out=length(Indices)+indice_out;
        if mensaje_control
            waitbar(indice_out/cant);
        end
    end
    iteraciones=1;


    
    
end
if mensaje_control
    close(barra_progreso);
end





end%fin funcion principal




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCIONES AUXILIARES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

function taylor=monomios_taylor(t,orden,exp_ini)

    taylor=repmat(t,[1,orden-1]).^repmat(exp_ini:orden-(2-exp_ini),[length(t),1]);
end
