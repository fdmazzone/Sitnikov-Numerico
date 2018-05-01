function [posicion_out,velocidad_out]=multipaso_nb_implicito_noaut(funcion,funcion_datos,efemerides_epocas,posicion,velocidad,Jini,parametros)

orden=parametros.orden;
delta_t=parametros.paso;
delta_t=delta_t*sign(efemerides_epocas(end)-Jini);
mensaje_control=logical(0);
if length(fields(parametros))>2
    mensaje=parametros.mensaje;
    mensaje_control=logical(1);
end
    


cant=length(efemerides_epocas);


%%Pre-alocando
posicion_out=zeros(cant,size(posicion,2));
velocidad_out=zeros(cant,size(posicion,2));
indice_out=1;


     
 



%Matrices del metodo numerico
load('matrices');%colocacion, para procedimiento inicio
Matriz_diferencias=Matriz_diferencias(1:orden-1,1:orden-1);
Matriz_iteracion=Matriz_iteracion(1:orden-1,1:orden-1);%colocacion
l=[2:orden];
exponentes=repmat(l,[orden-1,1]);
j=(0:orden-2)';
matriz_epocas_escala=(orden-2)^(-2)*(repmat(j,[1,orden-1]).^exponentes)./repmat(l.*(l-1),[orden-1,1]);


load('matrices_multistep_explicito'); %multipaso implicito
Matriz_diferencias_backward=Matriz_diferencias_backward(1:orden-1,1:orden-1);
Matriz_iteracion_explicito=Matriz_iteracion_multistep(1:orden-1,1:orden-1);%multipaso
matriz_epoca=delta_t^2./((2:orden).*((2:orden)-1));
matriz_epoca_vel=delta_t./((2:orden)-1);
load('matrices_multistep_implicito','Matriz_iteracion_multistep'); %multipaso implicito
Matriz_iteracion_implicito=Matriz_iteracion_multistep(1:orden-1,1:orden-1);




%Procedimiento de inicio, calculando la solucion entre las epocas de
%colocacion entre t(1) y t(2) por un metodo de colocacion

delta_inicio=(orden-2)*delta_t;
colocacion_epocas_ini=Jini;
colocacion_epocas_fin=Jini+delta_inicio;
colocacion_epocas=(colocacion_epocas_ini:delta_inicio/(orden-2):colocacion_epocas_fin)';

taylor_2_term=ones(orden-1,1)*posicion+(colocacion_epocas-colocacion_epocas_ini)*velocidad;

estado_epocas=taylor_2_term...
    +.5*(colocacion_epocas-colocacion_epocas_ini).^2*funcion(Jini,posicion,funcion_datos);

matriz_fuerza=funcion(estado_epocas,funcion_datos);



for indi=2:1:(orden-1)
    D=Matriz_diferencias*matriz_fuerza;
    C=Matriz_iteracion*D;
    estado_epocas=taylor_2_term+delta_inicio^2*matriz_epocas_escala*C;
    matriz_fuerza=funcion(colocacion_epocas,estado_epocas,funcion_datos);
end
%%Si hay epocas de las efemerides dentro la almacenamos



Indices=find(sign(delta_t)*efemerides_epocas>=sign(delta_t)*colocacion_epocas_ini & sign(delta_t)*efemerides_epocas<sign(delta_t)*colocacion_epocas_fin);
if ~isempty(Indices)
    c2y=C_2_Y(orden,delta_inicio);
    Y=c2y*C;
    efemerides_epocas_aux=efemerides_epocas(Indices);
    
    posicion_out(indice_out:indice_out+length(Indices)-1,:)=ones(length(Indices),1)*posicion...
        +(efemerides_epocas(Indices)-colocacion_epocas_ini)*velocidad...
        +taylor_desarrollo(efemerides_epocas(Indices)-colocacion_epocas_ini,orden,2)*Y;

    velocidad_out(indice_out:indice_out+length(Indices)-1,:)=ones(length(Indices),1)*velocidad...
        +taylor_desarrollo(efemerides_epocas(Indices)-colocacion_epocas_ini,orden,1)*Y;
    
    indice_out=length(Indices)+indice_out;
    if indice_out>cant
        return
    end
end
 
%% Posicion y velocidad en la nueva epoca 
epoca=colocacion_epocas_fin;
posicion=estado_epocas(end,:);
velocidad=velocidad+delta_inicio*([(orden-2).^(0:orden-2)]./(1:orden-1))*C;
    

%%Reordenamos la matriz fuerza de manera descendente y calculamos las
%%diferencias divididas
matriz_fuerza=matriz_fuerza(orden-1:-1:1,:);
div_dif_fuerzas=Matriz_diferencias_backward*matriz_fuerza;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%COMIENZO DEL CICLO PRINCIPAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

if mensaje_control
    barra_progreso = waitbar(0,mensaje);
end

while indice_out<=cant
   D=Matriz_iteracion_explicito*div_dif_fuerzas;
   posicion_extrapolada=posicion+delta_t*velocidad+matriz_epoca*D;
   
   fuerza_extrapolada=funcion(posicion_extrapolada,funcion_datos);
   
   div_dif_fuerzas_extrapolacion(1,:)=fuerza_extrapolada-div_dif_fuerzas(1,:);
   for j=2:orden-2
        div_dif_fuerzas_extrapolacion(j,:)=div_dif_fuerzas_extrapolacion(j-1,:)-div_dif_fuerzas(j,:);
   end
   div_dif_fuerzas_extrapolacion=[fuerza_extrapolada;div_dif_fuerzas_extrapolacion(1:orden-2,:)];
   D=Matriz_iteracion_implicito*div_dif_fuerzas_extrapolacion;
   
   %Calculamos efemerides si caen en el intervalo.
   Indices=find(sign(delta_t)*efemerides_epocas>=sign(delta_t)*epoca & sign(delta_t)*efemerides_epocas<sign(delta_t)*(epoca+delta_t));
      if ~isempty(Indices)
        efmerides_epocas_aux=efemerides_epocas(Indices);
        posicion_out(indice_out:indice_out+length(Indices)-1,:)=ones(length(Indices),1)*posicion...
            +(efemerides_epocas(Indices)-epoca)*velocidad...
            +monomios_taylor(delta_t^(-1)*(efemerides_epocas(Indices)-epoca),orden,2)*diag(matriz_epoca)*D;

        velocidad_out(indice_out:indice_out+length(Indices)-1,:)=ones(length(Indices),1)*velocidad...
            +monomios_taylor(delta_t^(-1)*(efemerides_epocas(Indices)-epoca),orden,1)*diag(matriz_epoca_vel)*D;


        indice_out=length(Indices)+indice_out;
        waitbar(indice_out/cant);
      end
   
   epoca=epoca+delta_t;
   posicion=posicion+delta_t*velocidad+matriz_epoca*D;
   
   velocidad=velocidad+matriz_epoca_vel*D;

   
   fuerza_nueva_epoca=funcion(posicion,funcion_datos);
   div_dif_fuerzas(1,:)=fuerza_nueva_epoca-div_dif_fuerzas(1,:);
   for j=2:orden-2
        div_dif_fuerzas(j,:)=div_dif_fuerzas(j-1,:)-div_dif_fuerzas(j,:);
   end
   div_dif_fuerzas=[fuerza_nueva_epoca;div_dif_fuerzas(1:orden-2,:)];
   
end
if mensaje_control
    close(barra_progreso);
end


end%fin funcion principal




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCIONES AUXILIARES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

function taylor=taylor_desarrollo(t,orden,exp_ini)
    taylor=(repmat(t,[1,orden-1]).^repmat(exp_ini:orden-(2-exp_ini),[length(t),1]))*diag(factorial(exp_ini:1:orden-(2-exp_ini)).^(-1));
end
function taylor=monomios_taylor(t,orden,exp_ini)
    taylor=repmat(t,[1,orden-1]).^repmat(exp_ini:orden-(2-exp_ini),[length(t),1]);
end

function c2y= C_2_Y(orden,delta)
    M=diag(   [factorial(0:1:orden-2)].^(-1).*(delta/(orden-2)).^(0:1:orden-2) );
    warning off
    c2y=inv(M);
    warning on
end


