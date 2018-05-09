
function [pos,vel]=simulacion_sistema_ficticio(GM,efemerides_epocas,...
    posicion_ini,velocidad_ini)
%%%%%%%%%%%Ejemplo de Uso%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% m=[2,2,2];
% m=ones(1,15);
% X0=poligono(15);
% % X0=[1 0; -sqrt(3)/2 -1/2;-sqrt(3)/2 +1/2]';
% X=HallaCC(X0,m);
% X(3,:)=zeros([1,15]);
% J=[ 0 -1 0;1 0 0;0 0 0];
% V=J*X;
% m=[m,0];
% X(:,16)=[0 0 0];
% V(:,16)=[0,0,1];
% efemerides_epocas=(0:.1:5*2*pi)';
% [pos,vel]=simulacion_sistema_ficticio(m,efemerides_epocas,X(:)',V(:)');
% animacion(pos)

funcion=@fuerza_nb;


%%%Colocacion
Integrador=@colocacion_nb_adap; 
parametros_Integrador.paso=.001;
parametros_Integrador.tol=1e-27;
parametros_Integrador.orden=12;
parametros_Integrador.iter=3;

%multipaso
% Integrador=@multipaso_nb_implicito;
% parametros_Integrador.paso=.001;
% parametros_Integrador.orden=12;
% cantidad_cuerpos_menores=0;

parametros_Integrador.mensaje='Integrando';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%====================================================
funcion_datos= sistema_ficticio(GM);


[pos,vel]=Integrador(funcion,funcion_datos,efemerides_epocas,posicion_ini,velocidad_ini,efemerides_epocas(1),parametros_Integrador);
pos=pos';
l3=length(efemerides_epocas);
l2=funcion_datos.cantidad_planetas;
q=reshape(pos,[3,l2,l3]);
pos=permute(q,[2 1 3]);


vel=vel';
q=reshape(vel,[3,l2,l3]);
vel=permute(q,[2 1 3]);


%%==========Convertir coordenadas a baric�ntricas

% [pos0,posdot0]=helio2bari(pos,posdot0,pos,vel,GM1);
% pos=[pos,pos0];
% vel=[vel,posdot0];
