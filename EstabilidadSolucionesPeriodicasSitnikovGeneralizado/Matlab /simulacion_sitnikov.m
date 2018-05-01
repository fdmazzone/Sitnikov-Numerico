function [pos,vel]=simulacion_sitnikov(primarios_config,efemerides_epocas,...
    posicion_ini,velocidad_ini)
%%%%%%%%%%%Ejemplo de Uso%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % m=[2,2,2];
% m=ones(1,3);
%% X0=poligono(15);
% X0=[1 0; -sqrt(3)/2 -1/2;-sqrt(3)/2 +1/2]';
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

funcion=@fuerza_primarios;


%%%Colocacion
Integrador=@colocacion_nb_adap_noaut; 
parametros_Integrador.paso=.01;
parametros_Integrador.tol=5e-22;
parametros_Integrador.orden=12;
parametros_Integrador.iter=3;

%multipaso
% Integrador=@multipaso_nb_implicito_noaut;
% parametros_Integrador.paso=.001;
% parametros_Integrador.orden=12;
% cantidad_cuerpos_menores=0;

parametros_Integrador.mensaje='Integrando';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%====================================================





[pos,vel]=Integrador(funcion,primarios_config,efemerides_epocas,posicion_ini,velocidad_ini,efemerides_epocas(1),parametros_Integrador);
pos=pos';



vel=vel';



%%==========Convertir coordenadas a baric�ntricas

% [pos0,posdot0]=helio2bari(pos,posdot0,pos,vel,GM1);
% pos=[pos,pos0];
% vel=[vel,posdot0];
