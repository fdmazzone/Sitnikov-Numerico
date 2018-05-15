
function [pos,vel]=solve2body(t,e)

funcion=@Force2b;


%%%Colocacion
Integrador=@colocacion_nb_adap; 
parametros_Integrador.paso=.001;
parametros_Integrador.tol=1e-25;
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
m=1;
funcion_datos= m;
r=1;
a=r/(1-e);
n=a^(-1.5);
T=2*pi/n;
v_per=sqrt((1+e)/(1-e)/a);
X=[1; 0];
J=[ 0 -1 ;1 0 ];
V=v_per*J*X;


[pos,vel]=Integrador(funcion,funcion_datos,t,X(:)',V(:)',t(1),parametros_Integrador);



%%==========Convertir coordenadas a baric�ntricas

% [pos0,posdot0]=helio2bari(pos,posdot0,pos,vel,GM1);
% pos=[pos,pos0];
% vel=[vel,posdot0];