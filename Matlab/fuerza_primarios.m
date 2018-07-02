function estado_aceleracion=fuerza_primarios(t,estado,param)

% % 
% % %estado =  matriz de 1 x 3*cantidad.de.tiempos
% % %param.CC = 3x cantidad.de.primarios
% % %param.m = matriz de 1 x cantidad.de.primarios
% % %param.raizlambda= raiz.de.lambda de la CC.

% 
% t=[0,0.25,0.5,0.75, 1];
% estado=[0.1 0.2 1; 0 1 1; -0.1 -0.2 1.1; 0.2 0.1 1.2; -0.2 -0.1 1.3]';
% param.CC=[1 0 0; -1/2 -sqrt(3)/2 0;-1/2 sqrt(3)/2 0]';
% param.m=[1,2,3];
% param.raizlambda=pi;

t=t';




..................................................................................................................................................................................................
    ...........................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................0.0t=t';

%CALCULAMOS LA POSICION DE LOS PRIMARIOS EN LOS TIEMPOS t_i
lt=size(t,2);
R=[cos(param.raizlambda*t) -sin(param.raizlambda*t) zeros(1,lt) sin(param.raizlambda*t) cos(param.raizlambda*t) zeros(1,lt) zeros(1,lt) zeros(1,lt) ones(1,lt) ];

T=reshape(R,[lt,9])';
R=reshape(T,[3,3*lt])';
posprimarios=R*param.CC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%CALCULAMOS LAS FUERZAS SOBRE LA PART�CILA NO GRAVE
cantprimarios=size(param.CC,2);

estado=estado';
estado1=estado(:);
estado=repmat(estado1,[1,cantprimarios]);

diferencias=posprimarios-estado; % contiene la diferencia entre cada primario y la part�cula no grave, en cada tiempo t_i
diferencias2=diferencias(:)';
%diferencias2=diferencias1';
diferencias3=reshape(diferencias2,[3,lt*cantprimarios]);
diferencias4=diferencias3.^2;
normas2=sum(diferencias4);
normas32=normas2.^(-3/2);
masas=repmat(param.m,lt,1);
masas1=masas(:)';
mnormas=masas1.*normas32; % contiene el producto de las masas por el reciproco de la norma al cubo, en cad tiempo t_i

%ACOMODAMOS LA MATRIZ DE diferencias entre los primarios y la particula no
%grave para luego poder multiplicarla por la matriz de masas x el reciproco
%de la norma al cubo
diferenciasa=reshape(diferencias,[3,lt*cantprimarios]); % 

%%ACOMODAMOS LA MATRIZ DE  de masas x el reciproco
%de la norma al cubo para luego poder multiplicarla por la matriz  de diferencias entre los primarios y la particula no
%grave
mnormas2=reshape(mnormas,[lt,cantprimarios]);
mnormas3=mnormas2(:);
mnormas4=repmat(mnormas3,1,3);
mnormas5=mnormas4';


%%ARMAMOS EL VECTOR DE FUERZAS EN CADA t_i
fuerzast=diferenciasa.*mnormas5;
fuerzast2=reshape(fuerzast,[3*lt,cantprimarios]);
estado_aceleracion=sum(fuerzast2,2)';% vector fila donde est�n las coordenads de las fuerzas aplicads sobre la particula no grave en cada instante t_i

estado_aceleracion=reshape(estado_aceleracion,[lt,3]);









%diferencias7=repmat(diferencias6,[3,1]);
%diferencias8=reshape(diferencias7,[3,cantprimarios,lt]);



% n=size(estado,2)/3;% n�mero de cuerpos
% s_est=size(estado,1);% n�mero de �pocas
% R=permute(reshape(estado',[3,n,s_est]),[2,1,3]);% transformamos la matriz estado en una sucecion de matrices nx3 (hojas),
% %una por cada tiempo, que contienen la posicion de los cuerpos
% sR=size(R);% ya se conoce este dato,sR= [n,3,s_est]
% 
% % calcula la posicion del cuerpo faltante en termino de los anteriores
% % R0=reshape(([-param.GM(2,param.cantidad_planetas+2),param.GM(1,2:param.cantidad_planetas)]/param.GM(1,1))*R(:,:),[1, sR(2:end)]);
% % R=cat(1,R0,R);% agrega el cuerpo faltante a los n-1 dados
% diferencias=R(param.indices(:,1),:,:)-R(param.indices(:,2),:,:);
% norma_diferencias=(sum(diferencias.^2,2).^(-1.5));
% diferencias=repmat(norma_diferencias,[1,3,1]).*diferencias;
% diferencias=diferencias(param.Indicador,:,:);
% sdif=size(diferencias);
% A=-param.GM*diferencias(:,:);
% B=reshape(A,[n,sdif(2:end)]);
% estado_aceleracion=permute(B,[2,1,3]);
% estado_aceleracion=reshape(estado_aceleracion,[3*n,s_est])';