function estado_aceleracion=fuerza_nb(estado,param)
n=size(estado,2)/3;% número de cuerpos
s_est=size(estado,1);% número de épocas
R=permute(reshape(estado',[3,n,s_est]),[2,1,3]);% transformamos la matriz estado en una sucecion de matrices nx3 (hojas),
%una por cada tiempo, que contienen la posicion de los cuerpos
sR=size(R);% ya se conoce este dato,sR= [n,3,s_est]

% calcula la posicion del cuerpo faltante en termino de los anteriores
% R0=reshape(([-param.GM(2,param.cantidad_planetas+2),param.GM(1,2:param.cantidad_planetas)]/param.GM(1,1))*R(:,:),[1, sR(2:end)]);
% R=cat(1,R0,R);% agrega el cuerpo faltante a los n-1 dados
diferencias=R(param.indices(:,1),:,:)-R(param.indices(:,2),:,:);
norma_diferencias=(sum(diferencias.^2,2).^(-1.5));
diferencias=repmat(norma_diferencias,[1,3,1]).*diferencias;
diferencias=diferencias(param.Indicador,:,:);
sdif=size(diferencias);
A=-param.GM*diferencias(:,:);
B=reshape(A,[n,sdif(2:end)]);
estado_aceleracion=permute(B,[2,1,3]);
estado_aceleracion=reshape(estado_aceleracion,[3*n,s_est])';

 