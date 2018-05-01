function AX=MatrizA(X,m)
    l=size(X);
    Ind=nchoosek(1:l(2),2);
    Indicador=[];
    for indi=1:l(2)
        Iaux=find(Ind(:,1)==indi);
        Jaux=find(Ind(:,2)==indi);
        Indicador=[Indicador,[Jaux',Iaux']];
    end
    I=1:l(2)^2;
    I=reshape(I,[l(2),l(2)])';
    for j=2:l(2)
        I(j,:) = [circshift(I(j,1:j),[0,j-1]),I(j,j+1:end)];
    end
    diferencias=X(:,Ind(:,1))-X(:,Ind(:,2)); %se forman todas las diferencias
    r=(sum(diferencias.^2,1).^(-1.5)); %se calcula r_ij^(-3/2)
    AXb=r(Indicador); % [r_12,r_13,....,t_i1, r_21,r_23,.....] en un gran fila
    AXb=diag(m)*reshape(AXb, [ l(2)-1 l(2)])'; %multiplicamos 
    %por masas y agrupamos en en un matriz #cant_cuerpos x #cant_cuerpos-1 
    A=sum(AXb,2);
    AXb=[-A,AXb];
    AXb=AXb';
    AX=AXb(I')'; %permutacion para llevar los elemntos "negativos" a la diagonal
    AX=X+X*AX;
end