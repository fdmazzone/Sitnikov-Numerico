function Estructura= sistema_ficticio(GM1);
  





%%%% Definicion de parametros del sistema para la funcion fuerza

NroCuerpos=length(GM1);


I=1:1:NroCuerpos;
pares=nchoosek(I,2);


Indicador=[];
for indi=1:NroCuerpos
    Iaux=find(pares(:,1)==indi);
    Jaux=find(pares(:,2)==indi);
    Indicador=[Indicador,[Jaux',Iaux']];
end




GM=zeros( NroCuerpos,NroCuerpos*(NroCuerpos-1));
columna=1;
GM(1,columna:columna+NroCuerpos-2)=GM1(2:end);
columna=columna+NroCuerpos-1;
for j=2:NroCuerpos-1
    GM(j,columna:columna+NroCuerpos-2)=[-GM1(1:j-1),GM1(j+1:end)];
    columna=columna+NroCuerpos-1;
end
GM(NroCuerpos,columna:columna+NroCuerpos-2)=-GM1(1:end-1);






Estructura.GM=sparse(GM);

Estructura.cantidad_planetas=NroCuerpos;
Estructura.Indicador=Indicador;
Estructura.indices=pares;