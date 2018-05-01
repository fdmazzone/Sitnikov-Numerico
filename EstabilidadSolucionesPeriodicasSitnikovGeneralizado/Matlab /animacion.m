function animacion(pos_rot);

figure;
hold on
pos=squeeze(pos_rot(:,:,1));
l=size(pos);
color=rand([l(1),3]);

for j=1:l(1);
    
    h{j}=plot3(pos(j,1),pos(j,2),pos(j,3),...
        '.','MarkerSize',30,'Color',color(j,:));
    %set(h{j},'Erase','xor');
end

long=1.5;
axis([-long long -long long -2.1 2.1])
%plot3(0,0,0,'.','MarkerSize',10,'Color','y')

% r=norm(pos(1,:,1));
% theta=0:.01:2*pi;
% circ=r*[cos(theta);sin(theta);zeros(1,length(theta))];
% plot3(circ(1,:),circ(2,:),circ(3,:))
%  PosPri1=squeeze(pos_rot(1,:,:));
% PosPri2=squeeze(pos_rot(2,:,:));
% plot3(PosPri1(1,:),PosPri1(2,:),PosPri1(3,:))
% plot3(PosPri2(1,:),PosPri2(2,:),PosPri2(3,:))
% % 
% PosPri3=squeeze(pos_rot(3,:,1:50));
% PosPri4=squeeze(pos_rot(4,:,1:50));
% plot3(PosPri3(1,:),PosPri3(2,:),PosPri3(3,:))
% plot3(PosPri4(1,:),PosPri4(2,:),PosPri4(3,:))


set(gca,'Color','k')
axis square
view([45,45])
grid off
hold off
indice_cuadro=1;
%nombrejpg=num2str(indice_cuadro);
%falta=5-length(nombrejpg);

%nombrejpg=[repmat('0',[1,falta]),nombrejpg,'.jpg'];
%eval(['print -djpeg ',nombrejpg]) ;
nombrejpg=['ani-',num2str(indice_cuadro),'.jpg'];

%saveas(gcf,nombrejpg)

indice_cuadro=indice_cuadro+1;


for j=1:size(pos_rot,3);
   drawnow
   pos=squeeze(pos_rot(:,:,j));
   for k=1:l(1);
    set(h{k},'XData',pos(k,1),'YData',pos(k,2),'Zdata',pos(k,3),'Color',color(k,:));
   end
%     
% 
%   M(j) = getframe;
   nombrejpg=['ani-',num2str(indice_cuadro),'.jpg'];
   %saveas(gcf,nombrejpg)
    indice_cuadro=indice_cuadro+1;
end
 %movie(M,1,100);
%movie2avi(M,'JupTroy.avi','compression','Cinepak');