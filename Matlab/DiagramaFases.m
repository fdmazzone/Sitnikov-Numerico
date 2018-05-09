 
e=0.1;
r=1;
a=r/(1-e);
n=a^(-1.5);
T=2*pi/n;
v_per=sqrt((1+e)/(1-e)/a);

 X=[-.5,.5,0;0,0,0;0,0,0];
 J=[ 0 -1 0;1 0 0;0 0 0];
 m=.5*[1,1,0];
 V=J*X;
 V(2,1:2)=v_per*V(2,1:2);
 n=17;
 figure;
 hold on;
 efemerides_epocas=(0:T:1000*T)';
 for j=.1:.01:2.05
     
     X(3,3)=j;
     [pos,vel]=simulacion_sistema_ficticio(m,efemerides_epocas,X(:)',V(:)');
     part=squeeze(pos(3,3,:));
     vel_part=squeeze(vel(3,3,:));
     plot(part,vel_part,'.') 
 end