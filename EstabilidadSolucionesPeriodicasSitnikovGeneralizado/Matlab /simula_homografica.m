L=1.5;
e=0.9;
m2=4*(1+L^2)^(3/2)*(8*L^3-(1+L^2)^(3/2))/(64*L^3-(1+L^2)^3);
m1=4*L^3*(1+L^2)^(3/2)*(8-(1+L^2)^(3/2))/(64*L^3-(1+L^2)^3);
X=[L 0 -L 0;0 1 0 -1;0 0 0 0];
m=[m1 m2 m1 m2];
n=(1-e)^1.5;
T=2*pi/n;
efemerides_epocas=(0:.05:10*T)';
J=[ 0 -1 0;1 0 0;0 0 0];
V=sqrt(1+e)*J*X;
[pos,vel]=simulacion_sistema_ficticio(m,efemerides_epocas,X(:)',V(:)');
cue=squeeze(pos(1,1,:));
cos_test=cos(efemerides_epocas*2*pi/T);
plot(efemerides_epocas,[cue,cos_test])