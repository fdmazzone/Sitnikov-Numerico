function X0 = poligono(n)
ang=2*pi/n;
R=[cos(ang) -sin(ang);sin(ang) cos(ang)];
X0=[1;0];
for k=2:n
    X0=[X0,R*X0(:,end)];
end

end

