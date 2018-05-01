function ColinealInv(x)
  #x vector fila de posiciones
  l=size(x)[2]
  D=ones(l,1)*x
  C=D-D'
  A=((abs.(C)).^(-3)).*C
  for i=1:l
    A[i,i]=0
  end
  return A\(-x')
end
#x=.01:.01:.99
#y=.01:.01:.99
#A = [ ColinealInv([-1 -x1 -x1*y1 x1*y1 x1 1]) for x1 in x, y1 in y ];
#z=[Float64(~any(k->k<0,l)) for l in A]
