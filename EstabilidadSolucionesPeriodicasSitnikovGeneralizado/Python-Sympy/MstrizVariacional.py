zeta=symbols('zeta')
Phi=Rational(1,4)+zeta**2
F_1=1+Rational(3,4)/Phi**(Rational(5,2))-1/Phi**(Rational(3,2))
F_2=1-1/Phi**(Rational(3,2))
A=Matrix([[0, 0, 1, 0],[0, 0, 0, 1],[F_1, 0, 0, 2],[0, F_2, -2, 0]])
