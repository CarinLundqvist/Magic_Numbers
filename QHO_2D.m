function Psi = QHO_2D(m,w,n_x,n_y)
syms x y
 hbar=1;                     %
% hbar=1.0545718e-34;         %Plank's constant

f_x=factorial(n_x-1);                         %n!
f_y=factorial(n_y-1);
psi_0x=(m*w/(pi*hbar))^(1/8)*exp(-(m*w*x^2)/(2*hbar));  %Ground states
psi_0y=(m*w/(pi*hbar))^(1/8)*exp(-(m*w*y^2)/(2*hbar));

%In the x-dirextion
if n_x==1
    Psi_x= psi_0x;
else
Dn = diff(exp(-m*w*x^2/hbar),n_x-1); 
    H_n=(-1)^(n_x-1)*exp(m*w*x^2/hbar)*Dn;
    Psi_x=(m*w/(pi*hbar))^(1/8)*1/(sqrt(2^(n_x-1)*f_x))*H_n*exp(-m*w*x^2/(2*hbar)); %Psi_n
end

%In the y-direction
if n_y==1
    Psi_y= psi_0y;
else
Dn = diff(exp(-m*w*y^2/hbar),n_y-1); 
    H_n=(-1)^(n_y-1)*exp(m*w*y^2/hbar)*Dn;
    Psi_y=(m*w/(pi*hbar))^(1/8)*1/(sqrt(2^(n_y-1)*f_y))*H_n*exp(-m*w*y^2/(2*hbar)); %Psi_n
end
%Combining both directions
Psi=Psi_x*Psi_y;
end