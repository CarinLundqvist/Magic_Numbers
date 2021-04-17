function Psi = QHO(m,w,n)
syms x
 hbar=1;                     %
% hbar=1.0545718e-34;         %Plank's constant

f=factorial(n-1);                         %n!
psi_0=(m*w/(pi*hbar))^(1/4)*exp(-(m*w*x^2)/(2*hbar));  %Ground state
if n==1
    Psi= psi_0;
else
    Dn = diff(exp(-m*w*x^2/hbar),n-1); 
    H_n=(-1)^(n-1)*exp(m*w*x^2/hbar)*Dn;
    Psi=(m*w/(pi*hbar))^(1/4)*1/(sqrt(2^(n-1)*f))*H_n*exp(-m*w*x^2/(2*hbar)); %Psi_n
end
end