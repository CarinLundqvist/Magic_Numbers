%% Density operator 1D
%Calulates the single particle probability density in 1D and a harmonic
%ocillator.

clc
close all
m=1;
w=1;
n_up=6;             %# spin up, n=1 is the groundstate
Psi_up=sym('a', [n_up 1]);
for i=1:n_up
        Psi_up(i,1)=abs(QHO(m,w,i))^2;
end
%Fills the vectors with the probability functions

x=linspace(-10,10,1000);
Prob_up=sum(Psi_up);
P_up=matlabFunction(Prob_up);

figure;
plot(x,P_up(x),'Linewidth',2)
set(gca,'FontSize',20)
lgd=legend('$\langle{\Psi_{\uparrow}}|\hat{n}_{\uparrow}|\Psi_{\uparrow}\rangle $','Interpreter','latex','Location','northeast');
xl=xlabel('$\xi$','Interpreter','latex');
lgd.FontSize=25;
xl.FontSize=35;

