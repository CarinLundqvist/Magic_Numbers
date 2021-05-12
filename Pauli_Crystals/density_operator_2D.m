%% Probability density 2D
%Calculates the single particle probability density in 2D for 1 ,3 6 and 10
%particles in an harmonic oscillaotor 

clc
close all
m=1;
w=1;
N=1;             %Total number of paticles - Only complete shells: 1,3,6, 10
Psi=sym('a', [N 1]);      %a is a placeholder
%Creates empty "vectors" for the one particle wave functions with a and x
%as variables
for i=1:N
    if i==1
        n_x=1;
        n_y=1;
    elseif i==2
        n_x=2;
        n_y=1;
    elseif i==3
        n_x=1;
        n_y=2;
    elseif i==4
        n_x=3;
        n_y=1;
    elseif i==5
        n_x=2;
        n_y=2;
    elseif i==6
        n_x=1;
        n_y=3;
    elseif i==7
        n_x=4;
        n_y=1;
    elseif i==8
        n_x=3;
        n_y=2;
    elseif i==9
        n_x=2;
        n_y=3;
    elseif i==10
        n_x=1;
        n_y=4;
    end
    Psi(i,1)=abs(QHO_2D(m,w,n_x,n_y))^2;
end
%Fills the vectors with the probability functions

edge=3;
datapunkter=1e2;
x=linspace(-edge,edge,datapunkter);
y=linspace(-edge,edge,datapunkter);
Prob=sum(Psi);
P=matlabFunction(Prob);

figure
[X,Y]=meshgrid(x,y);
surf(X,Y,P(X,Y))
set(gca,'FontSize',20)
axis equal
xticklabels({-2 0 2})
yticklabels({-2 0 2})
view(0,90)
colorbar
shading interp
%colormap(jet)

%% Plots the single-particle wave functions in 2D
num=3; %Particle number
WF=matlabFunction(Psi(num,1));

figure;
surf(X,Y,WF(X,Y))