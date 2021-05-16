%% Kandidatarbete
clc
clf
y=@(x) 0.5*x.^2;           %Harmonic potental
u=@(x) 3-(x-3).^2;         %Tilted potential
x1=linspace(-3,2);
x2=linspace(2,4.5);

plot(x1,y(x1),'k',x2,u(x2),'k','Linewidth',2.0)
axis([-3 4.5 0 4.5])
hold on
%Draw arrows
        pos = get(gca, 'Position'); 
        xPosition = 0.5;        
        yPosition = 2;
        slut=(yPosition-min(ylim))/diff(ylim)*pos(4)+pos(2);
        topp=(yPosition+0.4-min(ylim))/diff(ylim)*pos(4)+pos(2);
        annotation('textarrow',[0.5 0.5],[slut topp],'FontSize',10,'Linewidth',2)
        annotation('textarrow',[0.88 0.88],[topp slut],'FontSize',10,'Linewidth',2)
grid on
set(gca,'FontSize',15)
yl=ylabel('$E/\hbar \omega$','Interpreter','latex');
xl=xlabel('$\xi$','Interpreter','latex');
%tl=title('Tilted QHO Potential $V(\xi)$','Interpreter','latex');
xl.FontSize=32;
yl.FontSize=32;
%tl.FontSize=30;

%% Visually representing the QHO.
clc
clf
clear all
m=1;                % Mass 
w=1;                %Angular frequency
hbar=1;             %Plancks constant

n_up=4;                   % # spin up, n=1 is the groundstate
n_down=0;                 % # spin down
y=@(x) 0.5*x.^2;           %Harmonic potental
Psi_n_up=cell(n_up,1);     %spin up wavefunction
Psi_n_down=cell(n_down,1); %spin down wavefunction
for i=1:n_up
    Psi_n_up{i,1}=matlabFunction(QHO(m,w,i));
end
for i=1:n_down
    Psi_n_down{i,1}=matlabFunction(QHO(m,w,i));
end
x=linspace(-4,4,10000);
E=ones(size(x));    %Energy
figure(1)           %Plot occupation 
for j=1:n_up
    if j<=n_down 
         plot(x,abs(Psi_n_up{j,1}(x)).^2+j-1/2,x,E*(j-1/2),'k',x,y(x),'k','Linewidth',2.0)
         hold on
         plot(x,E*(j-1/2),'k',x,y(x),'k','Linewidth',2.0)
         hold on
         axis([-4 4 0 n_up])
         set(gca,'FontSize',15)
         pos = get(gca, 'Position');
         xPosition = 0.5;
         yPosition = j-1/2;
         slut=(yPosition-min(ylim))/diff(ylim)*pos(4)+pos(2);
         topp=(yPosition+0.3-min(ylim))/diff(ylim)*pos(4)+pos(2);
         annotation('textarrow',[0.5 0.5],[slut topp],'FontSize',10,'Linewidth',2)
         annotation('textarrow',[0.53 0.53],[topp slut],'FontSize',10,'Linewidth',2)
    else
          plot(x,abs(Psi_n_up{j,1}(x)).^2+j-1/2,x,E*(j-1/2),'k',x,y(x),'k','Linewidth',2.0)
          hold on
          plot(x,E*(j-1/2),'k',x,y(x),'k','Linewidth',2.0)
          hold on
          axis([-4 4 0 n_up])
          set(gca,'FontSize',20)
          pos = get(gca, 'Position');
          xPosition = 0.5;
          yPosition = j-1/2;
          slut=(yPosition-min(ylim))/diff(ylim)*pos(4)+pos(2);
          topp=(yPosition+0.3-min(ylim))/diff(ylim)*pos(4)+pos(2);
          annotation('textarrow',[0.5 0.5],[slut topp],'FontSize',10,'Linewidth',2)
    end
end
hold on
grid on
 set(gca,'FontSize',22)
yl=ylabel('$E/\hbar \omega$','Interpreter','latex');
xl=xlabel('$\xi$','Interpreter','latex');
 tl=title('Quantum harmonic oscillator','Interpreter','latex');
% tl=title('Occupied 1D-QHO states','Interpreter','latex');
xl.FontSize=26;
yl.FontSize=30;
tl.FontSize=30;

%% Probaility density
clc
clear all
clf
m=1;
w=1;
n_up=8;             %spin up, n=1 is the groundstate
n_down=5;           %spin down
Psi_up=sym('a', [n_up 1]);
Psi_down=sym('a', [n_down 1]);
for i=1:n_up
        Psi_up(i,1)=abs(QHO(m,w,i))^2;
end

for i=1:n_down
        Psi_down(i,1)=abs(QHO(m,w,i))^2;
end

x=linspace(-10,10,1000);
Prob_up=sum(Psi_up);
P_up=matlabFunction(Prob_up);

Prob_down=sum(Psi_down);
P_down=matlabFunction(Prob_down);

figure(2)
plot(x,P_up(x),x,P_down(x),'Linewidth',2)
set(gca,'FontSize',20)
lgd=legend('$\langle{\Psi_{\uparrow}}|\hat{n}_{\uparrow}|\Psi_{\uparrow}\rangle $','$\langle{\Psi_{\downarrow}}|\hat{n}_{\downarrow}|\Psi_{\downarrow}\rangle $','Interpreter','latex','Location','northeast');
xl=xlabel('$\xi$','Interpreter','latex');
lgd.FontSize=25;
xl.FontSize=35;
%% Separation energy
clc
clf
clear all

m=1;            %Mass
w=1;            %Frequency
N=6;            %Total number of particles

E_int=zeros(N,1);    %Interaction energy 
E_0=zeros(N,1);      %unperturbated energy
E_sep=zeros(N,1);    %Separation energy 1st order
E_sep2=zeros(N,1);   %Separation energy 2nd order

for k=1:N
n_up=0;             %spin up, n=1 is the groundstate
n_down=0;           %spin down

%Build the ground state
for i=1:k
    if mod(i,2) == 0    
n_down=n_down+1;
E_0(k)=E_0(k)+1/2+n_down-1;      %Ground state energy
    else
        n_up=n_up+1; 
        E_0(k)=E_0(k)+1/2+n_up-1; %Ground state energy
    end
end

Psi_up=sym('a', [n_up 1]);       
Psi_down=sym('a', [n_down 1]);

%Build wave functions squared
for i=1:n_up
        Psi_up(i,1)=abs(QHO(m,w,i))^2; 
end

for i=1:n_down
        Psi_down(i,1)=abs(QHO(m,w,i))^2;
end

Psi2=Psi_up*Psi_down';  %Matrix containing the integrands
lim=10;                 % Because of numerics; good approximation psi=0 outside about abs(xi)=5

%Calculate E_int
for i=1:n_up
    for j=1:n_down
        fun=matlabFunction(Psi2(i,j));
 E_int(k) = E_int(k)+integral(fun,-lim,lim);
    end
end
end
  
%2nd order energy correction
E_2=zeros(6,9);
E_2_10=[0 -0.078547750268004 -0.068644384192299 -0.119344771160209 -0.106317920988807 -0.142817317979524]'; %Hilbert space 10
E_2(:,1)=E_2_10;
E_2_20=[0 -0.0886831273738413 -0.0850096889575458 -0.150318331094729 -0.146257091255700 -0.196904396236936]'; %Hilbert space 20
E_2(:,2)=E_2_20;
E_2_30=[0 -0.0929236700919974 -0.0917089283867832 -0.163005254466902 -0.163363112096064 -0.221638636800787]'; %Hilbert space 30
E_2(:,3)=E_2_30;
E_2_40=[0 -0.0953864629020749 -0.0955583492894288 -0.170237199633526 -0.173035455022967 -0.235578842776610]'; %Hilbert space 40
E_2(:,4)=E_2_40;
E_2_50=[0 -0.0970409250861239 -0.0981284685547293 -0.175044327852592 -0.179431644131896 -0.244757144089653]'; %Hilbert space 50
E_2(:,5)=E_2_50;
E_2_60=[0 -0.0982491001201812 -0.0999976186403667 -0.178530394930852 -0.184054459641103 -0.251372140676776]'; %Hiilbert space 60
E_2(:,6)=E_2_60;
E_2_70=[0 -0.0991806398911389 -0.101434504449601 -0.181204872843186 -0.187592500390182 -0.256424968511996]'; %Hilbert space 70
E_2(:,7)=E_2_70;
E_2_80=[0 -0.100028730112028 -0.102650264168645 -0.183444628067795 -0.190477356699716 -0.260525565701513]'; %Hilbert 80
E_2(:,8)=E_2_80;
E_2_90=[0 -0.100541961634859 -0.103527811195704 -0.185093219466811 -0.192723626430964 -0.263738610019445]'; %Hilbert space 90
E_2(:,9)=E_2_90;

E_2_fit=[0 -0.1046   -0.1097   -0.1967   -0.2081   -0.2857]';       %Extrapolated values
E_2(:,9)=E_2_fit;

%Find best alpha
[alpha1 alpha2]=Best_alpha(E_int,E_2(:,9))

%Used in Exact diagonalization (ED)
alpha11=-0.45;
alpha22=-0.45;

E_sep11=zeros(N,1);    %Separation energy 1st order alpha=-0.45
E_sep22=zeros(N,1);    %Separation energy 2nd order alpha=-0.45

%Calculate separation energy
for i=1:N
    if i==1
   E_sep(1,1)=alpha1*E_int(1,1);
   E_sep11(1,1)=alpha11*E_int(1,1);
   E_sep2(1,1)=alpha2*E_int(1,1)+(alpha2)^2*E_2(1,9);
   E_sep22(1,1)=alpha22*E_int(1,1)+(alpha22)^2*E_2(1,9);
    else
   E_sep(i,1)=alpha1*(E_int(i,1)-E_int(i-1,1));
   E_sep11(i,1)=alpha11*(E_int(i,1)-E_int(i-1,1));
   E_sep2(i,1)=(alpha2*E_int(i,1)+(alpha2)^2*E_2(i,9))-(alpha2*E_int(i-1,1)+(alpha2)^2*E_2(i-1,9));
   E_sep22(i,1)=(alpha22*E_int(i,1)+(alpha22)^2*E_2(i,9))-(alpha22*E_int(i-1,1)+(alpha22)^2*E_2(i-1,9));
    end
end



E_exp=[0 -0.146 -0.1037 -0.212 -0.14844 -0.2798];              %Experimental results
err=[0 -0.009636 -0.0123403 -0.018868 -0.02315664 -0.015389];  %Accompanied error
E_ED=[0 -0.19885135 -0.08874585 -0.2389512 -0.1449866 -0.2724738]'; %Calculated with ED
n=linspace(1,N,N);

figure(3)
hold on
   p(4)=plot(n,E_ED,'o','MarkerSize',12,'Linewidth',3,'Color','[0.4660, 0.8740, 0.1880]	');
   p(1)=plot(n,E_sep,'-_','color','[0, 0.4470, 0.7410]','MarkerSize',18,'Linewidth',1.7);
   p(2)=plot(n,E_sep2,'-_','color','[0.8500, 0.3250, 0.0980','MarkerSize',18,'Linewidth',1.7);
   p(3)=errorbar(1:1:6,E_exp,err,'Color','k','Linestyle','none','Linewidth',2.5);
   p(5)=plot(n,E_sep22,'x','MarkerSize',15,'color','[0.8500, 0.3250, 0.0980','Linewidth',2);
   %  p(6)=plot(n,E_sep11,'x','MarkerSize',15,'color','b','Linewidth',2);
grid on

yl=ylabel('$E_{\mathrm{sep}}\: [\hbar\omega]$','Interpreter','latex');
xl=xlabel('N','Interpreter','latex');
% tl=title('Separation energy as a function of particle number','Interpreter','latex');
  set(gca,'FontSize',24)
xl.FontSize=30;
yl.FontSize=30;
% tl.FontSize=28;
   lgd=legend([p(1) p(2) p(3) p(4)],'First order perturbation theory','Second order perturbation theory','Zuern et al. [2] ','D''Amico et al. [17]','Interpreter','latex','Location','northeast');
lgd.FontSize=23;
%% Comparing analytical solution with perturbation

E_tot_1st= @(alpha) E_0(2)+alpha*E_int(2); %E_(2) with 1st order perturbation theory;
E_tot_2nd= @(alpha) E_0(2)+alpha*E_int(2)+(alpha).^2*E_2_fit(2);  %E_(2) with 2nd order perturbation theory;
alpha=linspace(-1,1);               %Interaction strength            
beta = @(E) -2.*sqrt(2).*gamma(-E/2+3/4)./gamma(-E/2+1/4);  %Analytical interaction strength
Es=linspace(-1,1);          %Energy
alphas=beta(Es);            %Interaction strength as a function of Energy

figure(2)
set(gca,'FontSize',24)
plot(alphas,Es+0.5,'color','[0.4660 0.6740 0.1880]','LineWidth',1.5)
hold on
plot(alpha,E_tot_1st(alpha),'b',alpha,E_tot_2nd(alpha),'r','LineWidth',1.5)
hold on 
axis([-1 1 0.4 1.4])
xticks(-1:0.2:1)
grid on
yl=ylabel('$E_0/\hbar\omega$','Interpreter','latex');
xl=xlabel('$\alpha/(l_{ho}\hbar\omega)$','Interpreter','latex');
% tl=title('Energy as a function of interaction strength','Interpreter','latex');
   set(gca,'FontSize',24)
xl.FontSize=28;
yl.FontSize=28;
% tl.FontSize=27;
  lgd=legend('Analytical solution','First order perturbation theory','Second order perturbation theory','Interpreter','latex','Location','southeast');
lgd.FontSize=28;

