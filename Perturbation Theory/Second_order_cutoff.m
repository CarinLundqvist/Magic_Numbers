%% Cut-off dependancy 
clc
clf
clear all

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
n=linspace(10,90,9);    %To plot E_2
n1=linspace(40,150);    %To plot power fit laws
B=zeros(2,5);           %Store power fit laws for N=2 to N=6                                       %Fitted parameters

for i=2:6
pwrlaw = @(b,n) b(1)+b(2)./n;                        %Data fit E(L)=a+b/L
Diff = @(b) norm(E_2(i,4:end)-pwrlaw(b,n(1,4:end))); %To be minimized, start from L=40
B(:,i-1) = fminsearch(Diff, [-1; 1]);                 %Start with -1,1 and optimize parameters
figure(1)
plot(n,E_2(i,:),'-+','MarkerSize',15,'Linewidth',1.5)
hold on
end
for i=1:5
plot(n1,pwrlaw(B(:,i),n1),'--','Linewidth',1.5,'color','k')
hold on
end


yl=ylabel('$E_0^{(2)}$','Interpreter','latex');
xl=xlabel('L: Highest QHO level included','Interpreter','latex');
  set(gca,'FontSize',24)
xl.FontSize=28;
yl.FontSize=28;
 tl.FontSize=28;
lgd=legend('N=2','N=3','N=4','N=5','N=6','Interpreter','latex','Location','southwest');
lgd.FontSize=24;

%Plot convergance with differences between Ns

% n1=linspace(10,160);
% 
% figure(2)
% for i=1:5
%     if i==1
%         plot(n1,pwrlaw(B(:,1),n1),'-',n,E_2(2,:),'+','Linewidth',1.5)
%     else
% plot(n1,pwrlaw(B(:,i),n1)-pwrlaw(B(:,i-1),n1),'-',n,E_2(i+1,:)-E_2(i,:),'+','Linewidth',1.5)
%     end
% hold on
% end
% lgd=legend('$E_2(L)-E_1(L)$','$E_3(L)-E_2(L)$','$E_4(L)-E_3(L)$','$E_5(L)-E_4(L)$','$E_6(L)-E_5(L)$','Interpreter','latex','Location','southwest');
% lgd.FontSize=10;
% yl=ylabel('$\Delta E^{(2)}/\hbar\omega$','Interpreter','latex');
% xl=xlabel('Highest QHO level included','Interpreter','latex');
%  tl=title('$E_N(L)-E_{N-1}(L)$ as a function of cut-off','Interpreter','latex');
%   set(gca,'FontSize',10)
% xl.FontSize=22;
% yl.FontSize=22;
%  tl.FontSize=22;
%  
%  figure(3)
% for i=2:5
% plot(n1,pwrlaw(B(:,i),n1)-pwrlaw(B(:,1),n1),'-','Linewidth',1.5)
% hold on
% end
% lgd=legend('$E_3(L)-E_2(L)$','$E_4(L)-E_2(L)$','$E_5(L)-E_2(L)$','$E_6(L)-E_2(L)$','Interpreter','latex','Location','southwest');
% lgd.FontSize=20;
% yl=ylabel('$\Delta E^{(2)}/\hbar\omega$','Interpreter','latex');
% xl=xlabel('Highest QHO level included','Interpreter','latex');
%  tl=title('$E_N(L)-E_{2}(L)$ as a function of cut-off','Interpreter','latex');
%   set(gca,'FontSize',25)
% xl.FontSize=32;
% yl.FontSize=32;
%  tl.FontSize=32;