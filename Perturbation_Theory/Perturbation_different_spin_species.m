%% Perturbation different spin species
clc
clf
clear all

m=1;            %Mass
w=1;            %Frequency
alpha=-0.444;     %Interaction strength

N=11;            %Total number of particles

E_int=zeros(N,3);  %Interaction energy
E_0=zeros(N,3);    %unperturbated energy


col = [0, 0, 1;
    0, 1, 0;
    1, 0, 0;
    0, 0, 0;
    1, 0, 1;
    0, 1, 1;
    0.6350 0.0780 0.184];        %Colors


for p=2:4       

for k=1:N
    n_up=0;             %spin up
    n_down=0;           %spin down
    n_left=0;           %third spin
    n_right=0;          %fourth spin
    
    %Build the ground state and ground state energy
        for i=1:k
            if p==2         %Two spin components
                 if mod(i,p) == 0
                    n_down=n_down+1;
                    E_0(k,p-1)=E_0(k,p-1)+1/2+n_down-1;
                 else
                    n_up=n_up+1;
                    E_0(k,p-1)=E_0(k,p-1)+1/2+n_up-1;
                 end
        
        	elseif p==3      %Three spin components
                if mod(i,p) == 2
                    n_down=n_down+1;
                    E_0(k,p-1)=E_0(k,p-1)+1/2+n_down-1;
                elseif mod(i,p) == 1
                    n_up=n_up+1;
                    E_0(k,p-1)=E_0(k,p-1)+1/2+n_up-1;
                else 
                    n_left=n_left+1;
                    E_0(k,p-1)=E_0(k,p-1)+1/2+n_left-1;
                end
            else
                    %Two Four components
                if mod(i,p) == 2
                  n_down=n_down+1;
                  E_0(k,p-1)=E_0(k,p-1)+1/2+n_down-1;
            
                elseif mod(i,p) == 1
                  n_up=n_up+1;
                  E_0(k,p-1)=E_0(k,p-1)+1/2+n_up-1;
            
                elseif mod(i,p) == 3
                 n_left=n_left+1;
                 E_0(k,p-1)=E_0(k,p-1)+1/2+n_left-1;
            
                else
                 n_right=n_right+1;
                 E_0(k,p-1)=E_0(k,p-1)+1/2+n_right-1;
                end
            end
        end
    
    
    Psi_up=sym('a', [n_up 1]);
    Psi_down=sym('a', [n_down 1]);
    Psi_left=sym('a', [n_left 1]);
    Psi_right=sym('a', [n_right 1]);
    
    %Buld wave functions squared
    for i=1:n_up
        Psi_up(i,1)=abs(QHO(m,w,i))^2;
    end
    
    for i=1:n_down
        Psi_down(i,1)=abs(QHO(m,w,i))^2;
    end
    
    for i=1:n_left
        Psi_left(i,1)=abs(QHO(m,w,i))^2;
    end
    for i=1:n_right
        Psi_right(i,1)=abs(QHO(m,w,i))^2;
    end
    
    Psi2_1=Psi_up*Psi_down';  %Matrix containing the integrands
    Psi2_2=Psi_up*Psi_left';  %Matrix containing the integrands
    Psi2_3=Psi_up*Psi_right';  %Matrix containing the integrands
    Psi2_4=Psi_down*Psi_left';  %Matrix containing the integrands
    Psi2_5=Psi_down*Psi_right';  %Matrix containing the integrands
    Psi2_6=Psi_left*Psi_right';  %Matrix containing the integrands
    
    lim=10;             % Because of numerics; good approximation psi=0 outside about abs(x/l_ho)=5
    
    if  p==2        %Two spin components
           for i=1:n_up
                for j=1:n_down
                fun=matlabFunction(Psi2_1(i,j));
                E_int(k,p-1) = E_int(k,p-1)+alpha*integral(fun,-lim,lim);
                end
           end
        
    elseif p==3     %Three spin components
        for i=1:n_up
            for j=1:n_down
              fun=matlabFunction(Psi2_1(i,j));
              E_int(k,p-1) = E_int(k,p-1)+alpha*integral(fun,-lim,lim);
            end
        end
        for i=1:n_up
            for l=1:n_left
             fun=matlabFunction(Psi2_2(i,l));
             E_int(k,p-1) = E_int(k,p-1)+alpha*integral(fun,-lim,lim);
            end
        end
        for j=1:n_down
            for l=1:n_left
            fun=matlabFunction(Psi2_4(j,l));
            E_int(k,p-1) = E_int(k,p-1)+alpha*integral(fun,-lim,lim);
            end
        end
    else
            %Four spin components
             for i=1:n_up
                for j=1:n_down
                fun=matlabFunction(Psi2_1(i,j));
                E_int(k,p-1) = E_int(k,p-1)+alpha*integral(fun,-lim,lim);
                end
             end
            for i=1:n_up
             for l=1:n_left
               fun=matlabFunction(Psi2_2(i,l));
              E_int(k,p-1) = E_int(k,p-1)+alpha*integral(fun,-lim,lim);
             end
            end
            for i=1:n_up
               for l=1:n_right
                fun=matlabFunction(Psi2_3(i,l));
               E_int(k,p-1) = E_int(k,p-1)+alpha*integral(fun,-lim,lim);
               end
            end
                for j=1:n_down
                 for l=1:n_left
                    fun=matlabFunction(Psi2_4(j,l));
                    E_int(k,p-1) = E_int(k,p-1)+alpha*integral(fun,-lim,lim);
                 end
                end
            for j=1:n_down
             for l=1:n_right
            fun=matlabFunction(Psi2_5(j,l));
            E_int(k,p-1) = E_int(k,p-1)+alpha*integral(fun,-lim,lim);
             end
            end
             for j=1:n_left
                 for l=1:n_right
                  fun=matlabFunction(Psi2_6(j,l));
                  E_int(k,p-1) = E_int(k,p-1)+alpha*integral(fun,-lim,lim);
                 end
             end
    
    end
end

end
E_sep=zeros(N,3);     %Separation energy

for j=1:3
for i=1:N
    if i==1
        E_sep(1,j)=E_int(1,j);
    else
        E_sep(i,j)=E_int(i,j)-E_int(i-1,j);
    end
end
end

%Separation energy multiple spin components
n=linspace(1,N,N);
figure(4)
plot(n,E_sep,'-+','Markersize',12,'Linewidth',1.7)
grid on
yl=ylabel('$E_{\mathrm{sep}}\: [\hbar\omega]$','Interpreter','latex');
xl=xlabel('N','Interpreter','latex');
%tl=title('Separation energy as a function of particle number','Interpreter','latex');
set(gca,'FontSize',24)
axis([1 11 -0.665 0])
xl.FontSize=30;
yl.FontSize=30;
%tl.FontSize=22;
lgd=legend('Two spin components','Three spin components','Four spin components','Interpreter','latex','Location','southwest');
lgd.FontSize=23;
 
 
%Energy with and without interaction 
 
n=linspace(0,N,N+1);
E_tot=zeros(N+1,3);

E_tot(2:end,:)=E_int+E_0;    %E_(N): total energy
E_00=zeros(N+1,3);
E_00(2:end,:)=E_0;           %E_(N): Energy without interactions   


for i=1:3
figure(i)
hold on
p(1)=plot(n,E_00(:,i),'-+','Linewidth',2.0,'MarkerSize',14,'Color','r');
hold on
    for j=1:4
        figure(i)
    p(2)=plot([(i+1)*(j-1) (i+1)*j],[E_tot((i+1)*(j-1)+1,i) E_tot((i+1)*(j)+1,i)],'Linewidth',2,'Color','[0.9290, 0.6940, 0.1250]');
    hold on
    end
axis([0 7 0 12.5]);
xticks(0:1:7)
grid on

p(3)=plot(n,E_tot(:,i),'-+','Linewidth',2.0,'MarkerSize',14,'Color','k');

hold on
plot(n,ones(size(n))*0, '-k',n, ones(size(n))*12.5,'-k', 'LineWidth',1)  % Left & Lower Axes
xline(0,'-k', 'LineWidth',2)
xline(7,'-k', 'LineWidth',2)
hold on
set(gca,'FontSize',32)
yl=ylabel('$E_0/\hbar\omega$','Interpreter','latex');
xl=xlabel('N','Interpreter','latex');
tl=title(['$E_0(N)$, ' num2str(i+1), ' Spin components'],'Interpreter','latex');
lgd=legend(p([1 3]),'Without interaction','Interaction present','Interpreter','latex','Location','northwest');

lgd.FontSize=38;
xl.FontSize=38;
yl.FontSize=38;
tl.FontSize=38;

end
