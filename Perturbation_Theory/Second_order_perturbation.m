%% Second order energy correction different Hilbert
clear all
clc
format long
m=1;                 %Mass
w=1;                 %Frequency
N=6;                 % Number of particles 
E_n2=zeros(N-1,8);   %Second order energy correction
lim=14;              %Integration limit
Hilbert=10;          %Set size of Hilbert space  

for h=1:1
    Hilbert=10*h;
for u=2:N
[n_up,n_down]=GroundState(u);

%Not really necessairy to create a matrix but represent the idea of the code
Occ=zeros(Hilbert,2);         %Possible occupied states
Occ(1:n_up,1)=1;              %Groundstate
Occ(1:n_down,2)=1;            %Groundstate


Psi_up=sym('a', [n_up 1]);      %Create wave function spin up      
Psi_down=sym('a', [n_down 1]);  %Create wave function spin down

%Build wave functions squared
for i=1:n_up
        Psi_up(i,1)=abs(QHO(m,w,i))^2; 
end

for i=1:n_down
        Psi_down(i,1)=abs(QHO(m,w,i))^2;
end

DeltaE=0;            %Denominator in perturbation

%If spin up and spin down differs by one occupation each
for i=1:n_down
    for j=1:n_up
        excited_down=i;
        excited_up=j;
        X = ['First for-loop: ex down = ',num2str(i),', ex up = ',num2str(j)];  %Keep track of which loop it is
        disp(X)
            for k=(n_down+1):length(Occ)
                for l=(n_up+1):length(Occ)
                    DeltaE=(k-i)+(l-j);     %E_excited-E_ground
                        fun=matlabFunction(QHO(m,w,i)*QHO(m,w,k)*QHO(m,w,j)*QHO(m,w,l)); %Integrand
                        E_n2(u-1,h)=E_n2(u-1,h)-1/(DeltaE)*(integral(fun,-lim,lim))^2; %Overlap
                end
            end
    end
end

%If spin up equals
for i=1:n_down
        excited_down=i;
        X = ['Second for-loop: ex down = ',num2str(i)];  %Keep track of which loop it is
        disp(X)
            for k=(n_down+1):length(Occ)
                    DeltaE=(k-i);     %E_excited-E_ground
                    fun=matlabFunction(QHO(m,w,i)*QHO(m,w,k)*(sum(Psi_up)));           %Integrand
                    E_n2(u-1,h)=E_n2(u-1,h)-1/(DeltaE)*(integral(fun,-lim,lim))^2;     %Overlap
            end
end

%If spin down equals
for j=1:n_up
        excited_up=j;
        X = ['Third for-loop: ex up = ',num2str(j)];  %Keep track of which loop it is
        disp(X)
            for l=(n_up+1):length(Occ)
                    DeltaE=(l-j);     %E_excited-E_ground
                    fun=matlabFunction(QHO(m,w,l)*QHO(m,w,j)*(sum(Psi_down)));           %Integrand
                    E_n2(u-1,h)=E_n2(u-1,h)-1/(DeltaE)*(integral(fun,-lim,lim))^2;       %Overlap
            end
end
end
end
