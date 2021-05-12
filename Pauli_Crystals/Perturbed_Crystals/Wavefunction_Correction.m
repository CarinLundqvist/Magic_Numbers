function Psi_1=Wavefunction_Correction(N)
m=1;
w=1;
lim=5;              %Integration limit
Hilbert=20;          %Set size of Hilbert space
x=sym('x', [N 1]);   %Create particle coordinates
Psi_1=0;             % Wavefunction correction
DeltaE=0;            %Denominator in perturbation

[n_up,n_down]=GroundState(N);   %Gives the number of spin up and down particles if N particles are in the groundstate
%Not really necessairy to create a matrix but represent the idea of the code
Occ=zeros(Hilbert,2);         %Possible occupied states
Occ(1:n_up,1)=1;              %Groundstate
Occ(1:n_down,2)=1;            %Groundstate


Psi_up=sym('a', [n_up 1]);
Psi_down=sym('a', [n_down 1]);
slater_matrix_up=sym('b', [n_up n_up]);      %b is a placeholder
slater_matrix_down=sym('b', [n_down n_down]);      %b is a placeholder

%Build wave functions squared
for i=1:n_up
    Psi_up(i,1)=abs(QHO(m,w,i))^2;
end

for i=1:n_down
    Psi_down(i,1)=abs(QHO(m,w,i))^2;
end

%If spin up and spin down differs by one occupation each - AKA thers's is
%one excited particle per spin state
for i=1:n_down
    for j=1:n_up
        excited_down=i;
        excited_up=j;
        for k=(n_down+1):length(Occ)    %length(Occ)=Hilbert
            for l=(n_up+1):length(Occ)
                DeltaE=(k-i)+(l-j);     %E_excited-E_ground - Isn't it the other way around?
                fun=matlabFunction(QHO(m,w,i)*QHO(m,w,k)*QHO(m,w,j)*QHO(m,w,l)); %Integrand
                Correction=integral(fun,-lim,lim);   %Overlap
                %Create slater determinant
                
                %Spin up particles
                for p=1:n_up
                    %Add not excited states to slater
                    if p~=j
                        one_particle_psi=matlabFunction(QHO(m,w,p));
                        for q=1:n_up
                            slater_matrix_up(p,q)=one_particle_psi(x(q));
                        end
                    end
                    %Add the excited state to slater, l here
                    if p==j %DETTA VAR FEL INNAN!!!
                        one_particle_psi=matlabFunction(QHO(m,w,l));
                        for q=1:n_up
                            slater_matrix_up(p,q)=one_particle_psi(x(q));
                        end
                    end
                end
                %Spin down particles
                %Add not excited states to slater
                for p=1:n_down
                    if p~=i
                        one_particle_psi=matlabFunction(QHO(m,w,p));
                        for q=1:n_down
                            slater_matrix_down(p,q)=one_particle_psi(x(n_up+q));
                        end
                    end
                    %Add the excited state to slater, k here
                    if p==i
                        one_particle_psi=matlabFunction(QHO(m,w,k));
                        for q=1:n_down
                            slater_matrix_down(p,q)=one_particle_psi(x(n_up+q));
                        end
                    end
                end
                %Add the corrected wavefunction
                slater_det_up=1/sqrt(factorial(n_up))*det(slater_matrix_up);
                slater_det_down=1/sqrt(factorial(n_down))*det(slater_matrix_down);
                Psi_1=Psi_1-Correction/(DeltaE)*slater_det_up*slater_det_down;
            end
        end
    end
end

%If spin up equals

%Spin up particles
%Add states to slater, will be the same throughout the loop below
for p=1:n_up
    one_particle_psi=matlabFunction(QHO(m,w,p));
    for q=1:n_up
        slater_matrix_up(p,q)=one_particle_psi(x(q));
    end
end
%Slater determinant for spin up
slater_det_up=1/sqrt(factorial(n_up))*det(slater_matrix_up);


for i=1:n_down
    excited_down=i;
    for k=(n_down+1):length(Occ)
        DeltaE=(k-i);     %E_excited-E_ground
        fun=matlabFunction(QHO(m,w,i)*QHO(m,w,k)*(sum(Psi_up))); %Integrand
        Correction=integral(fun,-lim,lim);
        %Create slater determinant
        %Spin down particles
        for p=1:n_down
            %Add not excited states to slater
            if p~=i
                one_particle_psi=matlabFunction(QHO(m,w,p));
                for q=1:n_down
                    slater_matrix_down(p,q)=one_particle_psi(x(n_up+q));
                end
            end
            %Add the excited state to slater, k here
            if p==i
                one_particle_psi=matlabFunction(QHO(m,w,k));
                for q=1:n_down
                    slater_matrix_down(p,q)=one_particle_psi(x(n_up+q));
                end
            end
        end
        %Add the corrected wavefunction
        slater_det_down=1/sqrt(factorial(n_down))*det(slater_matrix_down);
        Psi_1=Psi_1-Correction/(DeltaE)*slater_det_up*slater_det_down;
    end
end


%If spin down equals

%Spin down particles
%Add states to slater, will be the same throughout the loop below
for p=1:n_down
    one_particle_psi=matlabFunction(QHO(m,w,p));
    for q=1:n_down
        slater_matrix_down(p,q)=one_particle_psi(x(n_up+q));
    end
end
%Slater determinant for spin up
slater_det_down=1/sqrt(factorial(n_down))*det(slater_matrix_down);

for j=1:n_up
    excited_up=j;
    for l=(n_up+1):length(Occ)
        DeltaE=(l-j);     %E_excited-E_ground
        fun=matlabFunction(QHO(m,w,l)*QHO(m,w,j)*(sum(Psi_down))); %Integrand
        Correction=integral(fun,-lim,lim);                      %Overlap
        %Create slater determinant
        
        %Spin up particles
        %Add not excited states to slater
        for p=1:n_up
            if p~=j
                one_particle_psi=matlabFunction(QHO(m,w,p));
                for q=1:n_up
                    slater_matrix_up(p,q)=one_particle_psi(x(q));
                end
            end
            %Add the excited state to slater, l here
            if p==j
                one_particle_psi=matlabFunction(QHO(m,w,l));
                for q=1:n_up
                    slater_matrix_up(p,q)=one_particle_psi(x(q));
                end
            end
        end
        %Add the corrected wavefunction
        slater_det_up=1/sqrt(factorial(n_up))*det(slater_matrix_up);
        Psi_1=Psi_1-Correction/(DeltaE)*slater_det_up*slater_det_down;
    end
end

end

