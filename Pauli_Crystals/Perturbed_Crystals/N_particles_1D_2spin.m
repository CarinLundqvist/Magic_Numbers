function Psi_N_particles=N_particles_1D_2spin(N)
%Calculates the many-body wave function with two spin components. 

x=sym('x', [N 1]);  %Symbolic positions of the N particles

m=1;   
w=1;
[n_up,n_down]=GroundState(N);   %Number of particles in the spin up and down states

slater_matrix_up=sym('b', [n_up n_up]);      %b is a placeholder
slater_matrix_down=sym('b', [n_down n_down]);      %b is a placeholder
for i=1:n_up
    one_particle_psi=matlabFunction(QHO(m,w,i));    %All single particle states for spin up particles
    for j=1:n_up
       slater_matrix_up(i,j)=one_particle_psi(x(j));    %A Slater matrix with all the spin up states and postions
    end
end
for i=1:n_down
    one_particle_psi=matlabFunction(QHO(m,w,i));    %All single particle states for spin down particles
    for j=1:n_down
       slater_matrix_down(i,j)=one_particle_psi(x(n_up+j)); %A Slater matrix with all the spin down states and postions
    end
end
%The complete many-body wave function, formed by multiplying the spin up and down components 
slater_det_up=1/sqrt(factorial(n_up))*det(slater_matrix_up);
slater_det_down=1/sqrt(factorial(n_down))*det(slater_matrix_down);
Psi_N_particles=slater_det_up*slater_det_down;
end