function Psi_N_particles=N_particles_2D(N)
%Takes 1, 3 and 6 particles and creates the many-body fermionic wave
%function using the Slater determinant

%Creates N varibles for postions in both the x and y direction 
x=sym('x', [N 2]);
m=1;   
w=1;

%Constructs all the possible states for 10 particles in the ground state of
%a 2D harmonic oscillator. Even the degenerate ones.
one_particle_psi=sym('b', [N 1]);      %b is a placeholder
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
    one_particle_psi(i,1)=QHO_2D(m,w,n_x,n_y);
end

%Creates the Slater matrix using the just constructed single particle wave
%functions and inserting the position variables
slater_matrix=sym('b', [N N]);      %b is a placeholder
for j=1:N
    one_psi=matlabFunction(one_particle_psi(j,1));
    for k=1:N
       slater_matrix(j,k)=one_psi(x(k,1),x(k,2)); 
    end
end
%Makes the Slater determinant and adds a normalization constant
slater_det=1/sqrt(factorial(N))*det(slater_matrix);
Psi_N_particles=slater_det;
end