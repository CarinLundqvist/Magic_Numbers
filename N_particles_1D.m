function Psi_N_particles=N_particles_1D(N)
%Takes 1-10 particles and creates the many-body fermionic wave
%function using the Slater determinant
x=sym('x', [N 1]);

m=1;   
w=1;

slater_matrix=sym('b', [N N]);      %b is a placeholder
for i=1:N
    one_particle_psi=matlabFunction(QHO(m,w,i));
    for j=1:N
       slater_matrix(i,j)=one_particle_psi(x(j)); 
    end
end
slater_det=1/sqrt(factorial(N))*det(slater_matrix);
Psi_N_particles=slater_det;
end