function [n_up,n_down] = GroundState(N)
%Calculate the ground state of 2 spin components
n_up=0;
n_down=0;
for i=1:N
    if mod(i,2) == 0    
n_down=n_down+1;
    else
        n_up=n_up+1; 
    end
end


end
