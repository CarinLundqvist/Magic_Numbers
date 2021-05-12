function d_tot=pauli_distance(N,pos_cm)
%Transforms the particle positions from cartesian to polar coordinates
pos_pol=zeros(N,2);
[theta,rho]=cart2pol(pos_cm(:,1),pos_cm(:,2));
pos_pol(:,1)=theta;
pos_pol(:,2)=rho;

%Defines the Pauli crystal that the particles are aligned with.
if N==3
    pauli_pos=[2*pi/3,0.816;0,0.816;-2*pi/3,0.816];
elseif N==6
    pauli_pos=[4*pi/5,1.265;2*pi/5,1.265;0,1.265;-2*pi/5,1.265;-4*pi/5,1.265];
    %Removes the particle in the middle.
    [min_rho,index]=min(pos_pol(:,2));
    pos_pol(index,:)=[];
    N=5;
end
%Sorts the angles in descending order as to match pauli_pos.
pos_pol(:,1)=sort(pos_pol(:,1),'descend');

%Makes sure the smallest angle between pauli_pos and pos_pol is used to
%calculate d.
for i=1:N
    if pos_pol(i,1)-pauli_pos(i,1)<-pi
        pos_pol(i,1)=pos_pol(i+1)+2*pi;
    elseif pos_pol(i,1)-pauli_pos(i,1)>pi
        pos_pol(i,1)=pos_pol(i,1)-2*pi;
    end
end

%Assigns evey particle in a measurement an equvalent position for the ideal
%Pauli crystal. It then checks the difference squared between the angles of
%the measurement and the ideal postion for each particle and sums them up.
%If d is minimised, the measurments should be rotated in the way that most
%resemble the Pauli crystal.
d=(pos_pol(:,1)-pauli_pos(:,1)).^2;
d_tot=sum(d);

end