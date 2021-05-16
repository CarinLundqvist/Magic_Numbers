function [alpha1 alpha2]=Best_alpha(E_int1,E_int2)

 E_exp=[0 -0.146 -0.1037 -0.212 -0.14844 -0.2798]'; %Experimental data

E_sep=zeros(6,1);                       %1st order seperation energy
E_sep2=zeros(6,1);                      %2nd order seperation energy  

alpha1=-0.5;                            %First guess for 1st order
alpha2=-0.5;                            %First guess for 2nd order

E_1=E_int1*alpha1;                      %1st order energy correction
E_2=E_int1*alpha2+E_int2*(alpha2)^2;    %2nd order energy correction

w1=1;                                   %Weight for 1st order
w2=1;                                   %Weight for 2nd order

%Calculate separation energy
for i=1:6
    if i==1
     E_sep(1,1)=E_1(1,1);
     E_sep2(1,1)=E_2(1,1);
    else
     E_sep(i,1)=E_1(i,1)-E_1(i-1,1);
     E_sep2(1,1)=E_2(i,1)-E_2(i,1);
    end
end
   
%Square difference between calculation and experimental values
diff1=(E_exp-E_sep).^2;                 
diff1(2,1)=w1*diff1(2,1);               %Add weight for N=2

diff2=(E_exp-E_sep2).^2;
diff2(2,1)=w2*diff2(2,1);               %Add weight for N=2

%Sum the square errors
Delta_1_old=sum(diff1);                 %To be minimized 1st order
Delta_2_old=sum(diff2);                 %To be minimized 2nd order

iterations=1e5;         %Number of iterations
lambda=1e-1;            %Start value when incrementing alpha
j=0;                    %Used to decrease |lambda| in iteration

%first order
for k=1:iterations
        alpha1_new=alpha1+lambda;       %New alpha
        E_1=E_int1*alpha1_new;          %New 1st order energy correction
    
        %Calculate separation energy
        for i=1:6
            if i==1
             E_sep(1,1)=E_1(1,1);
            else
             E_sep(i,1)=E_1(i,1)-E_1(i-1,1);
            end
        end
        
        %Square difference between new calculations and experimental values
        diff1=(E_exp-E_sep).^2;
        diff1(2,1)=w1*diff1(2,1);       %Weight
        Delta_1_new=sum(diff1);         %New square error

            if Delta_1_new<Delta_1_old
                alpha1=alpha1_new;          %Assign alpha1 the new value
                Delta_1_old=Delta_1_new;    %New square error replaces old
            else
                lambda=-lambda;             %If new alpha worse, change direction
                j=j+1;                      %Count number of times this occurs
            end
                if j>2
                    lambda=0.99*lambda;     %If direction is changed more than twice, decrease increment of alpha
                    j=0;                    %Reset j
                end
end

lambda=1e-1;            %Increment alpha
j=0;                    %Used to decrease |lambda| in iteration

%2nd order 
for k=1:iterations
        alpha2_new=alpha2+lambda;                       %New alpha
        E_2=E_int1*alpha2_new+E_int2*(alpha2_new)^2;    %New 2nd order energy correction
    
        %Calculate separation energy
        for i=1:6
            if i==1
             E_sep2(1,1)=E_2(1,1);
            else
             E_sep2(i,1)=E_2(i,1)-E_2(i-1,1);
            end
        end
        
        %Square difference between new calculations and experimental values
        diff2=(E_exp-E_sep2).^2;
        diff2(2,1)=w2*diff2(2,1);           %weight
        Delta_2_new=sum(diff2);             %New square error

            if Delta_2_new<Delta_2_old
                alpha2=alpha2_new;           %Assign alpha2 the new value
                Delta_2_old=Delta_2_new;     %New square error replaces old
            else
                lambda=-lambda;             %If new alpha worse, change direction
                j=j+1;                      %Count number of times this occurs
            end
                if j>2
                    lambda=0.99*lambda;     %If direction is changed more than twice, decrease increment of alpha
                    j=0;
                end
end
end