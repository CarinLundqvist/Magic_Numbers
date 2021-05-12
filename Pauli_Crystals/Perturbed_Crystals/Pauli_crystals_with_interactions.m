% Pauli crystals with interactions
%% Creates unperturbed wave function and correction in 1D
%Run once for a certain N
clc
N=2;    %Number of particles 1-10
[n_up,n_down]=GroundState(N);   %Number of spin up and down particles in the ground state

Psi_0=N_particles_1D_2spin(N);      %Constructs the unperturbed many-body wave function
Psi_1=Wavefunction_Correction(N);   %First order correction to the wave function

%% Combines ground state and correction
%Run once after running the code above

%Select an appropriate interaction strength alpha. Generally with a magnitude less than 1
alpha=-0.5;
Prob=abs(Psi_0+alpha*Psi_1)^2;      %The perturbed probability distribution. Set alpha=0 for the unperturbed
Prob_function=matlabFunction(Prob); %Turns the probability function into a function handle
%% Finds most probable position
%Run as many times as you like after running the two sections above
clc
clf
%close all

iterations=1e3;     %Number of iterations in the Monte Carlo Algorithm
measurements=1e0;   %Number of times the algorithm is repeated
lambda=5e-1;        %The interval from the previous position where a new position can be randomized

x_start=zeros(N,1);
x_new=zeros(N,1);
x_saved=zeros(N,iterations);    %Final x-coordinates saved for the last iteration. Used in the animation
x_finals=zeros(N,measurements); %Final position of each measurement
P_finals=zeros(1,measurements); %Final probability of each measurement

for i=1:measurements
    %N random positions within the intervall [-square_radius,square_radius]
    for j=1:N
        square_radius=5;
        x_start(j)=-square_radius+(square_radius-(-square_radius)).*rand(1,1);
    end
    for k=1:iterations
        %N new random positions within an interval lambda of the previous
        %positions
        for l=1:N
            xi_1 = -1 + (1-(-1)).*rand(1,1);
            x_new(l)=x_start(l)+lambda*xi_1;
        end
        %Change the number of input arguments to be the same as N!
        
        %Calculates the probability of the new and old positions
        P_start=Prob_function(x_start(1),x_start(2));%,x_start(3),x_start(4));%,x_start(5));
        P_new=Prob_function(x_new(1),x_new(2));%,x_new(3),x_new(4));%,x_new(5));
        
        %If the new positions are more probable, the new positions will be set
        %as the starting postions for the next iteration
        %If not, a new postion is randomized within the interval lambda
        if P_new>P_start
            x_start=x_new;
        else
            x_start=x_start;
        end
        x_saved(:,k)=x_start;
    end

    
    %Plots the probability density. See density_operator_1D for more details
    Psi_up=sym('a', [n_up 1]);
    Psi_down=sym('a', [n_down 1]);
    m=1;
    w=1;
    for p=1:n_up
        Psi_up(p,1)=abs(QHO(m,w,p))^2;
    end
    for q=1:n_down
        Psi_down(q,1)=abs(QHO(m,w,q))^2;
    end
    x=linspace(-10,10,1000);
    Prob_up=sum(Psi_up);
    Prob_down=sum(Psi_down);
    P_up=matlabFunction(Prob_up);
    P_down=matlabFunction(Prob_down);   
    
    deltaY=1/(measurements+1);
    height=i*deltaY;
    if measurements==1
        height=0.6;
    end
    transparency=1;
    
    figure(1)
    %clf(1)
    hold on
    plot(x,P_down(x),'-.','Color',[0.8500, 0.3250, 0.0980],'Linewidth',2)
    plot(x,P_up(x),'--','Color',[0, 0.4470, 0.7410],'Linewidth',2)
    plot_up=plot(x_start(1:n_up),height*ones(n_up,1),'.','Color',[0, 0.4470, 0.7410],'Markersize',34);
    plot_down=plot(x_start(n_up+1:N),height*ones(n_down,1),'.','Color',[0.8500, 0.3250, 0.0980],'Markersize',30);
    axis([-5 5 0 0.8])
end
    
    %% Allows plots with multiple selected measurements above
    figure(2)
    height=0.4;
    hold on
    %Plots the most probable particle postions together with the probability density
    plot(x,P_down(x),'-.','Color',[0.8500, 0.3250, 0.0980],'Linewidth',2)
    plot(x,P_up(x),'--','Color',[0, 0.4470, 0.7410],'Linewidth',2)
    
    %Plots the particles
    scatter2 = scatter(x_start(n_up+1:N),height*ones(n_down,1),100,'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerEdgeColor',[0.8500, 0.3250, 0.0980]);
    scatter1 = scatter(x_start(1:n_up),height*ones(n_up,1),60,'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerEdgeColor',[0, 0.4470, 0.7410]); 
    %Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1.MarkerFaceAlpha = transparency;
    scatter1.MarkerEdgeAlpha = transparency;
    scatter2.MarkerFaceAlpha = transparency;
    scatter2.MarkerEdgeAlpha = transparency;

    axis([-5 5 0 0.8])
    set(gca,'FontSize',20)
    xl=xlabel('$\xi$','Interpreter','latex');
    xl.FontSize=35;
    
    texten=sprintf('\\alpha=%g',alpha);%,'Interpreter','latex');
    text(2.5,height+0.05,texten,'Fontsize',18)
    plot(x,height*ones(length(x),1),'k--')
    
    %grid minor
    %box on
    
    x_finals(:,i)=x_start;
    P_finals(i)=P_start;
%end
hold off

%% Animation of a single measurement
figure;
for i=1:2e2
    plot(x_saved(1:n_up,i),0.3*ones(n_up,1),'.','Color',[0, 0.4470, 0.7410],'Markersize',34)
    hold on
    plot(x_saved(n_up+1:N,i),0.3*ones(n_down,1),'.','Color',[0.8500, 0.3250, 0.0980],'Markersize',30)
    %plot(x,P_up(x),x,P_down(x),'Linewidth',2)
    axis([-5 5 0 1])
    pause(1e-2)
    hold off
end
