%Pauli crystals
%% Most probable position in 1D
%This code calculates the most probable positions for 1-10 fermionic 
%particles in a 1D harmonic oscillator potential using the Monte Carlo 
%algorithm. 

%First it constructs an analytical expression for the antisymmetric 
%many-body wave function using the Slater determinant and symbolic 
%variables. This is done in N_particles_1D.m. 
%It then implements the Monte Carlo algorithm, which after randomising the
%positions of N particles, enters a loop where the particles are only
%allowed to move within a square 2 lambda of their previous position of the
%probability increases.
%Eventually the particles reach the most probable configuration.

clc
close all
N=3;            %Number of particles 1-10

Psi=N_particles_1D(N);      %Constructs the many body wave function
Prob=abs(Psi)^2;            %Creates the probability distribution
Prob_function=matlabFunction(Prob);     %Turns the probability function into a function handle

iterations=1e4;     %Number of iterations in the Monte Carlo Algorithm
lambda=1e-2;        %The interval from the previous position where a new postion can be randomized


x_start=zeros(N,1);
x_new=zeros(N,1);

%N random positions within the intervall [-square_radius,square_radius]
for i=1:N
    square_radius=10;
    x_start(i)=-square_radius+(square_radius-(-square_radius)).*rand(1,1);
end

for j=1:iterations
    %N new random positions within an interval lambda of the previous
    %positions
    for k=1:N
        xi_1 = -1 + (1-(-1)).*rand(1,1);
        x_new(k)=x_start(k)+lambda*xi_1;
    end
    %Change the number of input arguments to be the same as N
    %Calculates the probability of the new and old positions
    P_start=Prob_function(x_start(1),x_start(2),x_start(3));
    P_new=Prob_function(x_new(1),x_new(2),x_new(3));
    
    %If the new positions are more probable, the new positions will be set
    %as the starting postions for the next iteration
    %If not, a new postion is randomized within the interval lambda
    if P_new>P_start
        x_start=x_new;
    else
        x_start=x_start;
    end
end

%Plots the probability density. See density_operator_1D for more details
Psi_up=sym('a', [N 1]);
m=1;
w=1;
for i=1:N
        Psi_up(i,1)=abs(QHO(m,w,i))^2;
end
x=linspace(-10,10,1000);
Prob_up=sum(Psi_up);
P_up=matlabFunction(Prob_up);

%Plots the most probable particle postions together with the probability
%density
figure;
plot(x,P_up(x),'Linewidth',2)
hold on
plot(x_start(:),0.5*ones(N,1),'.','Markersize',30)
axis([-5 5 0 1])
set(gca,'FontSize',20)
xl=xlabel('$\xi$','Interpreter','latex');
xl.FontSize=35;
grid minor
hold off
%% Most probable position 2D
%This code calculates the most probable positions for 1,3 and 6
%fermionic particles in a 2D harmonic oscillator potential using the Monte 
%Carlo algorithm.

%First it constructs an analytical expression for the antisymmetric 
%many-body wave function using the Slater determinant and symbolic 
%variables. This is done in N_particles_2D.m. 
%It then implements the Monte Carlo algorithm, which after randomising the
%positions of N particles, enters a loop where the particles are only
%allowed to move within a square 2 lambda of their previous position of the
%probability increases.
%Eventually the particles reach the most probable configuration.

tic
clc
close all
N=6;       %Number of particles 1,3 and 6 (10 particles takes too much memory)

Psi=N_particles_2D(N);      %Constructs the many body wave function
Prob=abs(Psi)^2;            %Creates the probability distribution
Prob_function=matlabFunction(Prob);     %Turns the probability function into a function handle

iterations=1e4;     %Number of iterations in the Monte Carlo algorithm
measurements=1e0;   %Allows multiple measurements as to study the rotation symmetry
lambda=1e-2;        %Forms a square with the sides 2*lambda where new positions can be randomized 

%Empty vectors for the positions
x_start=zeros(N,1);
y_start=zeros(N,1);
pos_new=zeros(N,2);
all_positions=zeros(N,2,measurements);
for l=1:measurements
    %Randomize starting postions within the square
    %[-square_radius,square_radius]x[-square_radius,square_radius]
    for i=1:N
        square_radius=5;
        x_start(i)=-square_radius+(square_radius-(-square_radius)).*rand(1,1);
        y_start(i)=-square_radius+(square_radius-(-square_radius)).*rand(1,1);
    end
    pos_start=[x_start,y_start];
    %N new random positions within an interval lambda of the previous
    %positions
    for j=1:iterations
        for k=1:N
            xi_1 = -1 + (1-(-1)).*rand(1,1);
            xi_2 = -1 + (1-(-1)).*rand(1,1);
            pos_new(k,:)=pos_start(k,:)+lambda*[xi_1,xi_2];
        end
        %Calculates the probabilities for the new and old positions depending
        %on particle number
        if N==1
            P_start=Prob_function(pos_start(1,1),pos_start(1,2));
            P_new=Prob_function(pos_new(1,1),pos_new(1,2));
        elseif N==3
            P_start=Prob_function(pos_start(1,1),pos_start(1,2),pos_start(2,1),pos_start(2,2),pos_start(3,1),pos_start(3,2));
            P_new=Prob_function(pos_new(1,1),pos_new(1,2),pos_new(2,1),pos_new(2,2),pos_new(3,1),pos_new(3,2));
        elseif N==6
            P_start=Prob_function(pos_start(1,1),pos_start(1,2),pos_start(2,1),pos_start(2,2),pos_start(3,1),pos_start(3,2),pos_start(4,1),pos_start(4,2),pos_start(5,1),pos_start(5,2),pos_start(6,1),pos_start(6,2));
            P_new=Prob_function(pos_new(1,1),pos_new(1,2),pos_new(2,1),pos_new(2,2),pos_new(3,1),pos_new(3,2),pos_new(4,1),pos_new(4,2),pos_new(5,1),pos_new(5,2),pos_new(6,1),pos_new(6,2));
        elseif N==10
            P_start=Prob_function(pos_start(1,1),pos_start(1,2),pos_start(2,1),pos_start(2,2),pos_start(3,1),pos_start(3,2),pos_start(4,1),pos_start(4,2),pos_start(5,1),pos_start(5,2),pos_start(6,1),pos_start(6,2),pos_start(7,1),pos_start(7,2),pos_start(8,1),pos_start(8,2),pos_start(9,1),pos_start(9,2),pos_start(10,1),pos_start(10,2));
            P_new=Prob_function(pos_new(1,1),pos_new(1,2),pos_new(2,1),pos_new(2,2),pos_new(3,1),pos_new(3,2),pos_new(4,1),pos_new(4,2),pos_new(5,1),pos_new(5,2),pos_new(6,1),pos_new(6,2),pos_new(7,1),pos_new(7,2),pos_new(8,1),pos_new(8,2),pos_new(9,1),pos_new(9,2),pos_new(10,1),pos_new(10,2));
        end
        %If the new positions are more probable, the new positions will be set
        %as the starting postions for the next iteration
        %If not, a new postion is randomized within the interval lambda
        if P_new>P_start
            pos_start=pos_new;
        else
            pos_start=pos_start;
        end
    end
    %Plots the most probable particle positions
    hold on
    plot(pos_start(:,1),pos_start(:,2),'.b','Markersize',25)
    
    all_positions(:,:,l)=pos_start;
end
axis equal
set(gca,'XLim',[-2.5 2.5],'XTick',-2:1:2)
set(gca,'YLim',[-2.5 2.5],'YTick',-2:1:2)
set(gca,'FontSize',18)
xl=xlabel('$\xi$','Interpreter','latex');
x2=ylabel('$\eta$','Interpreter','latex');
xl.FontSize=25;
x2.FontSize=25;
texten=sprintf('N=%i',N);
text(1,2,texten,'Fontsize',24)
box on

hold off
toc


%% Loading and saving files
%The cell Recreating Pauli crystals in 1D/2D takes a while to run. Because of
%this it is recommended to instead load a file that was previously done.
%Pick one "load" and proceed to the cells called Center of Mass in the 
%appropriate dimension after running this cell. 

%If you have done a long calculation and want to save it copy the line
%below into the command window and give it a new file name.
%save('filename','all_positions','-mat')

%The first number is for the number of particles and the second for the
%number of measurements. Pick one.
%structure_array=load('Figurer/1D/Actual_crystals/3_e4_1D','-mat');
%structure_array=load('Figurer/2D/Actual_Pauli_crystals/3_e4','-mat');
%structure_array=load('Figurer/2D/Actual_Pauli_crystals/3_e4_circular','-mat');
%structure_array=load('Figurer/2D/Actual_Pauli_crystals/6_e4_circular','-mat');
structure_array=load('Figurer/2D/Actual_Pauli_crystals/6_5e4_circular','-mat');
%structure_array=load('Figurer/2D/Actual_Pauli_crystals/6_e3','-mat');
%structure_array=load('Figurer/2D/Actual_Pauli_crystals/6_e4','-mat');

all_positions=structure_array.all_positions;
measurements=length(all_positions(1,1,:));
N=length(all_positions(:,1,1));

%% Recreating Pauli crystals in 1D
%This code is spit up into several parts and plots.

%This first part, mostly uses the same algorithm as the cell
%above, with one main difference. The particles do not always move to the
%most probable position. Instead they sometimes move to a less probable
%case. This leads to more spread out measurements, similar to what one
%would observe in the laboratory. The algorithm is called the Metropolis
%algorithm.

%After the spread out measurements have been created, each measurements is
%moved so that the particles centre of mass is located at the origin. 
%Hopefully this reveals the probability density. 

%After each step there will be the possibility to view the plot as a
%histogram. This will reveal much more than the ordinary plots since the
%particles often overlap.


tic
clc
close all
N=3;            %Number of particles 1-10

Psi=N_particles_1D(N);      %Constructs the many body wave function
Prob=abs(Psi)^2;            %Creates the probability distribution
Prob_function=matlabFunction(Prob);     %Turns the probability function into a function handle

iterations=1e4;     %Number of iterations in the Monte Carlo Algorithm
lambda=1e-2;        %The interval from the previous position where a new postion can be randomized

measurements=1e0;   %Needs to be >1e4 for clear histograms

x_start=zeros(N,1);
x_new=zeros(N,1);
all_positions=zeros(N,1,measurements);
hist=[];

for l=1:measurements
    %N random positions within the intervall [-square_radius,square_radius]
    for i=1:N
        square_radius=10;
        x_start(i)=-square_radius+(square_radius-(-square_radius)).*rand(1,1);
    end
    
    for j=1:iterations
        %N new random positions within an interval lambda of the previous
        %positions
        for k=1:N
            xi_1 = -1 + (1-(-1)).*rand(1,1);
            x_new(k)=x_start(k)+lambda*xi_1;
        end
        %Change the number of input arguments to be the same as N
        %Calculates the probability of the new and old positions
        P_start=Prob_function(x_start(1),x_start(2),x_start(3));%,x_start(4));
        P_new=Prob_function(x_new(1),x_new(2),x_new(3));%,x_new(4));
        
        %If the new positions are more probable, the new positions will be set
        %as the starting postions for the next iteration
        %If not, there are two cases. If a random variable xi_3 between 0
        %and one is smaller than p, the new position is set as the starting
        %position for the next iteration. Otherwise, the same starting
        %position is used again.
        p=P_new/P_start;
        if p>1
            x_start=x_new;
        else
            xi_3=rand(1,1);
            if xi_3<=p
                x_start=x_new;
            else
                x_start=x_start;
            end
        end
    end
    hold on
    plot(x_start(:),0.5*ones(N,1),'.b','Markersize',25)
    
    all_positions(:,:,l)=x_start;
    hist=[hist;x_start];
end
%Makes plot pretty
figure(1)
axis([-5 5 0 1])
set(gca,'FontSize',20)
xl=xlabel('$\xi$','Interpreter','latex');
xl.FontSize=35;
grid minor
box on
hold off
toc
%% Histogram without image processing
%Needed to be able to extrapolate anything.
figure(2)
clf(2)

side_length=5;
box_length=0.1;

histogram(hist,(-side_length:box_length:side_length))

%% Center of mass
clc
all_pos_cm=zeros(N,1,measurements);
hist_cm=[];
%Finds the center of mass and moves it to the origin.
for i=1:measurements
    x_cm=1/N*sum(all_positions(:,:,i));
    all_pos_cm(:,:,i)=all_positions(:,:,i)-x_cm;
    hist_cm=[hist_cm;all_pos_cm(:,:,i)];
end
figure(3)
clf(3)
%Plots all measurements
for j=1:measurements
    hold on
    plot(all_pos_cm(:,1,j),0.5*ones(N,1),'.b','Markersize',25)
end
axis([-10 10 0 1])
set(gca,'FontSize',20)
xl=xlabel('$\xi$','Interpreter','latex');
xl.FontSize=35;
grid minor
box on
hold off

%% Histogram with center of mass correction
%View the plot generated above as a histogram.
clc
figure(4);
clf(4)

side_length=5;
box_length=0.1;

histogram(hist_cm,(-side_length:box_length:side_length))

%% Recreating Pauli crystals in 2D
%This code is spit up into several parts and plots and does the same as
%above but in 2D and with one addition.

%Because of rotational symmentry the particles are rotated such that each 
%measurement aligns as closely as possible with the expected Pauli crystlas. 

%TO REDUCE COMPUTING TIME: Remove (comment) the plots in the cells that
%actually compute stuff. As they don't plot histogram there won't be any
%information in them when you exceed like 5 particles.

tic
clc
close all
N=6;       %Number of particles: 3 and 6

Psi=N_particles_2D(N);      %Constructs the many body wave function
Prob=abs(Psi)^2;            %Creates the probability distribution
Prob_function=matlabFunction(Prob);     %Turns the probability function into a function handle

iterations=1e4;     %Number of iterations in the Metropolis algorithm
lambda=1e-2;        %Forms a square with the sides 2*lambda where new positions can be randomized 

measurements=5e4;   %The more measurements, the better histograms. 1e3 is good for a quick look. 1e4 takes a while but gives good looking histograms.


%Empty vectors for the positions
x_start=zeros(N,1);
y_start=zeros(N,1);
pos_new=zeros(N,2);
all_positions=zeros(N,2,measurements);
hist=[];

figure(1);
for l=1:measurements
    %Randomize starting postions within the square
    %[-square_radius,square_radius]x[-square_radius,square_radius] or 
    %a circle of radius circle_radius
    for i=1:N
%         square_radius=5;
%         x_start(i)=-square_radius+(square_radius-(-square_radius)).*rand(1,1);
%         y_start(i)=-square_radius+(square_radius-(-square_radius)).*rand(1,1);
        circle_radius=5;
        circle_angle=2*pi;
        r_start=circle_radius*rand(1,1);
        theta_start=circle_angle*rand(1,1);
        [x_start(i),y_start(i)]=pol2cart(theta_start,r_start);
    end
    pos_start=[x_start,y_start];
    %N new random positions within an interval lambda of the previous
    %positions
    for j=1:iterations
        for k=1:N
            xi_1 = -1 + (1-(-1)).*rand(1,1);
            xi_2 = -1 + (1-(-1)).*rand(1,1);
            pos_new(k,:)=pos_start(k,:)+lambda*[xi_1,xi_2];
        end
        %Calculates the probabilities for the new and old positions depending
        %on particle number
        if N==1
            P_start=Prob_function(pos_start(1,1),pos_start(1,2));
            P_new=Prob_function(pos_new(1,1),pos_new(1,2));
        elseif N==3
            P_start=Prob_function(pos_start(1,1),pos_start(1,2),pos_start(2,1),pos_start(2,2),pos_start(3,1),pos_start(3,2));
            P_new=Prob_function(pos_new(1,1),pos_new(1,2),pos_new(2,1),pos_new(2,2),pos_new(3,1),pos_new(3,2));
        elseif N==6
            P_start=Prob_function(pos_start(1,1),pos_start(1,2),pos_start(2,1),pos_start(2,2),pos_start(3,1),pos_start(3,2),pos_start(4,1),pos_start(4,2),pos_start(5,1),pos_start(5,2),pos_start(6,1),pos_start(6,2));
            P_new=Prob_function(pos_new(1,1),pos_new(1,2),pos_new(2,1),pos_new(2,2),pos_new(3,1),pos_new(3,2),pos_new(4,1),pos_new(4,2),pos_new(5,1),pos_new(5,2),pos_new(6,1),pos_new(6,2));
        elseif N==10
            P_start=Prob_function(pos_start(1,1),pos_start(1,2),pos_start(2,1),pos_start(2,2),pos_start(3,1),pos_start(3,2),pos_start(4,1),pos_start(4,2),pos_start(5,1),pos_start(5,2),pos_start(6,1),pos_start(6,2),pos_start(7,1),pos_start(7,2),pos_start(8,1),pos_start(8,2),pos_start(9,1),pos_start(9,2),pos_start(10,1),pos_start(10,2));
            P_new=Prob_function(pos_new(1,1),pos_new(1,2),pos_new(2,1),pos_new(2,2),pos_new(3,1),pos_new(3,2),pos_new(4,1),pos_new(4,2),pos_new(5,1),pos_new(5,2),pos_new(6,1),pos_new(6,2),pos_new(7,1),pos_new(7,2),pos_new(8,1),pos_new(8,2),pos_new(9,1),pos_new(9,2),pos_new(10,1),pos_new(10,2));
        end
        %If the new positions are more probable, the new positions will be set
        %as the starting postions for the next iteration
        %If not, there are two cases. If a random variable xi_3 between 0
        %and one is smaller than p, the new position is set as the starting
        %position for the next iteration. Otherwise, the same starting
        %position is used again.
        p=P_new/P_start;
        if p>1
            pos_start=pos_new;
        else
            xi_3=rand(1,1);
            if xi_3<=p
                pos_start=pos_new;
            else
                pos_start=pos_start;
            end
        end
    end
    %Plots the particle positions for each measurement
    hold on
    plot(pos_start(:,1),pos_start(:,2),'.b','Markersize',25)
    
    all_positions(:,:,l)=pos_start;
    %Saves values for the histogram
    hist=[hist;pos_start];
end
%Makes the plot pretty
axis equal
set(gca,'XLim',[-5.5 5.5],'XTick',-5:2:5)
set(gca,'YLim',[-5.5 5.5],'YTick',-5:2:5)
set(gca,'FontSize',18)
xl=xlabel('$\xi$','Interpreter','latex');
x2=ylabel('$\eta$','Interpreter','latex');
xl.FontSize=25;
x2.FontSize=25;
texten=sprintf('N=%i',N);
text(3,4,texten,'Fontsize',24)
box on

hold off
toc

%% Histogram without image processing
%View the plot generated above as a histogram.
clc
figure(2);
clf(2)

side_length=5;
box_length=0.1;

hist3(hist,'Ctrs',{-side_length:box_length:side_length -side_length:box_length:side_length},'CDataMode','auto','FaceColor','interp','EdgeColor','interp')
axis equal
xl=xlabel('$\xi$','Interpreter','latex');
x2=ylabel('$\eta$','Interpreter','latex');
xl.FontSize=20;
x2.FontSize=20;
colorbar
view(2)

%% Center of mass
%Aligns all the measurements so that the center of mass of each measurments
%ends up at the origin.

clc
all_pos_cm=zeros(N,2,measurements);
%Finds the center of mass and moves it to the origin.
for i=1:measurements
    x_cm=1/N*sum(all_positions(:,:,i));
    all_pos_cm(:,:,i)=all_positions(:,:,i)-x_cm;
end
figure(3)
clf(3)
%Plots all measurements
for j=1:measurements
    hold on
    plot(all_pos_cm(:,1,j),all_pos_cm(:,2,j),'.b','Markersize',25)
end
%Makes plot pretty
axis equal
set(gca,'XLim',[-5.5 5.5],'XTick',-5:2:5)
set(gca,'YLim',[-5.5 5.5],'YTick',-5:2:5)
set(gca,'FontSize',18)
xl=xlabel('$\xi$','Interpreter','latex');
x2=ylabel('$\eta$','Interpreter','latex');
xl.FontSize=25;
x2.FontSize=25;
texten=sprintf('N=%i',N);
text(3,4,texten,'Fontsize',24)
box on

%Saves values for the histogram
hist_cm=[];
for k=1:measurements
    hist_cm=[hist_cm; all_pos_cm(:,:,k)];
end
hold off
%% Histogram with center of mass correction
%View the plot generated above as a histogram.
clc
figure(4);
clf(4)

side_length=5;
box_length=0.1;

hist3(hist_cm,'Ctrs',{-side_length:box_length:side_length -side_length:box_length:side_length},'CDataMode','auto','FaceColor','interp','EdgeColor','interp')
axis equal
xl=xlabel('$\xi$','Interpreter','latex');
x2=ylabel('$\eta$','Interpreter','latex');
xticks([-5 -3 -1 1 3 5])
xticklabels({'-5','-3','-1','1','3','5'})
yticks([-5 -3 -1 1 3 5])
yticklabels({'-5','-3','-1','1','3','5'})
ax = gca;
ax.FontSize = 16; 
xl.FontSize=20;
x2.FontSize=20;
colorbar
%set(gca,'YTickLabel',[]);
%set(gca,'XTickLabel',[]);
view(2)

%% Rotation
%Rotates each measurements so that the particles align as closely as
%possible with the predicted Pauli crystal.

clc

if N==3
    k=3;
elseif N==6
    k=5;
end
%Set angle_max to 2*pi/k for three particles. Still trying to figure out
%the best option for six particles.
angle_max=2*pi/k;

%Rotation matrix
R=@(angle) [cos(angle),-sin(angle);sin(angle),cos(angle)];

%Matrixes for the optimal rotation angle and particle positions for each
%measurement
vinkel=zeros(measurements,1);
all_pos_rot=zeros(N,2,measurements);

for k=1:measurements
    %Pick one measurment
    pos_cm=all_pos_cm(:,:,k);
    pos_rot=zeros(N,2);
    
    %Sets the not-rotated values as starting values
    pos_rot_opt=pos_cm;
    d=pauli_distance(N,pos_cm);
    
    %Steps through angles between zero and angle_max to find the angle that
    %gives the smallest d. d is defined in pauli_distance.m.
    for i=0:1e-2:angle_max
        %Rotates each particle in a measurement.
        for j=1:N
            pos_rot(j,:)=(R(i)*(pos_cm(j,:))')';
        end
        %Calculated the d for the rotation.
        d_new=pauli_distance(N,pos_rot);
        
        %If d is smaller, the position is set as the new optimal position
        if d_new<d
            d=d_new;
            pos_rot_opt=pos_rot;
            vinkel_opt=i;
        end
    end
    %Saves the optimal angles and positions for each measurement.
    vinkel(k)=vinkel_opt;
    all_pos_rot(:,:,k)=pos_rot_opt;
end
%Plots everything
figure(5)
clf(5)
for i=1:measurements
    hold on
    plot(all_pos_rot(:,1,i),all_pos_rot(:,2,i),'.b','MarkerSize',25);
end
%Makes the plot pretty
axis equal
set(gca,'XLim',[-5.5 5.5],'XTick',-5:2:5)
set(gca,'YLim',[-5.5 5.5],'YTick',-5:2:5)
set(gca,'FontSize',18)
xl=xlabel('$\xi$','Interpreter','latex');
x2=ylabel('$\eta$','Interpreter','latex');
xl.FontSize=25;
x2.FontSize=25;
texten=sprintf('N=%i',N);
text(3,4,texten,'Fontsize',24)
box on
hold off
hist_rot=[];
for k=1:measurements
    hist_rot=[hist_rot; all_pos_rot(:,:,k)];
end
%% Histogram with rotation
%View the plot generated above as a histogram.
clc
figure(6);
clf(6)

side_length=5;
box_length=0.1;

hist3(hist_rot,'Ctrs',{-side_length:box_length:side_length -side_length:box_length:side_length},'CDataMode','auto','FaceColor','interp','EdgeColor','interp')
axis equal
xl=xlabel('$\xi$','Interpreter','latex');
x2=ylabel('$\eta$','Interpreter','latex');
xticks([-5 -3 -1 1 3 5])
xticklabels({'-5','-3','-1','1','3','5'})
yticks([-5 -3 -1 1 3 5])
yticklabels({'-5','-3','-1','1','3','5'})
ax = gca;
ax.FontSize = 16; 
xl.FontSize=20;
x2.FontSize=20;
colorbar
%set(gca,'YTickLabel',[]);
%set(gca,'XTickLabel',[]);
view(2)

%box on
%set(gca,'linewidth',5)
