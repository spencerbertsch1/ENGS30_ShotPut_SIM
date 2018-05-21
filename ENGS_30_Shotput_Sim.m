%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGS 30 - Biological Physics
%~~~~~~~~~~~~~~~~~~~
% Shot Put Trajectory for different launch angles with quadratic drag 

% This script uses Euler stepping to simulate the launch of a shotput at
% different angles, showing the ideal launch angle of 38.5 degrees above
% zero. 

%
% Author: Spencer Bertsch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all variables
clear
clc
clf
%%%%%%%%%%%%%%%%%%%%
%% Variables

radius = 0.05; %(cm)
rho_ball = 7900; %(Kg/m^3)
m = (4/3)*pi*(radius^3) * rho_ball; 
A = pi*(radius^2); %area in cm^2
Cd = 0.47; %drag coefficient (unitless)
rho_air = 1.2; %(Kg/m^3)
g = 9.8; %(m/sec^2)

%Launch Angle
theta = 45; %launch angle above 0 rads
V0 = 10; %m/sec

Vx0=V0*cosd(theta);
Vy0=V0*sind(theta); 
x0 = 0;
y0 = 3;  


%% Numerical Solution
% Euler method is used as an integration scheme to calculate projectile motion

  %Time vector
    dt = 0.001; 
    t = 0:dt:5; 

Vx(1,1) = Vx0;
Vy(1,1) = Vy0;

V_total = ones(size(t));
V_total(1) = V0; 

x(1,1) =0; 
y(1,1) =0; 

%Compute drag constant outside the loop so we only have to do it once!
%Saves run time
Drag_Constants = 0.5*rho_air*Cd*A;

i=2; 
while y>=(0)
        
        %%%DRAG%%%
    % Find total drag on the ball at every dt 
    F_drag_total = ( Drag_Constants *( ( (Vx(i-1))^2 + (Vy(i-1))^2) ));
    
    V_total(i) = sqrt( (Vx(i-1))^2 + (Vy(i-1))^2);
    
    % Find x and y components of drag 
    A_Drag_x = (1/m) * -F_drag_total * ((Vx(i-1)) / V_total(i-1)); 
    A_Drag_y = (1/m) * -F_drag_total * ((Vy(i-1)) / V_total(i-1));
    
    % Calculate x component of acceleration, velocity, and position
    gravity_x = 0;  %acceleration due to gravity is zero in x direction
    
    Ax = gravity_x + A_Drag_x;
    Vx(i) = Vx(i-1) + dt * (Ax);
    x(i) = x(i-1) + dt * (Vx(i));
    
    % Calculate y component of acceleration, velocity, and position
    gravity_y = (-g); %negative acceleration due to gravity
    
    Ay = gravity_y + A_Drag_y;   
    Vy(i) = Vy(i-1) + dt * (Ay);
    y(i) = y(i-1) + dt * (Vy(i-1));
    
    i=i+1; 
end




%% Plots 
% comet(x,y) %<--- UNCOMMENT to see a quick movie of the cannon ball's path
% plot(x,y,'b')
% xlabel('Horizontal Distance (m)') % Label axes with units
% ylabel('Height (m)')
% title('Cannonball Trajectory')
% legend('Numerical (with drag)','Analytical (without drag)')
% grid on 


                          %%%% Dynamic Plot %%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',[97   22   1216   776]);
for pn = 1:10:(length(x))
    %figure(1)
    clf
    plot(x(pn),y(pn),'r.','markersize',70); %size of ball is proportional to its mass
    hold on 
    plot(x(1:pn),y(1:pn),'b--','linewidth',2)
    grid minor
    %axis square
    %h=gca; 
    axis([0 (10) 0 (2.2)]);
    title('Trajectory of Shot','fontSize',30)
    pause(0.003)

end

xlabel('Distance (meters)','fontSize',25)
ylabel('Height (Meters)','fontSize',25)
h = legend('Position of Shot');
set(h,'FontSize',15);

