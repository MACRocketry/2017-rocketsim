%Rocket Main


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This file will contain functions for all stages of Rocket Flight

%So far this contains the following:
    %basic Rocket Launch (treated as point particel
    %Rocket Decent with Parachoots


%future Additions:
    %more realistion model of rocket launch with center of pressure and
    %center of gravity
    
    %better model of Parachoot deployment
    
    %may turn into a GUI if apropreate
    %this may help with testing many parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Rokcet Parameters
m = 68/2.205;%mass in Kilograms
g = 9.81;  %m/s^2
r = 0.15/2; %15cm diameter
A = pi*r^2; %m^2
C = 0.75; %Drag Coefficient of a sphere
rho = 1.2; %kg/m^3 (density of air)
D_Rocket = rho*C*A/2;
%%Shoot Parameters
Cd = 1.75; %https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/rktvrecv.html
rho = 1.2; %assuming same every where.
A_Droge = 0.5; %m^2
A_Main = 4; %m^2
D_Droge = Cd*rho*A_Droge/2;
D_Main = Cd*rho*A_Main/2;

%%%%%%%Rocket Launch%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial conditions
delta_t = 0.01; %s
delta_m = 0.405; %kg/s
x0 = 0; %m
z0 = 0; %m
y0 = 0; %m
v0 = 0; %m/s
theta0 = 90; %deg
phi0 = 90; %deg
%basic Rocket launch function 
[xl,yl,tl,vl,theta,phi,i] = RocketLaunch(m,g,D_Rocket,delta_t,delta_m,x0,y0,z0,v0,theta0,phi0);
PlotTrogectory('Rocket Launch',yl,tl,m,delta_t);

%%%%%%%%Rocket Decent%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial conditions
delta_t = 0.01; %s
delta_m = 0.405; %kg/s
x0 = xl(i); %m
y0 = yl(i); %m
v0 = vl; %m/s
theta0 = theta; %deg
deployment_hight_D = 2800;
deployment_hight_M = 450;

%Rocket Decent with parachoots
[xd,yd,td] = RocketDecent(m,g,D_Rocket,D_Droge, D_Main,delta_t,x0,y0,v0,theta0,deployment_hight_D, deployment_hight_M);
[V,A,F] = PlotTrogectory('Rocket Decent',yd,td,m,delta_t);


%%%%%%Check to see if the Parachoots and eye bolts withstood deplayment%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



