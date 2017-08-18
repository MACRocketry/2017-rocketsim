clear
clc

% Written By Curtis Graham for McMaster Rocketry Team
% August 15th 2017

%Calculation of oscillations caused by trag forces for various Cp & Cg

%Radial symetry assumed

%Variables:

theta=[];               %Angle from y axis CCW
omega=[];               %Angular velocity of theta
alpha=[];               %Angular acceleration of theta


I=[1/12*30*2^2];        %Moment of inertia

t=[500];                %time
delt=[0.01];            %timestep

Calib=[1];              %Seperation of Cp Cg

Fd=[300];               %Drag Forces

% Equations

%domega/dt=Calib*Fd*sin(theta)/I;
%dtheta/dt=omega;


% Runge-Kutta Solution
%Solves domega/dt for each timestep and uses a dummy to solve dtheta/dt

N=t/delt;


time=transpose([0:delt:t]); 

theta1=transpose(zeros(1,N+1));
theta2=transpose(zeros(1,N));
theta3=transpose(zeros(1,N));
theta4=transpose(zeros(1,N));

omega1=transpose(zeros(1,N+1));         
omega2=transpose(zeros(1,N));
omega3=transpose(zeros(1,N));
omega4=transpose(zeros(1,N));

domegadt1=transpose(zeros(1,N));        
domegadt2=transpose(zeros(1,N));
domegadt3=transpose(zeros(1,N));
domegadt4=transpose(zeros(1,N));

dthetadt1=transpose(zeros(1,N));        
dthetadt2=transpose(zeros(1,N));
dthetadt3=transpose(zeros(1,N));
dthetadt4=transpose(zeros(1,N));

Komega=transpose(zeros(1,N));  
Ktheta=transpose(zeros(1,N)); 
theta1(1)=0.1;

for i=1:N;
    
    domegadt1(i)=-Calib*Fd*sin(theta1(i))/I;
    omega2(i)=omega1(i)+domegadt1(i)*(delt/2);
    dthetadt2(i)=omega2(i);
    theta2(i)=theta1(i)+dthetadt1(i)*(delt/2);
    
    domegadt2(i)=-Calib*Fd*sin(theta2(i))/I;
    omega3(i)=omega1(i)+domegadt2(i)*(delt/2);
    dthetadt3(i)=omega3(i);
    theta3(i)=theta2(i)+dthetadt2(i)*(delt/2);
    
    domegadt3(i)=-Calib*Fd*sin(theta3(i))/I;
    omega4(i)=omega1(i)+domegadt3(i)*delt;
    dthetadt4(i)=omega4(i);
    theta4(i)=theta3(i)+dthetadt3(i)*delt;
    
    domegadt4(i)=-Calib*Fd*sin(theta4(i))/I;
    
    Komega(i)=1/6*(domegadt1(i)+2*domegadt2(i)+2*domegadt3(i)+domegadt4(i));
    Ktheta(i)=1/6*(dthetadt1(i)+2*dthetadt2(i)+2*dthetadt3(i)+dthetadt4(i));
    omega1(i+1)=omega1(i)+Komega(i)*delt;
    theta1(i+1)=theta1(i)+Ktheta(i)*delt;
    
end
plot(time, theta1*(180/pi))


