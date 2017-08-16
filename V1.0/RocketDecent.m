function [x,y,t] = RocketDecent(Mass,Gravity,D_Rocket,D_Droge, D_Main,delta_t,x0,y0,v0,theta0,deployment_hight_D, deployment_hight_M)

%%rocket Perameters
m = Mass;%mass in Kilograms
g = Gravity; %gravitational constant

%%initial conditions of Rocket motion
x(1)= x0;
y(1)= y0;
v = v0; %m/s
theta = theta0; %deg
vx = v*cosd(theta);
vy = v*sind(theta);
t(1) = 0;

%shoot opinings
%Timed based release
Tau1 = 15;%shoot 1
Tau2 = 200;%shoot 2
 
%Start Loop
i = 1;
while min(y)> -0.001
    if (y(i) < deployment_hight_M)%could make an or stament for the two as redundantcy
    %if (t(i) > Tau2)
        ax = -(D_Rocket/m)*v*vx;
        ay = -g - (D_Rocket/m)*v*vy - (D_Droge/m*v*vy) - (D_Main/m*v*vy);
        vx = vx + ax*delta_t;
        vy = vy + ay*delta_t;
        v = sqrt(vx^2 + vy^2);
    elseif (y(i) < deployment_hight_D)
    %elseif (t(i) > Tau1)
        ax = -(D_Rocket/m)*v*vx;
        ay = -g - (D_Rocket/m)*v*vy - (D_Droge/m*v*vy);
        vx = vx + ax*delta_t;
        vy = vy + ay*delta_t;
        v = sqrt(vx^2 + vy^2);
    else
        ax = -(D_Rocket/m)*v*vx;
        ay = -g - (D_Rocket/m)*v*vy;
        vx = vx + ax*delta_t;
        vy = vy + ay*delta_t;
        v = sqrt(vx^2 + vy^2);
    end
    
    x(i+1) = x(i) + vx*delta_t ;%+ 0.5*ax*delta_t^2;
    y(i+1) = y(i) + vy*delta_t ;%+ 0.5*ay*delta_t^2;
    t(i+1) = t(i) + delta_t;        % the above will make it more acurate but increase computation time
    i = i + 1;
end

end