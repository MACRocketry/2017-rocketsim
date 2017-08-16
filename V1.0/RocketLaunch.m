function [x,y,t,va,theta_f,phi_f,i] = RocketLaunch(Mass,Gravity,D_Rocket,delta_t,delta_m,x0,y0,z0,v0,theta0, phi0)
m = Mass;
g = Gravity;
%initial conditions of Rocket motion
x(1) = x0;
z(1) = z0;
y(1) = y0;
v = v0; %m/s
theta = theta0; %deg
phi = phi0; %deg
vx(1) = v*cosd(theta)*sind(phi);
vz(1) = v*cosd(theta)*cosd(phi);
vy(1) = v*sind(theta);
t(1) = 0;
ax(1) = 0;
ay(1) = 0;

%Start Loop
i = 1;
while (min(y)> -0.001)
    if (max(y)==y(i))
        
        %Calculate Force
        %Thrust at time steps,     also assume constant burn rate
        if (t(i) <= 3.5)
            Thrust = 2000;
            m = m - delta_m*delta_t;
        elseif ((t(i)>3.5)&&(t(i)<5))
            Thrust =  - 333.33*(t(i) - 9.5);
            m = m - delta_m*delta_t;
        elseif((t(i)>5)&&(t(i)<6))
            Thrust =  - 1500*(t(i) - 6);
            m = m - delta_m*delta_t;
        else
            Thrust = 0;
        end
        %Thrust
        Thrustx(i+1) = Thrust*cosd(theta(i))*sind(phi(i));
        Thrustz(i+1) = Thrust*cosd(theta(i))*cosd(phi(i));
        Thrusty(i+1) = Thrust*sind(theta(i));
        %Drag
        Drag(i+1) = -D_Rocket*(vx(i)^2+vy(i)^2+vz(i)^2);
        %Gravity
        G = -g*m;
        %Total Force
        Fx(i+1) = Drag(i+1)*cosd(theta(i))*sind(phi(i)) + Thrustx(i+1);
        Fz(i+1) = Drag(i+1)*cosd(theta(i))*cosd(phi(i)) + Thrustz(i+1);
        Fy(i+1) = G + Drag(i+1)*sind(theta(i)) + Thrusty(i+1);
        
        %Calculate Acceleration
        ax(i+1) = Fx(i+1)/m;
        az(i+1) = Fz(i+1)/m;
        ay(i+1) = Fy(i+1)/m;
        
        %Calculate Velocity
        vx(i+1) = vx(i) + ax(i+1)*delta_t;
        vz(i+1) = vz(i) + az(i+1)*delta_t;
        vy(i+1) = vy(i) + ay(i+1)*delta_t;
        %v = sqrt(vx^2 + vy^2);
        
        %Calculate Position
        x(i+1) = x(i) + vx(i+1)*delta_t ;%+ 0.5*ax(i-1)*delta_t^2;
        z(i+1) = z(i) + vz(i+1)*delta_t ;%+ 0.5*az(i-1)*delta_t^2;
        y(i+1) = y(i) + vy(i+1)*delta_t ;%+ 0.5*ay(i-1)*delta_t^2;
        t(i+1) = t(i) + delta_t;
        
        %Calculate angel of inclination
        phi(i+1) = phi(i);  %Phi will not change in this basic system
        theta(i+1) = atand(vy(i+1)/(vx(i+1)^2 + vz(i+1)^2)); %sqrt(vx(i+1)^2+vx(i+1)^2));
        i = i + 1;
    else
        break
    end
    
end
va  = sqrt(vx(i)^2+vz(i)^2+vy(i)^2);
theta_f = theta(i);
phi_f = phi(i);
% x(i)
% y(i)
% vx(i)
% vy(i)

end

