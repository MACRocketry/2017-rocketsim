function [V,A,F] = PlotTrogectory(Figure_Name,y,t,m,delta_t)
 
%Altitude over time
figure('Name',Figure_Name)
subplot(2,2,1)
plot(t,y)
title('Altitude')
xlabel('Time (s)')
ylabel('Altitude (m)')

%Velocity over time
V = diff(y)/delta_t;
V(length(t)) = 0;
subplot(2,2,3)
plot(t,V)
title('Velocity')
xlabel('Time (s)')
ylabel('Vellocity (m/s)')

%Acelleration over time
A = diff(V)/delta_t;
A(length(t)) = 0;
%A(length(t)-1) = 0;
subplot(2,2,2)
plot(t,A)
xlim([0,t(length(t))-2*delta_t])
title('Acelleration')
xlabel('Time (s)')
ylabel('Acelleration (m/s^2)')

%Force over time
F = m*A;
subplot(2,2,4)
plot(t,F)
xlim([0,t(length(t))-2*delta_t])
title('Force')
xlabel('Time (s)')
ylabel('Force (N)')
end