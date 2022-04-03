clc;

% Run only after Task1.m file
E = evalin('base', 'E');        % in Newton/Coulomb
B = evalin('base', 'B');        % in Tesla
c = 3*10^8;         % in m/s
m = 9.109*10^-31;   % in kgs
q = -1.602*10^-19;   % in Coulomb

% Initial Conditions
v = [0 0.9*c 0];    % initial 3D velocity
y_prev = zeros([1,3]);
f = @(v) q*(E + cross(v,B))/m;    

% Run only after Task1.m file
Tp = evalin('base', 'Tp');
t_step = Tp/4;
T = [10*Tp 100*Tp 1000*Tp 10000*Tp];
y = zeros([2+int32(T(end)/t_step), 6]);   % x, y, z, px, py, pz
y(1,4:6) = m*v;

for t=0:t_step:T(end-1)
    k1 = t_step*f(v);
    k2 = t_step*f(v+(k1/2));
    k3 = t_step*f(v+(k2/2));
    k4 = t_step*f(v+k3);
    v = v + (k1 + 2*k2 + 2*k3 + k4)/6;
    y(2+int32(t/t_step), 1:3) = y_prev + v*t_step;   % x, y, z position after each time step
    y_prev = y(2+int32(t/t_step), 1:3);   %Setting previous timestep positioning of electron
    y(2+int32(t/t_step), 4:6) = m*v;    % x, y, z momentum after each time step
%     plot3(y(1+int32(t/t_step), 1), y(1+int32(t/t_step), 2), y(1+int32(t/t_step), 3), '.r')
%     hold on;
%     pause(0.0002)
end

% Position and Momentum after corresponding 10, 100, 1000, 10000 turns
Pos = y(2+(T/t_step), 1:3);
Momen = y(2+(T/t_step), 4:6);

%% Analytical solution


V0 = [0 0.9*c 0]';

Y = zeros([length(T), 3]);
i = 1;

for t=T
    [Time,V] = ode45(@(time,v) q*(E' + cross(v,B'))/m, [0 t], V0);
    vel = V(end,:);   
    Y(i, :) = m*vel;
    i = i + 1
end

err_percent = 100*(Momen-Y)./Y;   % Percentage error between integrator method and analytical integration method  (for 0 values, it will show NaN)
