%% Runge Kutta 6th Order

clc;clear

tStep = 10^-12;
T = 5*10^-10;   % 5*10^-11 with tStep=10^-12 gives 2 revolutions of electron
E = [0 0 0];        % in Newton/Coulomb
B = [0 0 1];        % in Tesla
c = 3*10^8;         % in m/s
m = 9.109*10^-31;   % in kgs
q = -1.602*10^-19;   % in Coulomb
i = 0;
% Initial Conditions
v = [0 0.9*c 0];    % initial 3D velocity
y = zeros([2+int32(T/tStep), 6]);   % x, y, z, px, py, pz
y_prev = zeros([1,3]);
y(1,4:6) = m*v;
f = @(v) q*(E + cross(v,B))/m;

for t=0:tStep:T
    % Runge Kutta 6th Order Numerical Integration
    k1 = tStep*f(v);
    k2 = tStep*f(v+k1);
    k3 = tStep*f(v+(3*k1+k2)/8);
    k4 = tStep*f(v+(8*k1+2*k2+8*k3)/27);
    k5 = tStep*f(v+((-21+9*21^(1/2))*k1-8*(7-21^(1/2))*k2+48*(7-21^(1/2))*k3-3*(21-21^(1/2))*k4)/392);
    k6 = tStep*f(v+(-5*(231+51*21^0.5)*k1-40*(7+21^0.5)*k2-320*21^0.5*k3+3*(21+121*21^0.5)*k4+392*(6+21^0.5)*k5)/1960);
    k7 = tStep*f(v+(15*(22+7*21^0.5)*k1+120*k2+40*(-5+7*21^0.5)*k3-63*(-2+3*21^0.5)*k4-14*(49+9*21^0.5)*k5+70*(7-21^0.5)*k6)/180);
    v = v + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180;
    y(2+int32(t/tStep), 1:3) = y_prev + v*tStep;   % x, y, z position after each time step
    y_prev = y(2+int32(t/tStep), 1:3);   %Setting previous timestep positioning of electron
    y(2+int32(t/tStep), 4:6) = m*v;    % x, y, z momentum after each time step
    if t~=0 && ~isequal(sign(y(2+int32(t/tStep)-1, 2)), sign(y(2+int32(t/tStep), 2)))
        i = i+1
        Time(i) = t;
    end
    plot3(y(2+int32(t/tStep),1), y(2+int32(t/tStep),2), y(2+int32(t/tStep),3), '.r');
    pause(0.0002)
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold on;
end

Tp = Time(2);

%% Analytical solution


V0 = [0 0.9*c 0]';

Y = zeros([2+int32(T/tStep), 6]);
Y_prev = zeros([1,3]);
Y(1,4:6) = m*V0;

for t=tStep:tStep:T
    [T,V] = ode45(@(time,v) q*(E' + cross(v,B'))/m, [0 t], V0);
    vel = V(end,:);
    Y(1+int32(t/tStep), 1:3) = Y_prev + vel*tStep;   
    Y_prev = Y(1+int32(t/tStep), 1:3);   
    Y(1+int32(t/tStep), 4:6) = m*vel;
    
end

err_percent = 100*(y-Y)./Y;   % Percentage error between integrator method and analytical integration method  (for 0 values, it will show NaN)