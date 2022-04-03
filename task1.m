clc;clear;

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
    % Runge Kutta 4th Order Numerical Integration
    k1 = tStep*f(v);
    k2 = tStep*f(v+(k1/2));
    k3 = tStep*f(v+(k2/2));
    k4 = tStep*f(v+k3);
    v = v + (k1 + 2*k2 + 2*k3 + k4)/6;
    y(2+int32(t/tStep), 1:3) = y_prev + v*tStep;   % x, y, z position after each time step
    y_prev = y(2+int32(t/tStep), 1:3);   %Setting previous timestep positioning of electron
    y(2+int32(t/tStep), 4:6) = m*v;    % x, y, z momentum after each time step
    if t~=0 && ~isequal(sign(y(2+int32(t/tStep)-1, 2)), sign(y(2+int32(t/tStep), 2)))  %Works perfectly for Electric field & Magnetic field other than in y-direction 
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

% plot3(y(:,1), y(:,2), y(:,3));
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
%xlim([-5e-8,5e-8])
% xlim([0,1]);
% ylim([-6e-5, 6e-5]);
% zlim([0, 1e-6]);

%% Calculating radius using Centrifugal Force equation

%Since Magnetic field (B) is the only responsible field for electron to
%move in a perfect circle and Electric field only makes the electron
%accelerate in the direction of E, centrifugal force has only one
%contributor of Magnetic force and not due to electric field.
Fe = q*E;
Fb = q*cross(v, B);
F = norm((Fb + Fe.*Fb),2);
V = norm(v,2); 
R = m*V^2/F;


%% Analytical solution


% V0 = [0 0.9*c 0]';
% 
% Y = zeros([2+int32(T/tStep), 6]);
% Y_prev = zeros([1,3]);
% Y(1,4:6) = m*V0;
% 
% for t=tStep:tStep:T
%     [T,V] = ode45(@(time,v) q*(E' + cross(v,B'))/m, [0 t], V0);
%     vel = V(end,:);
%     Y(1+int32(t/tStep), 1:3) = Y_prev + vel*tStep;   
%     Y_prev = Y(1+int32(t/tStep), 1:3);   
%     Y(1+int32(t/tStep), 4:6) = m*vel;
%     
% end
% 
% err_percent = 100*(y-Y)./Y;   % Percentage error between integrator method and analytical integration method  (for 0 values, it will show NaN)