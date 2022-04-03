%% Dormant-Prince745 Method

clc;clear

tStep = 10^-14;
T = 5*10^-11;   % 5*10^-11 with tStep=10^-12 gives 2 revolutions of electron
E = [0 0 0];        % in Newton/Coulomb
B = [0 0 1];        % in Tesla
c = 3*10^8;         % in m/s
m = 9.109*10^-31;   % in kgs
q = -1.602*10^-19;   % in Coulomb
i = 0;
eps = 1e-7;
% Initial Conditions
v = [0 0.9*c 0];    % initial 3D velocity
y = zeros([2+int32(T/tStep), 6]);   % x, y, z, px, py, pz
zv = zeros([2+int32(T/tStep), 6]);
y_prev = zeros([1,3]);
zv_prev = zeros([1,3]);
y(1,4:6) = m*v;
f = @(v) q*(E + cross(v,B))/m;
z = v;

for t=0:tStep:T
    % DoPri5 Numerical Integration
    k1 = tStep*f(v);
    k2 = tStep*f(v+k1/5);
    k3 = tStep*f(v+(3*k1+9*k2)/40);
    k4 = tStep*f(v+44*k1/45-56*k2/15+32*k3/9);
    k5 = tStep*f(v+19372*k1/6561-25360*k2/2187+64448*k3/6561-212*k4/729);
    k6 = tStep*f(v+9017*k1/3168-355*k2/33-46732*k3/5247+49*k4/176-5103*k5/18656);
    k7 = tStep*f(v+35*k1/384+500*k3/1113+125*k4/192-2187*k5/6784+11*k6/84);
    v = v + (1921409*k1 + 9690880*k3 + 13122270*k4 - 5802111*k5 + 1902912*k6 + 534240*k7)/21369600;  %Using 4th order RK
    z = z + (5179*k1/57600 + 7571*k3/16695 + 393*k4/640 - 92097*k5/339200 + 187*k6/2100 + 1*k7/40);  %Using 5th order RK
    y(2+int32(t/tStep), 1:3) = y_prev + v*tStep;   % x, y, z position after each time step
    y_prev = y(2+int32(t/tStep), 1:3);   %Setting previous timestep positioning of electron
    y(2+int32(t/tStep), 4:6) = m*v;    % x, y, z momentum after each time step
    
    zv(2+int32(t/tStep), 1:3) = zv_prev + z*tStep;
    zv_prev = zv(2+int32(t/tStep), 1:3);
    zv(2+int32(t/tStep), 4:6) = m*z;
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
    % Update tStep using doPri5
    err = norm(abs(zv(1:3)-y(1:3)),2);   %z-y
    if norm((zv(1:3)-y(1:3)),2)<eps*norm(y(1:3),2) && norm((zv(4:6)-y(4:6)),2)<eps*norm(y(4:6),2) && err~=0  %position_error < eps_step * length ; and momentum error < eps_step * mag(momentum)
        tStep = ((eps*tStep/(2*err))^(1/4))*tStep;
    end
end

Tp = Time(2);