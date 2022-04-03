%% Leapfrog method, taking momentum as midpoint of intervals of position

clc;clear

% Actual substep advancement of position and used in Leapfrog
tStep = 10^-12;  %This timestep corresponds to less accuracy of leapfrog 
%tStep = 10^-14;  %This timestep of more precision has high accuracy
T = 5*10^-10;   % 5*10^-11 with tStep=10^-12 gives 1 revolution of electron
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
f = @(p) q*(E + cross(p,B))/m;
% Calculating momentum at first midpoint using initial momentum (Taylor series upto first order)
p_mid = y(1,4:6) + 0.5*tStep*f(y(1,4:6));

for t=0:tStep:T
    % LeapFrog Numerical Integration Method    
    y(2+int32(t/tStep), 1:3) = y_prev + p_mid*tStep/m;   % x, y, z position after each time step
    y_prev = y(2+int32(t/tStep), 1:3);   %Setting previous timestep positioning of electron 
    y(2+int32(t/tStep), 4:6) = p_mid + 0.5*tStep*f(p_mid);    % x, y, z momentum after each time step
    p_mid = p_mid + f(y(2+int32(t/tStep), 4:6))*tStep;   % Updating midpoint momentum for next timestep
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

Tp = Time(2);  %Time Period