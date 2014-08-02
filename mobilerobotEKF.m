% #########################################
% GMU 2013
% ECE Scholarly Paper
% Cesar Corzo
% Mobile Robot Filtering Solutions
% Extended Kalman Filter
% #########################################

% Cleaning up environment
clear all force;
close all force;

% Robot Parameters
gain = 0.5;
L = 2;
v = 1; %m/s
T = 15;
delta = 0.1;
xdes = 10;
ydes = 10;

% Measurement Covariances
r = 0.004;
R = [r 0;0 r];

% Process Covariances
q = 0.0025;

%Static parameters
H = [1 0 0;0 1 0];%; 0 0 1];

% Initial Conditions
x(:,1) = zeros(3,1);
Pplus = [r 0 0;0 r 0;0 0 0];
xest(:,1) = x(:,1);
z(:,1) = H*x(:,1) + [sqrt(r)*randn sqrt(r)*randn]';
xesterror(:,1) = xest(:,1) - x(:,1);
sensorerror(:,1) = z(:,1) - H*x(:,1);
Sides(1) = -atan2((ydes-x(2,1)),(xdes-x(1,1)));
angleerror(1) = Sides(1)-x(3,1);

% Discrete EKF
for k=2:T/delta
    Sides(k-1) = -atan2((xdes-xest(1,k-1)),(ydes-xest(2,k-1)));
    % Calculate alpha
    angleerror(k-1) = Sides(k-1)-xest(3,k-1);
    if(angleerror(k-1) > 2*pi)
       angleerror(k-1) = angleerror(k-1) - 2*pi;
    elseif(angleerror(k-1) < -2*pi)
       angleerror(k-1) = angleerror(k-1) + 2*pi;
    end

    if(gain*abs(angleerror(k-1)) > pi/4)
       alpha(k-1)= (pi/4)*sign(angleerror(k-1));
    else
       alpha(k-1)= gain*(angleerror(k-1));
    end
    
    G = delta*[0 0 (v/(L)*(1/(cos(alpha(k-1)))^2))]'
    %Bu = [0 0 (-v/(L)*(1/(cos(alpha(k-1)))^2))]';
    
    % Dynamic Model
    f = [-v*sin(x(3,k-1)); v*cos(x(3,k-1)); (v/L)*tan(alpha(k-1)+ sqrt(q)*randn)];
    fest = [-v*sin(xest(3,k-1)); v*cos(xest(3,k-1)); (v/L)*tan(alpha(k-1))];
    % True State
    x(:,k) = x(:,k-1) + delta*f;
    F = [1 0 -delta*v*cos(xest(3,k-1));0 1 -delta*v*sin(xest(3,k-1));0 0 1];

    % Measurement model
    z(:,k) = H*x(:,k) + [sqrt(r)*randn sqrt(r)*randn]';%x_N + sqrt(R)*randn(3,1);
    
    % EKF Steps
    Pminus = F*Pplus*F' + G*q*G';
    K = Pminus*H'*inv(H*Pminus*H' + R);
    %Pplus = (eye(3)- K*H)*Pminus;
    Pplus=[eye(3) - K*H]*Pminus*[eye(3) - K*H]' + K*R*K';
    
    xpredict = xest(:,k-1) + fest*delta;
    xest(:,k) = xpredict + K*(z(:,k-1) - H*xpredict);

    % errors
    xesterror(:,k) = xest(:,k) - x(:,k);
    sensorerror(:,k) = z(:,k) - H*x(:,k);
    rmserror(:,k) = rms((xest(:,k) - x(:,k)),3);
end

% Plot robot Trajectory
figure(1)
plot(x(1,:),x(2,:));
hold on;
%plot(z(1,:),z(2,:), 'or');
%hold on;
plot(xest(1,:),xest(2,:), 'g');
title('Robot Trajectory');
xlabel('x')
ylabel('y')

figure(2)
plot(x(3,:), 'b');
hold on;
plot(xest(3,:), 'g');
%hold on;
%plot(z(3,:), 'or');
title('Heading');
xlabel('Sample')
ylabel('Heading Angle')
legend('\Psi (k)','\Psi_e_s_t(k)',1);

figure(3)
scatter(xesterror(1,:),xesterror(2,:));
hold on;
scatter(sensorerror(1,:),sensorerror(2,:), 'r');
title('Error');
xlabel('x')
ylabel('y')
legend('\epsilon_e_s_t (k)','\epsilon_z (k)',1);

figure(4)
plot(rmserror(1,:))
title('RMS Error');
xlabel('k')
ylabel('RMS(x)')

