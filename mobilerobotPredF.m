% #########################################
% GMU 2013
% ECE Scholarly Paper
% Cesar Corzo
% Mobile Robot Filtering Solutions
% Predictive Filter Estimation
% #########################################

% Cleaning up environment
clear all force;
close all force;

% Robot Parameters
gain = 0.5;
L = 2;
v = 1; %m/s
T = 17;
delta = 0.1;
xdes = 10;
ydes = 10;

% Measurement Covariances
r = 0.004;
R = [r 0;0 r];

% Process Covariances
q = 0.0025;

% Static parameters
H = [1 0 0;0 1 0];

% Initial Conditions
w1 = 0.73;
% Randomly generate W
%A = randn(3);
%[U,ignore] = eig((A+A')/2); % (Don't really have to divide by 2)
%W = U*diag(abs(randn(3,1)))*U';
%W=w1*[0.3903 0.0012;0.0012 0.3903];
W=w1*[0.3903 0;0 0.3903];
%W = w1*W;
x(:,1) = zeros(3,1);
Pplus = R;
xest(:,1) = x(:,1);
y_sensor(:,1) = H*x(:,1) + sqrt(r)*randn(2,1);
yest(:,1) = H*xest(:,1) + sqrt(r)*randn(2,1);
xesterror(:,1) = xest(:,1) - x(:,1);
sensorerror(:,1) = y_sensor(:,1) - H*x(:,1);
Sides(1) = -atan2((ydes-x(2,1)),(xdes-x(1,1)));
angleerror(1) = Sides(1)-x(3,1);

% Predictive Filter
for k=2:T/delta
    % change trajectory
    if(xest(1,k-1) > 4)
        xdes = 0;
    end
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
    
    %G = delta*[0 0 -v/(L*(cos(alpha(k-1)))^2)]';
    G = delta*[1 0;0 1;0 0];
    %G = delta*[0 0 1]';
    
    % Dynamic Model
    f = [-v*sin(x(3,k-1)); v*cos(x(3,k-1)); (v/L)*tan(alpha(k-1)+ sqrt(q)*randn)];
    fest = [-v*sin(xest(3,k-1)); v*cos(xest(3,k-1)); (v/L)*tan(alpha(k-1))];
    %fest = [0; v*cos(xest(3,k-1)); (v/L)*tan(alpha(k-1))];
    
    % True State
    x(:,k) = x(:,k-1) + delta*f;
    Phi = [1 0 -delta*v*cos(xest(3,k-1));0 1 -delta*v*sin(xest(3,k-1));
            0 0 1-delta*v/(L*(cos(alpha(k-1)))^2)];

    % Measurement model
    y_sensor(:,k) = H*x(:,k) + [sqrt(r)*randn sqrt(r)*randn]';
    yest(:,k-1) = H*xest(:,k-1);
    
    % Calculating d(t)
    % pi = 1, then
    z = delta*H*fest;
    Lambda = delta*eye(2);
    S = H*G;
    part1 = -inv((Lambda*S)'*inv(R)*Lambda*S + W)*(Lambda*S)'*inv(R);
    part2 = z - y_sensor(:,k-1) + yest(:,k-1);
    %part1 = eye(2)-inv(W)*G'*H'*inv(H*G*inv(W)*G'*H' + R)*H*G;
    %part2 = inv(W)*G'*H'*inv(R)*(y_sensor(:,k-1) + H*Phi*xest(:,k-1));
    d(:,k-1) = delta*part1*part2;
    xest(:,k) = xest(:,k-1) + fest*delta + G*d(:,k-1);
    % errors
    xesterror(:,k) = xest(:,k) - x(:,k);
    sensorerror(:,k) = y_sensor(:,k) - H*x(:,k);
    rmserror(:,k) = rms((xest(:,k) - x(:,k)),3);
end

% Covariance constraint to help compute W
% ebar = [mean(y_sensor(1,1:149)-yest(1,:)) mean(y_sensor(2,1:149)-yest(2,:))]';
% e = y_sensor(:,1:149)-yest(:,:);
% sum = 0;
% for i=1:149
%     a = (e(:,i)- ebar)*(e(:,i)- ebar)';
%     sum = sum + a;
% end
% sum = sum/149;
% Plot robot Trajectory
figure(1)
plot(x(1,:),x(2,:));
hold on;
%plot(y_sensor(1,:),y_sensor(2,:), 'or');
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
%plot(y_sensor(3,:), 'or');
title('Heading');
xlabel('Sample')
ylabel('Heading Angle')
legend('\Psi (k)','\Psi_e_s_t(k)',1);

figure(3)
scatter(xesterror(1,:),xesterror(2,:));
hold on;
scatter(sensorerror(1,:),sensorerror(2,:), 'r');
title('Error Scatter Plot');
xlabel('x')
ylabel('y')
legend('\epsilon_e_s_t (k)','\epsilon_z (k)',1);

% figure(4)
% plot(d(1,:));

figure(5)
plot(rmserror(1,:))
title('RMS Error');
xlabel('k')
ylabel('RMS(x)')
