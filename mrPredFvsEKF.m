% #########################################
% GMU 2013
% ECE Scholarly Paper
% Cesar Corzo
% Mobile Robot Filtering Solutions
% Predictive Filter Estimation vs Extended Kalman Filter
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

% Static parameters
H = [1 0 0;0 1 0];

% Initial Conditions
w1 = 0.73;
% Randomly generate W
%A = randn(2);
%[U,ignore] = eig((A+A')/2); % (Don't really have to divide by 2)
%W = U*diag(abs(randn(3,1)))*U';
%W=w1*[0.3903 0.0012;0.0012 0.3903];
W=w1*[0.3903 0;0 0.3903];
x(:,1) = zeros(3,1);
Pplus = [r 0 0;0 r 0;0 0 0];
xest(:,1) = x(:,1);
xest_EKF(:,1) = x(:,1);
y_sensor(:,1) = H*x(:,1) + sqrt(r)*randn(2,1);
yest(:,1) = H*xest(:,1) + sqrt(r)*randn(2,1);
xesterror(:,1) = rms((xest(:,1) - x(:,1)),3);
xEKFerror(:,1) = rms((xest_EKF(:,1) - x(:,1)),3); 
sensorerror(:,1) = y_sensor(:,1) - H*x(:,1);
Sides(1) = -atan2((ydes-x(2,1)),(xdes-x(1,1)));
angleerror(1) = Sides(1)-x(3,1);
Sides_EKF(1) = -atan2((ydes-x(2,1)),(xdes-x(1,1)));
angleerror_EKF(1) = Sides_EKF(1)-x(3,1);

% Predictive Filter
for k=2:T/delta
    % change trajectory
    if(xest(1,k-1) > 4)
        xdes = 0;
    end
    % Calculate alpha for Predictive filter
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
    % Calculate alpha for EKF
    Sides_EKF(k-1) = -atan2((xdes-xest_EKF(1,k-1)),(ydes-xest_EKF(2,k-1)));
    
    angleerror_EKF(k-1) = Sides_EKF(k-1)-xest(3,k-1);
    if(angleerror_EKF(k-1) > 2*pi)
       angleerror_EKF(k-1) = angleerror_EKF(k-1) - 2*pi;
    elseif(angleerror_EKF(k-1) < -2*pi)
       angleerror_EKF(k-1) = angleerror_EKF(k-1) + 2*pi;
    end

    if(gain*abs(angleerror_EKF(k-1)) > pi/4)
       alpha_EKF(k-1)= (pi/4)*sign(angleerror_EKF(k-1));
    else
       alpha_EKF(k-1)= gain*(angleerror_EKF(k-1));
    end
    
    G_EKF = delta*[0 0 -v/(L*(cos(alpha(k-1)))^2)]';
    G = delta*[1 0;0 1;0 0];
    
    % Dynamic Model
    f = [-v*sin(x(3,k-1)); v*cos(x(3,k-1)); (v/L)*tan(alpha(k-1)+ sqrt(q)*randn)];
    fest = [-v*sin(xest(3,k-1)); v*cos(xest(3,k-1)); (v/L)*tan(alpha(k-1))];
    %fest = [0; v*cos(xest(3,k-1)); (v/L)*tan(alpha(k-1))];
    
    % True State
    x(:,k) = x(:,k-1) + delta*f;
    %Phi = [1 0 -delta*v*cos(xest(3,k-1));0 1 -delta*v*sin(xest(3,k-1));
    %       0 0 1-delta*v/(L*(cos(alpha(k-1)))^2)];

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
    d(:,k-1) = delta*part1*part2;
    xest(:,k) = xest(:,k-1) + fest*delta + G*d(:,k-1);
    
    %%%%%%%%%%%%%%%%%%%%% Extended Kalman Filter %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dynamic Model
    fest_EKF = [-v*sin(xest_EKF(3,k-1)); v*cos(xest_EKF(3,k-1)); (v/L)*tan(alpha_EKF(k-1))];
    F = [1 0 -delta*v*cos(xest_EKF(3,k-1));0 1 -delta*v*sin(xest_EKF(3,k-1));0 0 1];
    
    % EKF Steps
    Pminus = F*Pplus*F' + G_EKF*q*G_EKF';
    K = Pminus*H'*inv(H*Pminus*H' + R);
    %Pplus = (eye(3)- K*H)*Pminus;
    Pplus=[eye(3) - K*H]*Pminus*[eye(3) - K*H]' + K*R*K';
    
    xpredict =xest_EKF(:,k-1) + fest_EKF*delta;
    xest_EKF(:,k) = xpredict + K*(y_sensor(:,k-1) - H*xpredict);
    
    % errors
    xEKFerror(:,k) = xest_EKF(:,k) - x(:,k); 
    sensorerror(:,k) = y_sensor(:,k) - H*x(:,k);
    xesterror(:,k) = xest(:,k) - x(:,k);
end
ekfmse = [ mean(xEKFerror(1,:).^2), mean(xEKFerror(2,:).^2), mean(xEKFerror(3,:).^2) ]'
predfmse = [ mean(xesterror(1,:).^2), mean(xesterror(2,:).^2), mean(xesterror(3,:).^2) ]'
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
% figure(1)
plot(x(1,:),x(2,:), '--b');
hold on;
%plot(y_sensor(1,:),y_sensor(2,:), 'or');
%hold on;
plot(xest(1,:),xest(2,:), 'k');
hold on;
plot(xest_EKF(1,:),xest_EKF(2,:), 'g');
title('Robot Trajectory');
xlabel('x')
ylabel('y')
legend('True','PredF','EKF', 2);

figure(2)
plot(radtodeg(x(3,:)), '--b');
hold on;
plot(radtodeg(xest_EKF(3,:)), 'g');
hold on;
plot(radtodeg(xest(3,:)), 'k');
%hold on;
%plot(y_sensor(3,:), 'or');
title('Heading (\Psi)');
xlabel('k')
ylabel('Angle (deg)')
legend('True','EKF', 'PredF', 2);

figure(3)
scatter(xesterror(1,:),xesterror(2,:),'k');
hold on;
scatter(xEKFerror(1,:),xEKFerror(2,:), 'g');
hold on;
scatter(sensorerror(1,:),sensorerror(2,:), 'r');
title('Error Scatter Plot');
xlabel('x')
ylabel('y')
legend('\epsilon_P_r_e_d_F (k)', '\epsilon_E_K_F (k)','\epsilon_z (k)',1);

figure(4)
subplot(2,2,1)
plot(abs(xesterror(1,:)));
hold on;
plot(abs(xEKFerror(1,:)), 'g');
xlabel('k')
ylabel('Error (m)')
title('|\epsilon_x|')
legend('\epsilon_P_r_e_d (k)', '\epsilon_E_K_F (k)', 1);
subplot(2,2,2)
plot(abs(xesterror(2,:)));
hold on;
plot(abs(xEKFerror(2,:)), 'g');
xlabel('k')
ylabel('Error (m)')
title('|\epsilon_y|')
legend('\epsilon_P_r_e_d (k)', '\epsilon_E_K_F (k)', 1);
subplot(2,2,3.5)
plot(radtodeg(abs(xesterror(3,:))));
hold on;
plot(radtodeg(abs(xEKFerror(3,:))), 'g');
xlabel('k')
ylabel('Error (deg)')
title('|\epsilon_\Psi|')
legend('\epsilon_P_r_e_d_F (k)', '\epsilon_E_K_F (k)', 1);
