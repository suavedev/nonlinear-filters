% #########################################
% GMU 2013
% ECE Scholarly Paper
% Cesar Corzo
% Mobile Robot Filtering Solutions
% Particle Filter vs EKF
% #########################################

% Cleaning up environment
clear all force;
close all force;

% Increase this for slower animation.
ptime=0.05;

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
qsi = 0.0025;
Q = [q 0 0;0 q 0;0 0 qsi];

% Extended Kalman Filter Parameters
% Static parameters
H = [1 0 0;0 1 0];

% Initial Conditions
x(:,1) = zeros(3,1);
Pplus = [r 0 0;0 r 0;0 0 0];
xest(:,1) = x(:,1);
z(:,1) = H*x(:,1) + sqrt(r)*randn(2,1);
xesterror(:,1) = xest(:,1) - x(:,1);
sensorerror(:,1) = z(:,1) - H*x(:,1);

% SET UP PARTICLE FILTER
% Set the number of particles.
nPart = 500;
% Sample the particles from the prior (initial) distribution.
xPart = repmat(x(:,1),1,nPart)+(delta*sqrt(r)*randn(3,nPart));
xPartMean = sum(xPart,2)/nPart;
wPart = zeros(1,nPart);
xPart_N(:,:) = zeros(3,nPart);
z(:,1) = H*x(:,1) + sqrt(r)*randn(2,1);
xPartMeanerror(:,1) = xPartMean(:,1) - x(:,1);
%sensorerror(:,1) = z(:,1) - H*x(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP THE PLOTS
figure(1)
clf;
title('Robot Trajectory');
xlabel('x');
ylabel('y');

axis equal;
%axis([-1 11 -1 11]);
%box on;
hold on;
% True trajectory
%rTrue = line(x(1),x(2));
%set(rTrue,'linestyle','none','color','b','marker','.','markersize',20);
rTrueLine = line(x(1),x(2));
set(rTrueLine,'linestyle','--','color','b','linewidth',2);
% Particles
rPart = line(xPart(1,:),xPart(2,:));
set(rPart,'linestyle','none','color','r','marker','.');
% Particle Mean
%rPartMean = line(xPartMean(1),xPartMean(2));
%set(rPartMean,'linestyle','none','markeredgecolor','g','markerfacecolor','w','linewidth',2,'marker','o','markersize',12);
rPartMeanLine = line(xPartMean(1),xPartMean(2));
set(rPartMeanLine,'linestyle','-','color','r','linewidth',1);

% EKF
rEKFLine = line(xest(1), xest(2));
set(rEKFLine,'linestyle','-','color','g','linewidth',1);

% Make sure everything is drawn.
drawnow;
pause(ptime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause;

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
    
    G = delta*[0 0 -v/(L*(cos(alpha(k-1)))^2)]';
   
    % Nonlinear process model
    f = [-v*sin(x(3,k-1)); v*cos(x(3,k-1)); (v/L)*tan(alpha(k-1))];
    % True State
    x(:,k) = x(:,k-1) + delta*f;
    %x_N = x(:,k-1) + f*delta + delta*sqrt(q)*randn(3,1);
    % Measurement model
    z(:,k) = H*x(:,k) + [sqrt(r)*randn sqrt(r)*randn]';
    
    %% Particle Filter
    %% 1. (Simulate) Propagate each sample according to the dynamic model
    for i=1:nPart
        fPart = [-v*sin(xPart(3,i)); v*cos(xPart(3,i)); (v/L)*tan(alpha(k-1))];
        xPart_N(:,i) = xPart(:,i) + fPart*delta + G*sqrt(q)*randn;
        zPart = H*xPart_N(:,i);
        % Calculate the weights for each particle
        pdf = normpdf(z(:,k), zPart, sqrt(R*[1 1]'));
        % Multiply each probability
        wPart(i) = pdf(1)*pdf(2);%*pdf(3);
    end
    
    % Normalize to form a probability distribution (i.e. sum to 1).
    pPart = wPart/sum(wPart);
    
    %% 3. Resampling
    % From this new distribution, now we randomly sample from it 
    % to generate our new estimate particles
    cspPart = cumsum(pPart); % cumulative sum
    for i = 1:nPart
        u = rand;
        j = find(u <= cspPart,1);
        xPart(:,i) = xPart_N(:,j);
    end
    
    xPartMean(:,k) = mean(xPart, 2);   
    
    %%%%%%%%%%%%%%%%%%%%% Extended Kalman Filter %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dynamic Model
    fest = [-v*sin(xest(3,k-1)); v*cos(xest(3,k-1)); (v/L)*tan(alpha(k-1))];
    F = [1 0 -delta*v*cos(xest(3,k-1));0 1 -delta*v*sin(xest(3,k-1));0 0 1];
    
    % EKF Steps
    Pminus = F*Pplus*F' + G*q*G';
    K = Pminus*H'*inv(H*Pminus*H' + R);
    %Pplus = (eye(3)- K*H)*Pminus;
    Pplus=[eye(3) - K*H]*Pminus*[eye(3) - K*H]' + K*R*K';
    
    xpredict = xest(:,k-1) + fest*delta;
    xest(:,k) = xpredict + K*(z(:,k-1) - H*xpredict);
    
    % errors
    xPartMeanerror(:,k) = xPartMean(:,k) - x(:,k);
    sensorerror(:,k) = z(:,k) - H*x(:,k);
    xEKFerror(:,k) = xest(:,k) - x(:,k);
    
    %%%%%%%%%%%%%%%% DISPLAY THE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %set(rTrue,'xdata',x(1,k),'ydata',x(2,k));
    set(rTrueLine,'xdata',x(1,1:k),'ydata',x(2,1:k));
    set(rPart,'xdata',xPart(1,:),'ydata',xPart(2,:));
    %set(rPartMean,'xdata',xPartMean(1,k),'ydata',xPartMean(2,k));
    set(rPartMeanLine,'xdata',xPartMean(1,:),'ydata',xPartMean(2,:));
    set(rEKFLine,'xdata',xest(1,:),'ydata',xest(2,:));
    
    drawnow;
    legend('True','Particles','PF', 'EKF', 4);
    pause(ptime);
end
ekfmse = [ mean(xEKFerror(1,:).^2), mean(xEKFerror(2,:).^2), mean(xEKFerror(3,:).^2) ]';
pfmse = [ mean(xPartMeanerror(1,:).^2), mean(xPartMeanerror(2,:).^2), mean(xPartMeanerror(3,:).^2) ]';

figure(2)
plot(radtodeg(x(3,:)), '--b');
hold on;
plot(radtodeg(xPartMean(3,:)), 'k');
hold on;
plot(radtodeg(xest(3,:)), 'g');
%hold on;
%plot(z(3,:), 'or');
%title('Heading');
xlabel('k')
ylabel('Angle(deg)')
legend('\Psi (k)', '\Psi_P_a_r_t(k)','\Psi_E_K_F(k)', 'z_\Psi (k)',2);

figure(3)
scatter(xPartMeanerror(1,:),xPartMeanerror(2,:));
hold on;
scatter(xEKFerror(1,:),xEKFerror(2,:), 'g');
hold on;
scatter(sensorerror(1,:),sensorerror(2,:), 'r');
title('Error');
xlabel('x')
ylabel('y')
legend('\epsilon_P_a_r_t (k)', '\epsilon_E_K_F (k)','\epsilon_z (k)',1);

figure(4)
subplot(2,2,1)
plot(abs(xPartMeanerror(1,:)));
hold on;
plot(abs(xEKFerror(1,:)), 'g');
xlabel('k')
ylabel('Error(m)')
title('|\epsilon_x|')
legend('\epsilon_P_a_r_t (k)', '\epsilon_E_K_F (k)', 1);
subplot(2,2,2)
plot(abs(xPartMeanerror(2,:)));
hold on;
plot(abs(xEKFerror(2,:)), 'g');
xlabel('k')
ylabel('Error(m)')
title('|\epsilon_y|')
legend('\epsilon_P_a_r_t (k)', '\epsilon_E_K_F (k)', 1);
subplot(2,2,3.5)
plot(radtodeg(abs(xPartMeanerror(3,:))));
hold on;
plot(radtodeg(abs(xEKFerror(3,:))), 'g');
xlabel('k')
ylabel('Error(deg)')
title('|\epsilon_\Psi|')
legend('\epsilon_P_a_r_t (k)', '\epsilon_E_K_F (k)', 1);


