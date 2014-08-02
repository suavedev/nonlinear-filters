% #########################################
% GMU 2013
% ECE Scholarly Paper
% Cesar Corzo
% Mobile Robot Filtering Solutions
% Particle Filter
% #########################################

% Cleaning up environment
clear all force;
close all force;

% Increase this for slower animation.
ptime=0.1;

% Robot Parameters
L = 2;
v = 1; %m/s
T = 15;
delta = 0.1;
xdes = 10;
ydes = 10;
gain = 0.5;

% Measurement Covariances
r = 0.004;
rsi = 0.003;
R = [r 0;0 r];

% Process Covariances
q = 0.0025;
qsi = 0.0002;
Q = [q 0 0;0 q 0;0 0 qsi];

x(:,1) = zeros(3,1);
H = [1 0 0;0 1 0];

% Compute robot trajectory
% State vector (x) contains robot location (x, y) and Si variables.
%[t, x(:,:), Sides, perfect_alpha] = robotTrajectory(L, v, T, delta, xdes, ydes);

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
sensorerror(:,1) = z(:,1) - H*x(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP THE PLOTS
figure(1);
clf;
title('Robot Trajectory');
xlabel('x')
ylabel('y')
%axis equal;
%axis([-1 11 -1 11]);
box on;
hold on;
% True trajectory
rTrue = line(x(1),x(2));
set(rTrue,'linestyle','none','color','b','marker','.','markersize',20);
rTrueLine = line(x(1),x(2));
set(rTrueLine,'linestyle','-','color','b','linewidth',2);
% Particles
rPart = line(xPart(1,:),xPart(2,:));
set(rPart,'linestyle','none','color','r','marker','.');
% Particle Mean
%rPartMean = line(xPartMean(1),xPartMean(2));
%set(rPartMean,'linestyle','none','markeredgecolor','g','markerfacecolor','w','linewidth',2,'marker','o','markersize',12);
rPartMeanLine = line(xPartMean(1),xPartMean(2));
set(rPartMeanLine,'linestyle','-','color','r','linewidth',1);

% Make sure everything is drawn.
drawnow;
pause(ptime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause;

for k=2:T/delta
    Sides(k-1) = -atan2((xdes-xPartMean(1,k-1)),(ydes-xPartMean(2,k-1)));
    % Calculate alpha
    angleerror(k-1) = Sides(k-1)-xPartMean(3,k-1);
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
    f = [-v*sin(x(3,k-1)); v*cos(x(3,k-1)); (v/L)*tan(alpha(k-1)+ sqrt(q)*randn)];
     % True State
    x(:,k) = x(:,k-1) + delta*f;
    %x_N = x(:,k-1) + delta*f + G*sqrt(q)*randn;
    %x_N = x(:,k-1) + f*delta + delta*sqrt(Q)*randn(3,1);
    % Measurement model
    z(:,k) = H*x(:,k) + [sqrt(r)*randn sqrt(r)*randn]';
    
    %% Particle Filter
    %% 1. (Simulate) Propagate each sample according to the dynamic model
    for i=1:nPart
        fPart = [-v*sin(xPart(3,i)); v*cos(xPart(3,i)); (v/L)*tan(alpha(k-1))];
        xPart_N(:,i) = xPart(:,i) + fPart*delta + G*sqrt(q)*randn;
        % Update -  compute zPart estimate
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
    % errors
    xPartMeanerror(:,k) = xPartMean(:,k) - x(:,k);
    sensorerror(:,k) = z(:,k) - H*x(:,k);
    
    %%%%%%%%%%%%%%%% DISPLAY THE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(rTrue,'xdata',x(1,k),'ydata',x(2,k));
    set(rTrueLine,'xdata',x(1,1:k),'ydata',x(2,1:k));
    set(rPart,'xdata',xPart(1,:),'ydata',xPart(2,:));
    %set(rPartMean,'xdata',xPartMean(1,k),'ydata',xPartMean(2,k));
    set(rPartMeanLine,'xdata',xPartMean(1,:),'ydata',xPartMean(2,:));
    
    drawnow;
    pause(ptime);
end

figure(2)
plot(x(3,:), 'b');
hold on;
plot(xPartMean(3,:), 'g');
%hold on;
%plot(z(3,:), 'or');
title('Heading');
xlabel('Sample')
ylabel('Heading Angle')
legend('\Psi (k)','\Psi_e_s_t(k)',1);

figure(3)
scatter(xPartMeanerror(1,:),xPartMeanerror(2,:));
hold on;
scatter(sensorerror(1,:),sensorerror(2,:), 'r');
title('Error');
xlabel('x')
ylabel('y')
legend('\epsilon_P_a_r_t (k)','\epsilon_z (k)',1);

