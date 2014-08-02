function pfLandkid

% PRELIMINARIES
% Clear the display.
clear;
clc;
% Increase this for slower animation.
ptime=0.05;

% PARAMETERS
% Mean and variance of initial state.
muX = zeros(2,1);
varX = 1e-0*eye(2);
% Mean and variance of process noise.
muV = zeros(2,1);
varV = 2e-3*eye(2);
% Mean and variance of measurement noise.
muW = 0.001;
varW = 1e-3;
% Number of time steps.
T = 250;
% Sample the initial state (normally distributed).
x(:,1) = [0.8;0];
% Number and location of the landmarks.
nLandmarks = 2;
xLandmarks = -1+2*rand(2,nLandmarks);

% SET UP PARTICLE FILTER
% Set the number of particles.
nPart = 500;
% Sample the particles from the prior (initial) distribution.
xPart = repmat(muX,1,nPart)+(sqrt(varX)*randn(2,nPart));
xPartMean = sum(xPart,2)/nPart;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP THE PLOTS
figure(1);
clf;
axis equal;
axis([-1 1 -1 1]);
box on;
hold on;
rTrueLine = line(x(1),x(2));
set(rTrueLine,'linestyle','-','color','b','linewidth',2);
rPartMeanLine = line(xPartMean(1),xPartMean(2));
set(rPartMeanLine,'linestyle','-','color','r','linewidth',1);
hLandmarks = line(xLandmarks(1,:),xLandmarks(2,:));
set(hLandmarks,'linestyle','none','color','k','marker','*','markersize',20);
rPart = line(xPart(1,:),xPart(2,:));
set(rPart,'linestyle','none','color','r','marker','.');
rPartMean = line(xPartMean(1),xPartMean(2));
set(rPartMean,'linestyle','none','markeredgecolor','g','markerfacecolor','w','linewidth',2,'marker','o','markersize',12);
rTrue = line(x(1),x(2));
set(rTrue,'linestyle','none','color','b','marker','.','markersize',30);
% Make sure everything is drawn.
drawnow;
pause(ptime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pause;

% SIMULATE THE SYSTEM, MAINTAINING THE POSTERIOR (I.E., THE BELIEF)
% Iterate over each time step.
for i=1:T  
    %%%%%%%%%%%%%%%%%%% SIMULATE %%%%%%%%%%%%%%%%%%%%%
    % Simulate measurement
    for j=1:nLandmarks
        % Calculate difference between robot position and landmark
        dx = x(:,i)-xLandmarks(:,j);
        d = sqrt(sum(dx.^2));
        % make the distance + mean + variance of measurement noise
        % the simulated measurement
        y(j,i) = d + (muW+(sqrt(varW)*randn));
    end
    
    if(i>=51)&&(i<60)
        pause;
    end
    
    theta = 2*pi*(((i+1)-1)/(T-1));
    
    % kidnap the robot
    if(i>50)
        x(:,i+1) = 0.4*[cos(theta);sin(theta)];
    else
        x(:,i+1) = 0.85*[cos(theta);sin(theta)];
    end
    
    
    %%%%%%%%%%%%%%%%%% MAINTAIN THE POSTERIOR %%%%%%%%%%%%
        
    % PARTICLE FILTER
    % 1. (Simulate) Propagate each sample according to the dynamic model
    xPartBar = xPart + repmat(muV,1,nPart) + sqrt(varV)*randn(2,nPart);
    % 2. (Weight) See how likely each particle is, given the measurement
    wPart = ones(1,size(xPartBar,2));
    for j=1:nLandmarks
        % Calculate difference between each particle and landmark positions
        dx = xPart-repmat(xLandmarks(:,j),1,nPart);
        d = sqrt(sum(dx.^2,1));
        % Multiple particles weight for each landmark
        wPart = wPart.*normpdf(d,y(j,i)+muW,sqrt(varW));
    end
    % 3. (Resample) Sample nPart particles again, according to their weights
    pPart = cumsum(wPart)/sum(wPart);
    for j=1:nPart
        u=rand;
        k=find(pPart>=u,1,'first');
        xPart(:,j) = xPartBar(:,k);
    end
    
    % COMPUTE MEANS FOR DISPLAY
    % Particle filter (average position of particles).
    xPartMean(:,i) = sum(xPart,2)/nPart;   
    
    %%%%%%%%%%%%%%%% DISPLAY THE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(rTrue,'xdata',x(1,i),'ydata',x(2,i));
    set(rPart,'xdata',xPart(1,:),'ydata',xPart(2,:));
    set(rPartMean,'xdata',xPartMean(1,i),'ydata',xPartMean(2,i));
    set(rPartMeanLine,'xdata',xPartMean(1,:),'ydata',xPartMean(2,:));
    set(rTrueLine,'xdata',x(1,:),'ydata',x(2,:));
    drawnow;
    pause(ptime);
    
end
