function OneDMotion

% PRELIMINARIES
% Clear the display.
clear;
clc;
% Increase this for slower animation.
ptime=0.1;

% PARAMETERS
% Mean and variance of initial state.
muX = -0.5;
varX = 1e-1;
% Mean and variance of process noise.
muV = 1e-2;
varV = 5e-3;
% Mean and variance of measurement noise.
muW = 0.0;
varW = 5e-2;
% Number of time steps.
T = 100;
% Bounds on x, the state.
xmin = -1;
xmax = 1;
% Sample the initial state (normally distributed).
x(1) = sqrt(varX)*randn(1,1)+muX;

% SET UP HISTOGRAM FILTER
% Number of bins.
n = 31;
% Size of each bin.
dx = (xmax-xmin)/n;
% Compute prior (initial belief).
X = linspace(xmin,xmax,n);
PX = dx*normpdf(X,muX,sqrt(varX));
PX = (PX/sum(PX));
% Compute state transition model and measurement model.
for i=1:n
    xcur = X(i);
    PXX{i} = dx*normpdf(X,xcur+muV,sqrt(varV));
    PXX{i} = PXX{i}/sum(PXX{i});
    PYX{i} = dx*normpdf(X,xcur+muW,sqrt(varW));
    PYX{i} = PYX{i}/sum(PYX{i});
end

% SET UP KALMAN FILTER
xbar = muX;
rbar = varX;

% SET UP PARTICLE FILTER
% Set the number of particles.
nPart = 100;
% Sample the particles from the prior (initial) distribution.
xPart = muX + sqrt(varX)*randn(nPart,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP THE PLOTS (YOU CAN SAFELY IGNORE THIS)
% First window.
figure(1);
clf;
hMainPlot = subplot('Position',[0.2 0.45 0.75 0.5]);
title('probability density functions');
ymax=6*max(PX/dx);
%axis([xmin xmax 0 ymax]);
set(gca,'xtick',[]);
box on;
hold on;
hParticlePlot = subplot('Position',[0.2 0.375 0.75 0.05]);
axis([xmin xmax -1 1]);
box on;
set(gca,'xtick',[],'ytick',[0]);
set(gca,'yticklabel',{'particles'});
hMeansPlot = subplot('Position',[0.2 0.1 0.75 0.25]);
xlabel('position');
axis([xmin xmax 0 6]);
box on;
set(gca,'xtick',[xmin:(xmax-xmin)/10:xmax],'ytick',[]);
subplot(hMainPlot);
hHistogram=bar(X,PX/dx,'barwidth',1,'facecolor','r');
Xfine=linspace(-1,1,1000);
hKalman=line(Xfine,normpdf(Xfine,xbar,sqrt(rbar)));
set(hKalman,'linewidth',2,'color','b');
subplot(hParticlePlot);
tmpHeight=0.7;
tmpX=reshape([xPart';xPart';inf(size(xPart'))],[],1);
tmpY=reshape([tmpHeight*ones(size(xPart'));-tmpHeight*ones(size(xPart'));inf(size(xPart'))],[],1);
hParticle=line(tmpX,tmpY);
set(hParticle,'linestyle','-','color','k','linewidth',0.1);
subplot(hMeansPlot);
hTrueMeanLine = line([0 0],[0 6]);
set(hTrueMeanLine,'color','g','linewidth',2);
hTrueMean = line(0,5);
set(hTrueMean,'marker','o','markersize',12,'markeredgecolor','k','markerfacecolor','g');
hMeasMean = line(0,4);
set(hMeasMean,'marker','o','markersize',12,'markeredgecolor','k','markerfacecolor','w');
hHistMean = line(0,3);
set(hHistMean,'marker','o','markersize',12,'markeredgecolor','k','markerfacecolor','r');
hKalmMean = line(0,2);
set(hKalmMean,'marker','o','markersize',12,'markeredgecolor','k','markerfacecolor','b');
hPartMean = line(0,1);
set(hPartMean,'marker','o','markersize',12,'markeredgecolor','k','markerfacecolor','k');
set(gca,'ytick',1:5,'YTickLabel',{'particles mean';'kalman mean';'histogram mean';'measurement';'true mean'})
set(hTrueMean,'xdata',x(1));
set(hTrueMeanLine,'xdata',[x(1) x(1)]);
set(hMeasMean,'xdata',inf);
set(hHistMean,'xdata',sum(PX.*X));
set(hKalmMean,'xdata',xbar);
set(hPartMean,'xdata',mean(xPart));
% Second window.
figure(2);
clf;
axis([1 T xmin xmax]);
box on;
hold on;
hTrueMean2 = line(inf,inf);
set(hTrueMean2,'color','g','linewidth',4);
hMeasMean2 = line(inf,inf);
set(hMeasMean2,'marker','o','markersize',6,'markeredgecolor','k','markerfacecolor','w','linestyle','none');
hHistMean2 = line(inf,inf);
set(hHistMean2,'color','r','linewidth',2);
hKalmMean2 = line(inf,inf);
set(hKalmMean2,'color','b','linewidth',2);
hPartMean2 = line(inf,inf);
set(hPartMean2,'color','k','linewidth',2);
xlabel('Time');
ylabel('Position');
% Make sure everything is drawn.
drawnow;
pause(ptime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% SIMULATE THE SYSTEM, MAINTAINING THE POSTERIOR (I.E., THE BELIEF)
% Iterate over each time step.
for i=1:T
    
    
    %%%%%%%%%%%% SIMULATE
    
    % Get the current measurement (and enforce that it is in [xmin,xmax]).
    y(i) = x(i) + ((sqrt(varW)*randn)+muW);
    if (y(i)<xmin) y(i)=xmin+1e-6;
    elseif (y(i)>xmax) y(i)=xmax-1e-6;
    end

    % Get the next state (and enforce that it is in [xmin,xmax]).
    x(i+1) = x(i) + ((sqrt(varV)*randn)+muV);
    if (x(i+1)<xmin) x(i+1)=xmin;
    elseif (x(i+1)>xmax) x(i+1)=xmax;
    end
    

    %%%%%%%%%%%% MAINTAIN THE POSTERIOR
        
    % UPDATE HISTOGRAM FILTER
    % Get the bin into which the current measurement fell.
    iMeas = ceil((y(i)-xmin)/dx);
    % Update the posterior.
    % 1. Prediction.
    PXpred = zeros(size(PX));
    for j=1:n
        PXpred = PXpred + (PXX{j}.*PX(j));
    end
    % 2. Correction.
    PXcorr = PYX{iMeas}.*PXpred;
    % 3. Normalize the result.
    PX = PXcorr/sum(PXcorr);
    
    % KALMAN FILTER
    if (i~=1)
        % 1. Time update.
        xbar = xhat(i-1)+muV;
        rbar = rhat(i-1)+varV;
    end
    % 2. Measurement update.
    k = rbar/(rbar+varW);
    xhat(i) = xbar + k*(y(i)-xbar);
    rhat(i) = (1-k)*rbar;
    
    % PARTICLE FILTER
    % 1. (Simulate) Propagate each sample according to the dynamic model
    xPartBar = xPart + (muV + sqrt(varV)*randn(nPart,1));
    % 2. (Weight) See how likely each particle is, given the measurement
    wPart = normpdf(xPartBar,y(i)+muW,sqrt(varW));
    % 3. (Resample) Sample nPart particles again, according to their weights
    pPart = cumsum(wPart)/sum(wPart);
    for j=1:nPart
        u=rand;
        k=find(pPart>=u,1,'first');
        xPart(j) = xPartBar(k);
    end
    
    % COMPUTE MEANS FOR DISPLAY
    % Histogram (weighted average of bins).
    histMean = sum(PX.*X);
    % Kalman filter (just the current state estimate).
    kalmMean = xhat(i);
    % Particle filter (average position of particles).
    partMean = mean(xPart);
    
    
    %%%%%%%%%%%% DISPLAY THE RESULTS (YOU CAN SAFELY IGNORE THIS)
    set(hHistogram,'xdata',X,'ydata',PX/dx);
    set(hKalman,'xdata',Xfine,'ydata',normpdf(Xfine,xhat(i),sqrt(rhat(i))));
    tmpHeight=0.7;
    tmpX=reshape([xPart';xPart';inf(size(xPart'))],[],1);
    tmpY=reshape([tmpHeight*ones(size(xPart'));-tmpHeight*ones(size(xPart'));inf(size(xPart'))],[],1);
    set(hParticle,'xdata',tmpX,'ydata',tmpY);
    set(hTrueMean,'xdata',x(i));
    set(hTrueMeanLine,'xdata',[x(i) x(i)]);
    set(hMeasMean,'xdata',y(i));
    set(hHistMean,'xdata',histMean);
    set(hKalmMean,'xdata',kalmMean);
    set(hPartMean,'xdata',partMean);
    set(hTrueMean2,'xdata',1:i,'ydata',x(1:i));
    set(hMeasMean2,'xdata',1:i,'ydata',y);
    ytmp = get(hHistMean2,'ydata');
    set(hHistMean2,'xdata',0:i,'ydata',[ytmp histMean]);
    ytmp = get(hKalmMean2,'ydata');
    set(hKalmMean2,'xdata',0:i,'ydata',[ytmp kalmMean]);
    ytmp = get(hPartMean2,'ydata');
    set(hPartMean2,'xdata',0:i,'ydata',[ytmp partMean]);
    drawnow;
    pause(ptime);
    
end