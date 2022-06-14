% This script gives examples of using simDIRECT with single/multiple
% objectives, and with/without nonlinear constraints.  The test problems 
% (objectives and constraints) are also in this file, after the examples. 
%
% Version 1.1 (June 13, 2022).

% Start with a clean slate

close all;
clearvars;
clc;

%**************************************************************************
% Example 1 
% Test Problem:     Hartman 6  
% Characteristics:  1 objective, 0 constraints, 6 variables
%**************************************************************************

% Set inputs for optimizing the "Hartman 6" test problem

fun   = @(x) Hartman6(x);
con   = [];
xl    = zeros(1,6);
xu    = ones(1,6);
feps  = 0.0001;
fu    = inf;
maxfn = 150;

% Optimize with simDIRECT 

output = simDIRECT(fun,con,xl,xu,feps,fu,maxfn);

% Report the optimum and plot the convergence history

fprintf('Example 1:  Hartman 6\n');
if ~any(output.isPareto)
    fprintf('No feasible solution found!\n');
else
    fprintf('Found feasible solution:\n');
    % Get index of 1st Pareto pt (may be several if multiple tied minima)
    [~,i] = max(output.isPareto); 
    % Report the solution
    xmin = output.x(i,:);
    fmin = output.f(i);
    fprintf('xmin = [%f, %f, %f, %f, %f, %f]\n',...
                     xmin(1),xmin(2),xmin(3),xmin(4),xmin(5),xmin(6));
    fprintf('fmin = %f\n',fmin)
    % Plot convergence history
    fglobal = -3.3223680114155152; % known global min for Hartman6
    neval = 1:maxfn;
    fmin  = zeros(1,maxfn);
    fmin(1) = output.f(1);
    for i=2:maxfn
        fmin(i) = min(output.f(i),fmin(i-1));
    end
    figure
    plot(neval,fmin,'k-','LineWidth',2,'DisplayName','simDIRECT');
    hold on;
    plot([0,maxfn],[fglobal,fglobal],'k:','DisplayName','global optimum',...
        'LineWidth',2);
    axes = gca;
    axes.FontName = 'Times';
    axes.FontSize = 16;
    axes.PlotBoxAspectRatio = [ 1.5,1,1 ];
    axes.XLim = [0,150];
    axes.XTick = 0:25:150;   
    xlabel('Number Function Evaluations','Interpreter','latex',...
        'FontSize', 18);
    ylabel({'Best','Objective','Value'},'Interpreter','latex', ...
        'FontSize', 18, 'Rotation',0,'HorizontalAlignment','right',...
        'VerticalAlignment','middle');
    legend('EdgeColor','none','Interpreter','latex','FontSize',18); 
    title('Iteration History for Hartman 6','Interpreter','latex','FontSize',20);
end
keyboard


%**************************************************************************
% Example 2 
% Test Problem:     NASA speed reducer
% Characteristics:  1 objective, 11 constraints, 7 variables
%**************************************************************************

% Set inputs for the "NASA speed reducer" test problem

fun   = @(x) fun_NASA_speed_reducer(x);
con   = @(x) con_NASA_speed_reducer(x);
xl    = [2.6, 0.7, 17, 7.3, 7.8, 2.9, 5.0];
xu    = [3.6, 0.8, 28, 8.3, 8.3, 3.9, 5.5];
feps  = 0.01;
fu    = inf;
maxfn = 10000;

% Optimize with simDIRECT 

output = simDIRECT(fun,con,xl,xu,feps,fu,maxfn);

% Report the optimum and plot the convergence history

fprintf('Example 2:  NASA speed reducer\n');
if ~any(output.isPareto)
    fprintf('No feasible solution found!\n');
else
    fprintf('Found feasible solution:\n');
    % Get index of 1st Pareto pt (may be several if multiple tied minima)
    [~,i] = max(output.isPareto); 
    % Report the solution
    xmin = output.x(i,:);
    fmin = output.f(i);
    fprintf('xmin = [%f, %f, %f, %f, %f, %f]\n',...
                     xmin(1),xmin(2),xmin(3),xmin(4),xmin(5),xmin(6));
    fprintf('fmin = %f\n',fmin)
    % Plot convergence history
    fglobal =  2996.3481692405475769; % known global min 
    neval = 1:maxfn;
    fmin  = zeros(1,maxfn);
    BestFeasObj = NaN;
    for i=1:maxfn
        if all(output.g(i,:) <= 0) 
            if isnan(BestFeasObj)
                BestFeasObj = output.f(i);
            elseif output.f(i) < BestFeasObj
                BestFeasObj = output.f(i);
            end
        end
        fmin(i) = BestFeasObj;   
    end
    figure
    plot(neval,fmin,'k-','LineWidth',2,'DisplayName','simDIRECT');
    hold on;
    plot([0,maxfn],[fglobal,fglobal],'k:','DisplayName','global optimum',...
        'LineWidth',2);
    axes = gca;
    axes.FontName = 'Times';
    axes.FontSize = 12;
    axes.PlotBoxAspectRatio = [ 1.5,1,1 ];
    axes.XLim = [0,maxfn];
    axes.XTick = 0:maxfn/5:maxfn;   
    xlabel('Number Function Evaluations','Interpreter','latex',...
        'FontSize', 18);
    ylabel({'Best  ','Objective  ','Value  '},'Interpreter','latex', ...
        'FontSize', 18, 'Rotation',0,'HorizontalAlignment','right',...
        'VerticalAlignment','middle');
    legend('EdgeColor','none','Interpreter','latex','FontSize',18); 
    title('Iteration History for NASA Speed Reducer','Interpreter','latex','FontSize',20);
end
keyboard


%**************************************************************************
% Example 3 
% Test Problem:     SRN
% Characteristics:  2 objectives, 2 constraints, 2 variables
%**************************************************************************

% Set inputs for the "SRN" test problem

fun   = @(x) fun_SRN(x);
con   = @(x) con_SRN(x);
xl    = [ -20,  -20];
xu    = [  20,   20];
feps  = [0.01, 0.01];
fu    = [1000,  100];
maxfn = 5000;

% Optimize with simDIRECT 

output = simDIRECT(fun,con,xl,xu,feps,fu,maxfn);

% Plot the Pareto set in input & output spaces

figure % output space
color = {'#0072BD','#D95319'};
xplot = output.f(~output.isPareto,1);
yplot = output.f(~output.isPareto,2);
plot(xplot,yplot,'+','Color',color{1},'MarkerSize',10,'LineWidth',1,...
    'DisplayName','dominated');
hold on;
xplot = output.f(output.isPareto,1);
yplot = output.f(output.isPareto,2);
plot(xplot,yplot,'o','Color',color{2},'MarkerSize',10,'LineWidth',1,...
    'DisplayName','nondominated');
axes = gca;
axes.XLim = [-100, 1000];
axes.XTick = -100:200:900;
axes.YLim = [-600, 200];
axes.YTick = -600:100:100;
axes.FontName = 'Times';
axes.FontSize = 14;
title('Pareto Set for SRN Problem:  Output Space',...
    'Interpreter','latex','FontSize',20);
xlabel('$f_1$','Interpreter','latex', 'FontSize', 20);
ylabel('$f_2$','Interpreter','latex', 'FontSize', 20, ...
    'Rotation',0,'HorizontalAlignment','right');
legend('Interpreter','latex');
axes.PlotBoxAspectRatio = [ 1,1,1 ];

figure % input space
xdom    = output.x(~output.isPareto,:);
xnondom = output.x(output.isPareto,:);
plot(xdom(:,1),xdom(:,2), '+','MarkerSize',10,'LineWidth',1,...
    'Color',color{1},'DisplayName','dominated');
hold on
plot(xnondom(:,1),xnondom(:,2), 'o','MarkerSize',10,'LineWidth',1,...
    'Color',color{2},'DisplayName','nondominated');    
axes = gca;
axes.FontName = 'Times';
axes.FontSize = 14;
title('Pareto Set for SRN Problem:  Input Space',...
    'Interpreter','latex','FontSize',20);
axes.XLim = [-20,20];
axes.YLim = [-20,20];
axes.YTick = -20:5:15;
xlabel('$x_1$','Interpreter','latex', 'FontSize', 20);
ylabel('$x_2$','Interpreter','latex', 'FontSize', 20,'Rotation',0,...
    'HorizontalAlignment','right');
legend('Interpreter','latex');
axes.PlotBoxAspectRatio = [ 1,1,1 ];
keyboard


%**************************************************************************
% Example 4 
% Test Problem:     DTLZ2 with nobj=3, nvar=8, x*=sqrt(0.5)
% Characteristics:  3 objectives, 0 constraints, 8 variables
%**************************************************************************

% Set inputs for the "DTLZ2" test problem
% with nobj=3, nvar=8, x*=sqrt(0.5)

fun   = @(x) DTLZ2(x,3,8,sqrt(0.5));
con   = [];
xl    = zeros(1,8);
xu    = ones(1,8);
feps  = [0.0001, 0.0001, 0.0001];
fu    = [1.5,  1.5,  1.5];
maxfn = 5000;

% Optimize with simDIRECT 

output = simDIRECT(fun,con,xl,xu,feps,fu,maxfn);

% Plot the Pareto set in the output space
    
figure
color = {'#0072BD','#D95319'};
gray = 0.5*ones(1,3);
R = 1; % earth radius
latspacing = 10; 
lonspacing = 10; 
% lines of longitude: 
[lon1,lat1] = meshgrid(-180:lonspacing:180,linspace(-90,90,300)); 
[x1,y1,z1] = sph2cart(lon1*pi/180,lat1*pi/180,R); 
plot3(x1,y1,z1,'-','color',gray)
hold on
% lines of latitude: 
[lat2,lon2] = meshgrid(-90:latspacing:90,linspace(-180,180,300)); 
[x2,y2,z2] = sph2cart(lon2*pi/180,lat2*pi/180,R); 
plot3(x2,y2,z2,'-','color',gray)
hold on;
f1 = output.f(output.isPareto,1);
f2 = output.f(output.isPareto,2);
f3 = output.f(output.isPareto,3);
h = scatter3(f1,f2,f3,'MarkerEdgeColor',color{2});
axes = gca;
axes.XLim = [0,1.5];
axes.YLim = [0,1.5];
axes.ZLim = [0,1.5];
axes.CameraPosition = [-6.8885   -8.7843   -3.6659];
title('Pareto Set for DTLZ2 Problem (3 obj, 8 var):  Output Space',...
    'Interpreter','latex','FontSize',18);
keyboard

%**************************************************************************
% Example 5 
% Test Problem:     LandH2x2 
%**************************************************************************

% Set inputs for the "LandH2x2" test problem

fun   = @(x) LandH2x2(x);
con   = [];
xl    = [-0.75, -2.50];
xu    = [ 0.75,  0.12];
feps  = [0.0001, 0.0001];
fu    = [-0.8, -0.8];
maxfn = 500;

% Optimize with simDIRECT 

output = simDIRECT(fun,con,xl,xu,feps,fu,maxfn);

% Plot Pareto set in output and input space

figure % output space
color = {'#0072BD','#D95319'};
xplot = output.f(~output.isPareto,1);
yplot = output.f(~output.isPareto,2);
plot(xplot,yplot,'+','Color',color{1},'MarkerSize',10,'LineWidth',1,...
    'DisplayName','dominated');
hold on;
xplot = output.f(output.isPareto,1);
yplot = output.f(output.isPareto,2);
plot(xplot,yplot,'o','Color',color{2},'MarkerSize',10,'LineWidth',1,...
    'DisplayName','nondominated');
axes = gca;
axes.XLim = [-2.2, -0.8];
axes.XTick = -2.2:0.2:-0.8;
axes.YLim = [-2.2, -0.8];
axes.YTick = -2.2:0.2:-0.8;
axes.FontName = 'Times';
axes.FontSize = 14;
title('Pareto Set for LandH2x2 Problem:  Output Space',...
    'Interpreter','latex','FontSize',20);
xlabel('$f_1$','Interpreter','latex', 'FontSize', 20);
ylabel('$f_2$','Interpreter','latex', 'FontSize', 20, ...
    'Rotation',0,'HorizontalAlignment','right');
legend('Interpreter','latex');
axes.PlotBoxAspectRatio = [ 1,1,1 ];

figure % input space
xdom    = output.x(~output.isPareto,:);
xnondom = output.x(output.isPareto,:);
plot(xdom(:,1),xdom(:,2), '+','MarkerSize',10,'LineWidth',1,...
    'Color',color{1},'DisplayName','dominated');
hold on
plot(xnondom(:,1),xnondom(:,2), 'o','MarkerSize',10,'LineWidth',1,...
    'Color',color{2},'DisplayName','nondominated');    
axes = gca;
axes.FontName = 'Times';
axes.FontSize = 14;
title('Pareto Set for SRN Problem:  Input Space',...
    'Interpreter','latex','FontSize',20);
axes.XLim = [xl(1),xu(1)];
axes.YLim = [xl(2),xu(2)];
xlabel('$x_1$','Interpreter','latex', 'FontSize', 20);
ylabel('$x_2$','Interpreter','latex', 'FontSize', 20,'Rotation',0,...
    'HorizontalAlignment','right');
legend('Interpreter','latex','FontSize',14);
axes.PlotBoxAspectRatio = [ 1,1,1 ];
title('Pareto Set for LandH2x2 Problem:  Input Space',...
    'Interpreter','latex','FontSize',20);
keyboard


%**************************************************************************
% Example 6 
% Test Problem:     DTLZ2 with nobj=2, nvar=10, x*=sqrt(0.5)
% Characteristics:  2 objectives, 0 constraints, 10 variables
%**************************************************************************

% Set inputs for the "DTLZ2" test problem
% with nobj=2, nvar=10, x*=sqrt(0.5)

fun   = @(x) DTLZ2(x,2,10,sqrt(0.5));
con   = [];
xl    = zeros(1,10);
xu    = ones(1,10);
feps  = 1e-4*ones(1,2);
fu    = [1.5,  1.5];
maxfn = 5000;

% Optimize with simDIRECT 

output = simDIRECT(fun,con,xl,xu,feps,fu,maxfn);

% Plot the Pareto set in the output space

figure
color = {'#0072BD','#D95319'};
axes = gca;
gray = 0.8*ones(1,3);
theta = linspace(0,pi/2,101);
xplot = cos(theta);
yplot = sin(theta);
plot(xplot,yplot,'-','Color',gray,'LineWidth',1,...
    'DisplayName','Pareto front');
hold on;
xplot = output.f(~output.isPareto,1);
yplot = output.f(~output.isPareto,2);
plot(xplot,yplot,'+','Color',color{1},'MarkerSize',6,'LineWidth',1,...
    'DisplayName','dominated');
xplot = output.f(output.isPareto,1);
yplot = output.f(output.isPareto,2);
plot(xplot,yplot,'o','Color',color{2},'MarkerSize',10,'LineWidth',1,...
    'DisplayName','nondominated');
axes.FontName = 'Times';
axes.FontSize = 12;
axes.XTick = 0:0.2:1;
axes.YTick = 0:0.2:1;
title('Pareto Set for DTLZ2 Problem (2 obj, 10 var):  Output Space',...
    'Interpreter','latex','FontSize',18);
mytitle.Interpreter = 'latex';
mytitle.FontSize = 18;
myxlabel = xlabel('$f_1$',...
                  'Interpreter','latex', 'FontSize', 18);
myylabel = ylabel('$f_2$',...
                  'Interpreter','latex', 'FontSize', 18, ...
                  'Rotation',0,...
                  'HorizontalAlignment','right');
axes.PlotBoxAspectRatio = [ 1,1,1 ];
axes.XLim = [0,1.6];
axes.YLim = [0,1.6];
axes.XTick = 0:0.5:1.5;
axes.YTick = 0:0.5:1.5;
legend('Interpreter','latex','FontSize',14);
keyboard

% ************************* Test functions ********************************

function f = Hartman6(x)
    % This test function is from the "DIRECTLib" collection of global
    % optimization test problems compiled by Linas Stripinis and 
    % Remigijus Paulavicius.   The full collection can be found at:
    % https://zenodo.org/record/5830927#.Yp4EUy1h1t9
    % Some simplifications have been made.
    if size(x, 2) > size(x, 1)
        x = x'; 
    end
    a(1,1) = 10.0; a(1,2) = 3.0;  a(1,3) = 17.0; a(1,4) = 3.5;  a(1,5) = 1.7;  a(1,6) = 8.0;
    a(2,1) = 0.05; a(2,2) = 10.0; a(2,3) = 17.0; a(2,4) = 0.1;  a(2,5) = 8.0;  a(2,6) = 14.0;
    a(3,1) = 3.0;  a(3,2) = 3.5;  a(3,3) = 1.7;  a(3,4) = 10.0; a(3,5) = 17.0; a(3,6) = 8.0;
    a(4,1) = 17.0; a(4,2) = 8.0;  a(4,3) = 0.05; a(4,4) = 10.0; a(4,5) = 0.1;  a(4,6) = 14.0;
    c(1) = 1.0;    c(2) = 1.2;    c(3) = 3.0;    c(4) = 3.2;
    p(1,1) = 0.1312; p(1,2) = 0.1696; p(1,3) = 0.5569; p(1,4) = 0.0124;	p(1,5) = 0.8283; p(1,6) = 0.5886;
    p(2,1) = 0.2329; p(2,2) = 0.4135; p(2,3) = 0.8307; p(2,4) = 0.3736;	p(2,5) = 0.1004; p(2,6) = 0.9991;
    p(3,1) = 0.2348; p(3,2) = 0.1451; p(3,3) = 0.3522; p(3,4) = 0.2883;	p(3,5) = 0.3047; p(3,6) = 0.6650;
    p(4,1) = 0.4047; p(4,2) = 0.8828; p(4,3) = 0.8732; p(4,4) = 0.5743;	p(4,5) = 0.1091; p(4,6) = 0.0381;
    s = 0;
    for i = 1:4
       sm=0;
       for j = 1:6
          sm = sm + a(i, j)*(x(j) - p(i, j))^2;
       end
       s = s + c(i)*exp(-sm);
    end
    f = -s;
end

function f = fun_NASA_speed_reducer(x)
    % This test function is from the "DIRECTLib" collection of global
    % optimization test problems compiled by Linas Stripinis and 
    % Remigijus Paulavicius.   The full collection can be found at:
    % https://zenodo.org/record/5830927#.Yp4EUy1h1t9
    % Some simplifications have been made.
    if size(x, 2) > size(x, 1)
        x = x'; 
    end
    f = 0.7854*x(1)*x(2)^2*(3.3333*x(3)^2 + 14.9334*x(3)-43.0934) -...
        1.508*x(1)*(x(6)^2 + x(7)^2) + 7.4777*(x(6)^3 + x(7)^3) +...
        0.7854*(x(4)*x(6)^2 + x(5)*x(7)^2);
end

function g = con_NASA_speed_reducer(x)
    % This test function is from the "DIRECTLib" collection of global
    % optimization test problems compiled by Linas Stripinis and 
    % Remigijus Paulavicius.   The full collection can be found at:
    % https://zenodo.org/record/5830927#.Yp4EUy1h1t9
    % Some simplifications have been made.
    g(1) = (27/(x(1)*x(2)^2*x(3))) - 1;
    g(2) = (397.5/(x(1)*x(2)^2*x(3)^2)) - 1;  
    g(3) = ((1.93*x(4)^3)/(x(2)*x(3)*x(6)^4)) - 1;  
    g(4) = ((1.93*x(5)^3)/(x(2)*x(3)*x(7)^4)) - 1;  
    g(5) = ((sqrt((((745*x(4))/(x(2)*x(3)))^2) + 16.9 * 10^6)) * (1/(110 * x(6)^3))) - 1; 
    g(6) = (((((745*x(5))/(x(2)*x(3)))^2 + 157.5 * 10^6)^(1/2))/(85 * x(7)^3)) - 1 ; 
    g(7) = ((x(2)*x(3))/(40)) - 1;   
    g(8) = ((5*x(2))/(x(1))) - 1; 
    g(9) = ((x(1))/(12*x(2))) - 1; 
    g(10) = (((1.5*x(6)+1.9))/(x(4))) - 1; 
    g(11) = ((1.1*x(7) + 1.9)/(x(5))) - 1; 
end

function f = fun_SRN(x)
    % Source:  K. Deb, P. Roy, R. Hussein, "Surrogate Modeling Approaches
    % for Multiobjective Optimization:  Methods, Taxonomy, and Results,"
    % Mathematical and Computational Applications," vol 26, no. 5 (April
    % 2021), p.5.   https://dx.doi.org/10.3390/mca26010005
    x1 = x(1);
    x2 = x(2);
    f1 = 2 + (x1-2)^2 + (x2-1)^2;
    f2 = 9*x1 - (x2-1)^2;
    f = [f1,f2];
end

function g = con_SRN(x)
    % Source:  K. Deb, P. Roy, R. Hussein, "Surrogate Modeling Approaches
    % for Multiobjective Optimization:  Methods, Taxonomy, and Results,"
    % Mathematical and Computational Applications," vol 26, no. 5 (April
    % 2021), p.5.   https://dx.doi.org/10.3390/mca26010005
    x1 = x(1);
    x2 = x(2);
    g1 = x1^2 + x2^2 - 225;
    g2 = x1 - 3*x2 + 10;
    g = [g1,g2];
end

function f = DTLZ2(x,nobj,nvar,xstar)
    % Source: https://github.com/fcampelo/DEMO/blob/master/Octave-Matlab/DTLZ
    % Some modifications and simplifications were done to the original.
    % Set number of objectives and k (which determines nvar=M+k-1)
    M = nobj;
    k = nvar - M + 1;
    % transpose x because it is a row vector and Deb expects column vector
    x = x';
    % compute the objectives
    n = (M-1) + k; 
    xm = x(n-k+1:end,:); %xm contains the last k variables
    g = sum((xm - xstar).^2, 1);
    f(1,:) = (1 + g).*prod(cos(pi/2*x(1:M-1,:)),1);
    for ii = 2:M-1
       f(ii,:) = (1 + g) .* prod(cos(pi/2*x(1:M-ii,:)),1) .* ...
          sin(pi/2*x(M-ii+1,:));
    end
    f(M,:) = (1 + g).*sin(pi/2*x(1,:));
    f = f';
end

function f = LandH2x2(x)
    a = -sqrt(2)/2;
    b = -sqrt(4*pi/65);
    c = -sqrt((90*pi)/112);
    d = .65;
    e = 2.8;
    x1 = x(1);
    x2 = x(2);
    f1 =  a*x1 + b*exp( - (x1^2 + x2^2)/d^2 ) + c * exp( - ( x1^2 + (x2+1.5)^2 )/e^2 );
    f2 = -a*x1 + b*exp( - (x1^2 + x2^2)/d^2 ) + c * exp( - ( x1^2 + (x2+1.5)^2 )/e^2 );
    f = [f1,f2];    
end
