%%%%%%%%%%%%%%%%%%%%%%% Rayleigh's Equation Solver %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This code investigates the hydrodynamic stability of a parallel,  %%%%
%%%% incompressible and inviscid shear flow. The equation is:          %%%%
%%%%                                                                   %%%%
%%%% (U-c)(v''-alpha^2*v)-U''*v=0                                      %%%%
%%%%                                                                   %%%%
%%%% The boundary conditions:                                          %%%%
%%%%                                                                   %%%%
%%%% v--> -Inf, v=0 / v--> +Inf, v=0                                   %%%%
%%%%                                                                   %%%%
%%%% Furkan Oz                                                         %%%%
%%%% Oklahoma State University                                         %%%%
%%%% foz@okstate.edu                                                   %%%%
%%%% 2/10/2022                                                         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; format long;

% Initialization
llim =-30*pi; % Left/Lower Domain Limit
rlim = 30*pi; % Right/Upper Domain Limit

options = optimoptions('fsolve','Display','off','StepTolerance',1e-9,'FunctionTolerance',1e-6,'MaxFunctionEvaluations',2000,'Algorithm','levenberg-marquardt'); % Fsolve options

dalpha = 0.1; % Wave number steping for plotting alpha vs. omega
alphaplot = size(1/dalpha+1,1); % Variable to store alpha values for plotting
omegaplot = size(1/dalpha+1,1); % Variable to store omega values for plotting
cplot = size(1/dalpha+1,1); % Variable to store omega values for plotting
j = 1; % iterator to prevent dynamic allocation
fprintf('alpha\treal(c)\timag(c)\treal(omega)\timag(omega)\n');
fprintf('----------------------------------------------------\n');
c(1) = 0.5; % Initial guess for real part of the phase speed
c(2) = 0.03; % Initial guess for imaginary part of the phase speed
for alpha=0.9:-dalpha:0.00    % Wave number in x-direction
    c = fsolve(@(c)fun1(alpha,llim,rlim,c),c,options); % Find ci that satisfies boundary conditions
    cr = c(1); % Take the real part
    ci = c(2); % Take the imaginary part
    
    % The rest is plotting and storing
    omegar = alpha*cr; % omega = alpha * c
    omegai = alpha*ci;
    alphaplot(j) = alpha; % Storing alpha for plotting
    omegaplot(j) = omegai; % Storing omega_i for plotting
    cplot(j) = ci;
    j = j+1; % iterator advancement
    fprintf('%.2f\t%.2f\t%.6f\t%.2f\t%.6f\n',alpha,cr,ci,omegar,omegai);
end
plot(alphaplot,omegaplot,'linewidth',2);
ax = gca; ax.FontSize = 11; 
xlabel('$$\alpha$$','FontSize',18,'interpreter','latex');
ylabel('$$\omega_i$$','FontSize',18,'interpreter','latex');
grid on
ylh = get(ax,'ylabel'); gyl = get(ylh); ylp = get(ylh, 'Position');
set(ylh,'Rotation',0,'Position',ylp,'HorizontalAlignment','right')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RK45 Function
function dydx = vdp1(y,v,alpha,c)
    U = 0.5*(1+tanh(y)); % Background Velocity
    Uder2 =-sech(y).^2.*tanh(y); % 2nd Derivative of Background Velocity
    dydx = [v(2)
       -(alpha^2*c-alpha^2*U-Uder2)*v(1)/(U-c)];
end

% fsolve function to satisfy
function F = fun1(alpha,llim,rlim,c)
    c = c(1) + 1i*c(2); % Phase speed
    [~,v] = ode45(@(y,v)vdp1(y,v,alpha,c),[llim rlim],[0; 1]); %RK45 for solution
    F = [(real(v(end,1))-0); imag(v(end,1))-0];% v(N)-rBC = 0 function
end
