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
llim =-15*pi; % Left/Lower Domain Limit
rlim = 15*pi; % Right/Upper Domain Limit
N = 1001; % Number of node points
y = llim:(rlim-llim)/(N-1):rlim; % Node Locations
v = zeros(N,1); % Initialize velocity vector
cr = 0.5; % Real part of Phase Speed

options = optimoptions('fsolve','Display','off','TolFun',1e-9); % Fsolve options

dalpha = 0.1; % Wave number stepping for plotting alpha vs. omega
alphaplot = size(1/dalpha+1,1); % Variable to store alpha values for plotting
omegaplot = size(1/dalpha+1,1); % Variable to store omega values for plotting
j = 1; % iterator to prevent dynamic allocation
fprintf('alpha\treal(c)\timag(c)\treal(omega)\timag(omega)\n');
fprintf('----------------------------------------------------\n');

for alpha=1:-dalpha:0    % Wave number in x-direction
    ci = 0.03; % Initial guess for imaginary part of phase speed
    ci = fsolve(@(ci)fun1(alpha,cr,llim,rlim,ci),ci,options); % Find ci that satisfies boundary conditions
    ci = real(ci); % Take the real part in case it returns imaginary number
    
    % The rest is plotting and storing
    omegar = alpha*cr; % omega = alpha * c
    omegai = alpha*ci;
    alphaplot(j) = alpha; % Storing alpha for plotting
    omegaplot(j) = omegai; % Storing omega_i for plotting
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
function F = fun1(alpha,cr,llim,rlim,ci)
    c = cr + 1i*ci; % Phase speed
    [~,v] = ode45(@(y,v)vdp1(y,v,alpha,c),[llim rlim],[0; 1]); %RK45 for solution
    F = real(v(end))-0; % v(N)-rBC = 0 function
end