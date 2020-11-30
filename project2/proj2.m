% Project 2 q1 
 
% Plot nulclines and quiver to generate phase portrait

a = 1;
n = 3;

x1 = 0:0.1:10;
y1 = a ./ (1 + x1 .^ n);
y2 = 0:0.1:10;
x2 = a ./ (1 + y2 .^ n);
figure(1)
plot(x1,y1,x2,y2,'linewidth',2);

x3 = fsolve(@switch_ss, [0.7;0.7]); % Solve for critical point

% Find Jacobian and eigensystem of Jacobian for equilibrium and stability analysis

syms x y

global a 
global n
a = 1;
n = 3;

J = [D

jacobian([-x+(a/(1+(y)^n)), -y+(a/(1+(x)^n))],[x,y])

x = 0;
y = 0;

% Simulation parameters
tspan = [0 100];
z0 = [0;0];
[t,z] = ode45(@de,tspan,z0);


% Plotting

figure(2)
subplot(2,1,1)
plot(t,z(:,1))
subplot(2,1,2)
plot(t,z(:,2))

figure(3)

plot(z(:,1),z(:,2))

function dzdt = de(t,z)

% ODE parameters
a = 1;
n = 3;

% DEs
dzdt = [-z(1)+(a/(1+(z(2))^n));
    -z(2)+(a/(1+(z(1))^n))
    ];

end

function [ F ] = switch_ss( v )
%SWITCH_SS look nullclines for switch model
x = v(1);
y = v(2);
a = 1;
n = 3;

F = [ x - a/(1+y^n);
 y - a/(1+x^n)];

end
