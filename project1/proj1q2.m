% Project 1

%% Question 2

% Simulation parameters

tspan = [0 100];
y0 = [0;0;0;1;1;1];
[t,y] = ode45(@de,tspan,y0);


% Plotting
figure(4)

subplot(6,1,1)
plot(t,y(:,1))
ylabel('Concentration of mRNA A')

subplot(6,1,2)
plot(t,y(:,4))
ylabel('Concentration of protein A')

subplot(6,1,3)
plot(t,y(:,2))
ylabel('Concentration of mRNA B')

subplot(6,1,4)
plot(t,y(:,5))
ylabel('Concentration of protein B')

subplot(6,1,5)
plot(t,y(:,3))
ylabel('Concentration of mRNA C')

subplot(6,1,6)
plot(t,y(:,6))
ylabel('Concentration of protein C')
xlabel('time')

function dydt = de(t,y)

% ODE parameters

a0 = 4;
a = 2.5 ;
b = 0.01;
tau = 0.005;

% DEs
dydt = [
    tau*(a0-y(1)-a*y(6))+y(1); 
    tau*(a0-y(2)-a*y(4))+y(2);
    tau*(a0-y(3)-a*y(5))+y(3); 
    (tau*b*(y(1)-y(4)))+y(4); 
    b*(y(2)-y(5)) + y(5); 
    (tau*b*(y(3)-y(6))) + y(6)
    ];

end
