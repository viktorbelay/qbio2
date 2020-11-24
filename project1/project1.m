% Qbio 2 Project 1
rng('default')
s = rng;
%% Question 1
u = random('norm',5,1,1000,1);
v = random('norm',5,1,1000,1);
x = 0.25 * u + 0.75 * v;
y = 25 * u - 75 * v + 700;
datam = [x y];
% Scatter plot of the data
figure(1)
plot(x,y,'b.')
xlabel('x data')
ylabel('y data') 
% Find covariance matrix and manually find eigensystem
A = cov(x,y);
[E,A1] = eigs(A);

U = 1./A1;
U1 = [U(1) 0; 0 U(2,2)];
U2 = sqrt(U1);
W = E*U2;
% Perform PCA and plot vectors on scatter plot
[pcs, trans, evs] = pca(datam);
figure(2)
plot(x,y,'.')
xlabel('x data')
ylabel('y data') 
v1 = E(:,1);
hold on 
plot([-v1(1), v1(1)],[-v1(2), v1(2)], 'k', 'linewidth', 3)
% Mean centering and scaling data based on standard deviation
mux = mean(x);
muy = mean(y);
sdx = std(x);
sdy = std(y);
xc = x - mux;
yc = y - muy;
xcs = xc/sdx;
ycs = yc/sdy;
datam1 = [xcs,ycs];
figure (3)
plot(xcs,ycs,'.')
xlabel('x data')
ylabel('y data') 
% Perform PCA on new dataset
A22 = cov(xcs,ycs);
[E1,A2] = eigs(A22);
[pcs1,trans1,evs1] = pca(datam1);

V = 1./A2;
V1 = [V(1) 0; 0 V(2,2)];
V2 = sqrt(V1);
W2 = E*V2;

if eq(mean(x),5.0195) % Test for randomness
    "rnd is default, data is not random"
else
    "Warning: data is random"
end


%% Question 2
tic % Start recording the time it takes for this simulation to run
% Simulation parameters

tspan = [0 200000];
y0 = [0;0;0;1;0;0];
[t,y] = ode45(@de,tspan,y0);


% Plotting
figure(4)

subplot(6,1,1)
plot(t,y(:,1),'k','linewidth',2)
ylabel('Concentration of mRNA A')

subplot(6,1,2)
plot(t,y(:,4),'k','linewidth',2)
ylabel('Concentration of protein A')

subplot(6,1,3)
plot(t,y(:,2),'k','linewidth',2)
ylabel('Concentration of mRNA B')

subplot(6,1,4)
plot(t,y(:,5),'k','linewidth',2)
ylabel('Concentration of protein B')

subplot(6,1,5)
plot(t,y(:,3),'k','linewidth',2)
ylabel('Concentration of mRNA C')

subplot(6,1,6)
plot(t,y(:,6),'k','linewidth',2)
ylabel('Concentration of protein C')
xlabel('time')

toc % Output the time it took to run this simulation
function dydt = de(t,y)

% ODE parameters

a0 = 4;
a = 1.8 ;
b = 0.01;
tau = 0.005;

% DEs
dydt = [
    a0-y(1)-a*y(6); 
    a0-y(2)-a*y(4);
    a0-y(3)-a*y(5); 
    b*(y(1)-y(4)); 
    b*(y(2)-y(5)); 
    b*(y(3)-y(6))
    ];

end

% Perform PCA and plot vectors on scatter plot
[pcs, trans, evs] = pca(datam);
figure(2)
plot(x,y,'.')
xlabel('x data')
ylabel('y data') 
v1 = E(:,1);
hold on 
plot([-v1(1), v1(1)],[-v1(2), v1(2)], 'k', 'linewidth', 3)
% Mean centering and scaling data based on standard deviation
mux = mean(x);
muy = mean(y);
sdx = std(x);
sdy = std(y);
xc = x - mux;
yc = y - muy;
xcs = xc/sdx;
ycs = yc/sdy;
datam1 = [xcs,ycs];
figure (3)
plot(xcs,ycs,'.')
xlabel('x data')
ylabel('y data') 
% Perform PCA on new dataset
A22 = cov(xcs,ycs);
[E1,A2] = eigs(A22);
[pcs1,trans1,evs1] = pca(datam1);

V = 1./A2;
V1 = [V(1) 0; 0 V(2,2)];
V2 = sqrt(V1);
W2 = E*V2;

if eq(mean(x),5.0195) % Test for randomness
    "rnd is default, data is not random"
else
    "Warning: data is random"
end


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





