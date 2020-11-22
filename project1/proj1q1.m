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
