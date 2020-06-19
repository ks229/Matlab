%using two parameter gama prior and posterior
clc;
%clear all;
%Reliabilty study using Bayes theorem to estimate system parameters
lamda_actual = 0.1;%failure rate
theta_actual = 0.01;
mu_actual = 2;
beta_actual = 0.05;
c = 0.95;%coverage factor
n = 1000;
w1 = 1;w2 = 0.1;w3 = 20;w4 = 0.5;
v1 = 1;v2 = 1;v3 = 1;v4 = 1;
%data needed given above
U1 = exprnd(1/lamda_actual,n,1);
U2 = exprnd(1/theta_actual,n,1);
U3 = exprnd(1/mu_actual,n,1);
U4 = exprnd(1/beta_actual,n,1);
T1= sum(U1);
T2 = sum(U2);
T3 = sum(U3);
T4 = sum(U4);
mcrun = 10000;
%Prior distribution
mean_lamda = w1/v1;
var_lamda = w1/v1^2;
lamda_o = linspace(0,50,n);
lamda_prior = gampdf(lamda_o,w1,v1);
%posterior distribution
pos = zeros(4,mcrun);
mean_lamda_pos = (n+w1)/(T1 + v1);
mean_theta_pos = (n+w2)/(T2 + v2);
mean_mu_pos = (n+w3)/(T3 + v3);
mean_beta_pos = (n+w4)/(T4 + v4);
%plot(lamda_prior);
lamda = Gamma(mcrun,(T1+v1),(n+w1));
theta = Gamma(mcrun,(T2+v2),(n+w2));
mu = Gamma(mcrun,(T3+v3),(n+w3));
beta = Gamma(mcrun,(T4+v4),(n+w4));
for i = 1:mcrun
    A(i) = (mu(i)*beta(i)*(2*lamda(i)+2*lamda(i)*theta(i)+theta(i)*mu(i)+c))/(mu(i)*beta(i)*(2*lamda(i)+2*lamda(i)*theta(i)+theta(i)*mu(i)+c)+2*lamda(i)*theta(i)*(lamda(i)*beta(i)+(1-c)*mu(i)^2));
    MTTF(i) = ((2*c*lamda(i) + theta(i))*(lamda(i) + mu(i))+ 2*c*lamda(i)*theta(i))/(2*lamda(i)*theta(i)*(lamda(i) + (1-c)*mu(i)));
end
histfit(A,n);
%hpdi_A = hpdi(A',4);
