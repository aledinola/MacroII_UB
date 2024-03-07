%% There are two ways to solve this problem. 
% METHOD 1: Work with the value U
% The relevant equation is
% U = b+beta*(1-lambda)*U+beta*lambda*E[max(V(w),U)], where V(w)=w/(1-beta)
% METHOD 2: Work using the two wage values given in the question (less
% general)
% The relevant equation is 
% R=b+beta*lambda/(1-beta)*[0.5*max(wL-R,0)+0.5*max(wH-R,0)]
clear
clc
close all

%% Parameters
beta   = 0.75;     % Discount factor
lambda = 1.0;      % Prob of receiving a job offer: either 1 or 0.5
b      = 2;        % Unemployment benefits
wL     = 2.7;      % Low-wage offer
wH     = 3.3;      % High-wage offer
pr     = 0.5;      % Distribution of wage offers
w_grid = [wL,wH]'; % Support of wage distrib
w_prob = [pr,pr]'; % Wage distrib (probability mass function)

%% Initial guess
R = 2;
weight_old = 0.9;

%% Fixed point iteration using eq.2
tol = 1e-8;
err = tol+1;
iter = 1;

while err>tol

    R_new = fun_rhs(R,pr,wL,wH,b,beta,lambda);

    err = abs(R-R_new);

    fprintf('iter = %d, err = %f \n',iter,err)

    iter = iter+1;
    R = (1-weight_old)*R_new+weight_old*R;

end

%% Initial guess for iteration over U
U = 2;

%% Fixed point iteration using eq.1
tol = 1e-8;
err = tol+1;
iter = 1;
weight_old = 0.5;

while err>tol

    U_new = fun_rhs2(U,b,beta,lambda,w_grid,w_prob);

    err = abs(U-U_new);

    fprintf('iter = %d, err = %f \n',iter,err)

    iter = iter+1;
    U = (1-weight_old)*U_new+weight_old*U;

end

R_anal = fun_anal(pr,wL,wH,b,beta,lambda);

disp('Numerical R, method 1:')
disp(R)

disp('Numerical R, method 2:')
disp(U*(1-beta))

disp('Analytical R:')
disp(R_anal)

disp('Low wage   High  wage   Reservation wage')
disp([wL,wH,R_anal])

% Probability of finding a job
% Compute 0-1 indicator for whether w'>=res_wage, for each w' in the
% support
accept = w_grid>=R_anal;
jfp    = lambda*sum(w_prob(accept));

disp('Prob. of finding a job')
disp(jfp)

function R_new = fun_rhs(R,pr,wL,wH,b,beta,lambda)

opt_val = pr*max(wL-R,0)+pr*max(wH-R,0);
R_new = b+(beta*lambda/(1-beta))*opt_val;

end

function U_new = fun_rhs2(U,b,beta,lambda,w_grid,w_prob)

val_accept = w_grid/(1-beta);
U_new = b+beta*(1-lambda)*U+beta*lambda*sum(max(val_accept,U).*w_prob);

end

function R_anal = fun_anal(pr,wL,wH,b,beta,lambda)
% We always assume that R<=wH

% CASE 1: Assume wL<R
R_1 = (2*(1-beta)*b+beta*lambda*wH)/(2*(1-beta)+beta*lambda);

% CASE 2: Assume R<wL 
E_wage = pr*wL+(1-pr)*wH;
R_2 = (b*(1-beta)+beta*lambda*E_wage)/(1-beta+beta*lambda);

if wL<R_1
    R_anal = R_1;
elseif R_1<=wL
    R_anal = R_2;
end

end %end function fun_anal
