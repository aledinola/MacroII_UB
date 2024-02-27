%% Problem set 1, Question 2
% Macroeconomics II, Master in Econ
% Author: Alessandro Di Nola
clear
clc 
close all

%% Economic parameters 

b      = 0.5;
beta   = 0.9;
phi    = 1.03;
lambda = 0.2;

% pdf of the exponential distribution
f_pdf = @(w) lambda*exp(-lambda*w);

%% Set numerical parameters
tol  = 1e-6; % Tolerance criterion
damp = 0.3;  % Dampening parameter, called "psi" in the problem set

%% Fixed point iteration: R=T(R)

% Initial guess
R_old = 1; % More or less arbitrary (as long as R_old>0!!)

% Test the function T()
R_new = T(R_old,f_pdf,b,phi,beta);

% If R_new is close to R_old, I have found a solution
% If not, keep on iterating. We automatize this using a "while" loop

dist = tol+1;
iter = 1;
R_store = []; % To store result of each iteration

while dist>=tol && iter<=10000

    R_new = T(R_old,f_pdf,b,phi,beta);

    % Store here iteration results, by appending new values
    R_store = [R_store;R_new];
    % Check distance between R new and R old
    dist = abs(R_new-R_old);
    
    fprintf('Iter = %d, dist = %f \n',iter,dist)

    % Update
    R_old = damp*R_new+(1-damp)*R_old;
    iter = iter+1;

end

disp('Reservation wage:')
disp(R_old)

figure
plot(1:iter-1,R_store,'LineWidth',2)
xlabel('Iteration n')
ylabel('Reservation wage at iter n')
yline(R_old)








