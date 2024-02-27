function [R_new] = T(R_old,f_pdf,b,phi,beta)

% Integrand function
myfun = @(w) (w-R_old).*f_pdf(w);

% Numerical integral
approx_integ = integral(myfun,R_old,inf);

R_new = b*(1-phi*beta)/(1-beta)+(beta/(1-beta))*approx_integ;


end