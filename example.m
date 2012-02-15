% Test the "lbfgs" Matlab interface on the Hock & Schittkowski test problem
% #38. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
% Nonlinear Programming Codes. Lecture Notes in Economics and Mathematical
% Systems Vol. 187, Springer-Verlag.

% The starting point.
x0  = [-3  -1  -3  -1];   % The starting point.
lb  = [-10 -10 -10 -10];  % Lower bound on the variables.
ub  = [+10 +10 +10 +10];  % Upper bound on the variables.

objHS038 = @(x) (100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 +  ...
			(1-x(3))^2 + 10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1));
gradHS038 = @(x) [	(-400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1))) ...
					(200*(x(2)-x(1)^2) + 20.2*(x(2)-1) + 19.8*(x(4)-1)) ...
					(-360*x(3)*(x(4)-x(3)^2) -2*(1-x(3))) ...
					(180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1))];

obj = @(x) deal(objHS038(x), gradHS038(x));
x = lbfgsb(obj, x0, lb, ub)