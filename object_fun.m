function f = object_fun(x)
% objective function
% consist of raw object function and penalty functions

% Constraints
c = zeros(1,2);
% Penalty factors
penalty = zeros(1,2);

of_raw = 10*(x(1)-1)^2+20*(x(2)-2)^2+(x(3)-3)^2;

% Define constraints
% C1: x1+x2+x3 <= 5
% C2: x1^2+2*x2 < x3
c(1) = x(1)+x(2)+x(3)-5;
c(2) = x(1)^2+2*x(2)-x(3);

% Define penalty factor for each constraint
penalty(c>0) = 1;

% Penalty factor
pf = 1000;

f = of_raw + pf*sum(penalty);
end

