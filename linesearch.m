function [lambda, No_of_iterations] = linesearch(func,x,d)

% Armijo's rule, init
epsilon = 10^-2;
alpha = 6;
No_of_iterations = 0;
lambda_curr = 0;
dir = d;

% Derivation of function
F = @(lambda) func(x + lambda * dir);
h = 10^-8;

% Find a proper lambda to start with
lambda_add = 10^-3;
while isinf(F(lambda_curr + h)) || isinf(F(lambda_curr - h)) || ...
    abs(F(lambda_curr + h) - F(lambda_curr - h)) > 10^2
    dir = dir./7;
    F = @(lambda) func(x + lambda * dir);
end

while isinf(F(lambda_curr + h)) || isinf(F(lambda_curr - h)) || ...
    abs(F(lambda_curr + h) - F(lambda_curr - h)) < 10^-12
    dir = dir*7;
    F = @(lambda) func(x + lambda * dir);
end

% Iterate until close enough
while norm(lambda_add,2) > 10^-7 && No_of_iterations < 100;
    % Derivative
    F_prim = (F(lambda_curr + h) - F(lambda_curr - h))/(2*h);
    T = @(lambda) F(lambda_curr) + epsilon * lambda * F_prim;
    if F_prim < 0
        lambda_add = 10^-3;
    else
        lambda_add = -10^-3;
    end
    
    % Get inside interval
    while F(lambda_curr + lambda_add) > T(lambda_add)
        lambda_add = lambda_add/alpha;
    end
    while F(alpha*(lambda_curr + lambda_add)) < T(alpha*lambda_add)
        lambda_add = lambda_add * alpha;
    end
    
    lambda_curr = lambda_curr + lambda_add;
    No_of_iterations = No_of_iterations + 1;
end

for k = 1:length(d)
    if d(k) ~= 0
        lambda = lambda_curr * dir(k)/d(k);
        break;
    end
end