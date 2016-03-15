function x = gaussnewton(phi,t,y,start,tol,use_linesearch,printout,plotout)

% Initiate
iteration = 1;
size_x = length(start);
x = start;
h = 1e-8;
lambda = 1;
ls_iter = 0;
d = 1;

if printout == 1
    fprintf('%s \t%s \t\t%s \t%s \t\t%s \t%s \t%s \t%s \t%s\n', ...
    'iter', 'x', 'step size', '  f(x)', 'max(abs(r))', 'norm(grad)', ... 
    'ls iter', '  lambda', 'grad`*d/norm(d)');
end

while norm(lambda*d) > tol && iteration < 100 % Stopping criterion
    
    % Calculate Jacobian
    J = zeros(length(t), size_x);
    I = eye(size_x);
    for i = 1:size_x
        theJ = (phi(x+h*I(:,i),t) - phi(x-h*I(:,i),t))/(2*h);
        J(:,i) = theJ(:,i);
        %J(:,i) = (phi(x+h*I(:,i),t) - phi(x-h*I(:,i),t))/(2*h);
    end
    thephi=phi(x,t');
    r = (thephi(:,1)'-y);
    %r = (phi(x,t) - y);
    JJ = J'*J;
    if isnan(norm(JJ)) || isinf(norm(JJ))
        error('Bad starting point, choose another');
    end
    
    % If J'J is close to singular, add some epsilon
    epsilon = 10^-6*norm(JJ);
    while rcond(JJ) < 1e-6 || isnan(rcond(JJ))
        JJ = JJ + eye(size_x)*epsilon;
    end
    d = (JJ)\(-J'*r');
    % Check that d is ok
    for m = 1:size_x
        if isnan(d(m)) || isinf(d(m))
            d = zeros(size_x,1);
            d(mod(iteration,size_x)+1)=1;
        end
    end
    
    f = @(alpha) sum((phi(alpha,t) - y).^2);     % Residuals
    
    if use_linesearch == 1
        [lambda, ls_iter] = linesearch(f, x, d);
    end
    x = x + lambda*d;
    
    grad_f = zeros(size_x, 1);
    for j = 1:size_x
        grad_f(j,1) = (f(x+I(:,j)*h) - f(x - I(:,j)*h))/(2*h);
    end
    
    if printout  == 1
        fprintf('%3d %8.4f \t%8.4f \t%7.4f \t%8.4f \t\t%8.4f \t%5d \t\t%8.4f \t%8.4f\n', iteration, x(1), norm(lambda*d), ...
        f(x), max(abs(r)), norm(grad_f), ls_iter, lambda, grad_f'*d/norm(d));
        for k = 2:size_x
            fprintf('\t%8.4f\n', x(k));
        end
    end
    iteration = iteration + 1;
end

if plotout == 1
    clf
    figure(1)
    plot(t,y,'x');
    theY=phi(x,t');
    hold on;
    %plot(t, phi(x, t));
    plot(t,theY(:,1));
end
end