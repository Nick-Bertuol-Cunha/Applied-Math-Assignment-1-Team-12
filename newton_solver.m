function [x, xi_values] = newton_solver(fun, x0)
    tol = 1e-9; max_iter = 100;
    x = x0; xi_values = [];
    for n = 1:max_iter
        [fval, dfdx] = fun(x);          
        if ~isfinite(fval) || abs(fval) < tol, break; end
        if ~isfinite(dfdx) || dfdx == 0,  break; end
        x = x - fval/dfdx;
        xi_values(end+1) = x;      
    end
end