function [fval,dfdx] = orion_test_func2(x)
    fval = x.^2-2;
    dfdx = 2*x;
end