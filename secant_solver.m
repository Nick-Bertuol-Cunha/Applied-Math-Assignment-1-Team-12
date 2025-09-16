function x = secant_solver(fun,x0,x1)
    test_func01 = fun;
    x1_values = [];
    while abs(test_func01(x1)) >= 0.00000001
        fx0 = test_func01(x0);
        fx1 = test_func01(x1);
        x2 = (x0 * fx1 - x1 * fx0) / (fx1 - fx0);
        x1_values = [x1_values, x1];
        x0 = x1;
        x1 = x2;
    end
    x = x1; 
    disp('Iterated x1 values:');
    disp(x1_values);
end