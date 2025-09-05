
%Definition of the test function and its derivative
test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
test_derivative01 = @(x) 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;

x = linspace(-15,40,1000);

fx = test_func01(x);

plot(x, fx)
grid on;
yline(0, '--k', 'LineWidth', 1);

% Function calls

bisection()
newton()
secant()
bisection_solver(test_func01,20,40);
newton_solver(@my_test_funct, 30);
secant_solver(test_func01, 30,40);

function [fout, dfdxout] = my_test_funct(x)
  fout = (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
  dfdxout = 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;
end

function bisection
    test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
    interval = [20,40];

    c = (interval(1,1)+interval(1,2))/2;
    
    fc = test_func01(c);

    while abs(fc) > 0.000001
        c = (interval(1,1)+interval(1,2))/2;
        fa = test_func01(interval(1,1));
        fb = test_func01(interval(1,2));
        fc = test_func01(c);

        if fa * fc < 0
            % f(a) and f(c) have opposite signs
            interval = [interval(1,1),c];
        else
            % f(a) and f(c) have same signs
            if fb * fc < 0
                % f(b) and f(c) have opposite signs
                interval = [c,interval(1,2)];
            else
                % f(b) and f(c) have same signs
            end
        end
    end
        x = interval(1,1);
end

function newton
    test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
    test_derivative01 = @(x) 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;
    x=30;
        while abs(test_func01(x))>=0.000000001
            x=x-test_func01(x)/test_derivative01(x);
        end
end

function secant
    test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
    x0=30;
    x1=40;
    while abs(test_func01(x1))>=0.00000001
        fx0=test_func01(x0);
        fx1=test_func01(x1);
        x2=(x0*fx1-x1*fx0)/(fx1-fx0);
        x0=x1;
        x1=x2;
    end
    x = x2;
end

%generalized bisection function
function x = bisection_solver(func,x_left,x_right)
    
    interval = [x_left,x_right];

    c = (interval(1,1)+interval(1,2))/2;
    
    fc = func(c);

    while abs(fc) > 1e-7
        c = (interval(1,1)+interval(1,2))/2;
        fa = func(interval(1,1));
        fb = func(interval(1,2));
        fc = func(c);

        if fa * fc < 0
            % f(a) and f(c) have opposite signs
            interval = [interval(1,1),c];
        else
            % f(a) and f(c) have same signs
            if fb * fc < 0
                % f(b) and f(c) have opposite signs
                interval = [c,interval(1,2)];
            else
                % f(b) and f(c) have same signs
            end
        end
    end
       x = interval(1,1);
end

%generalized newton function
function x = newton_solver(f_in,x0)

    x_val = x0;
    [f_val, dfdx_val] = f_in(x_val);
    x_prev=x0+1;

    while abs(dfdx_val) > 1e-14 && abs(x_val-x_prev)>0
        x_prev = x_val;
        x_val = x_val - f_val/dfdx_val;
        [f_val, dfdx_val] = f_in(x_val);
    end
    x = x_val;
end

%generalized secant function
function x = secant_solver(fun,x0,x1)
    test_func01 = fun;
    while abs(test_func01(x1)) >= 1e-14
        fx0 = test_func01(x0);
        fx1 = test_func01(x1);
        x2 = (x0 * fx1 - x1 * fx0) / (fx1 - fx0);
        x0 = x1;
        x1 = x2;
    end
    x = x2;
end