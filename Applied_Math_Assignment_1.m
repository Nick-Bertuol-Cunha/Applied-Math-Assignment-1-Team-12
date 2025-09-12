my_recorder = input_recorder();

%Day 1
test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
test_derivative01 = @(x) 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;

x = linspace(-15,40,1000);
fx = test_func01(x);
%plot(x, fx)
%grid on;

%bisection()
%newton()
%secant()
%bisection_solver(@(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6),20,40)
%newton_solver(@orion_test_func,30)
%secant_solver(@orion_test_func, 1,3)
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
        
end

function newton()
    test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
    test_derivative01 = @(x) 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;
    x=30;
        while abs(test_func01(x))>=0.000000001
            x=x-test_func01(x)/test_derivative01(x);
        end 
end

function secant()
    test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
    x0=1;
    x1=3;
    while abs(test_func01(x1))>=0.00000001
        fx0=test_func01(x0);
        fx1=test_func01(x1);
        x2=(x0*fx1-x1*fx0)/(fx1-fx0);
        x0=x1;
        x1=x2;
    end
end
%% 

%generalized bisection function
function x = bisection_solver(fun,x_left,x_right)
    func = fun;
    interval = [x_left,x_right];
    c_values = [];

    c = (interval(1,1)+interval(1,2))/2;
    
    fc = func(c);

    while abs(fc) > 0.000001
        c = (interval(1,1)+interval(1,2))/2;
        fa = func(interval(1,1));
        fb = func(interval(1,2));
        fc = func(c);
        
        c_values = [c_values, c]; % Save the current c value

      
        % disp(interval)

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
    x = c;
end

%generalized newton function
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


%generalized secant function
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

function [fval,dfdx] = orion_test_func(x)
    fval =  (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
    dfdx =  3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;
end

function [fval,dfdx] = orion_test_func2(x)
    fval = x.^2-2;
    dfdx = 2*x;
end

%% 

%Day 2

f_record = my_recorder.generate_recorder_fun(@orion_test_func);
x0 = 5;
x_root = newton_solver(f_record,x0);
input_list = my_recorder.get_input_list();
%semilogy(1:length(input_list),abs(input_list-x_root),'ko','markerfacecolor','k');
xn    = input_list(1:end-1);
xnp1  = input_list(2:end);
nVals = (1:numel(xn)).';   % iteration index for each pair

num_trials = 1000;
spread     = 2; 
x0_list    = linspace(x_root - spread, x_root + spread, num_trials);

xn_all = []; xnp1_all = []; n_all = []; trial_id = [];

for t = 1:numel(x0_list)
    rec = input_recorder();
    f_rec = rec.generate_recorder_fun(@orion_test_func);
        x_star = newton_solver(f_rec, x0_list(t));        
   

    xs = rec.get_input_list();                          
    if numel(xs) < 2, continue; end

    xn   = xs(1:end-1).';
    xnp1 = xs(2:end).';
    n    = (1:numel(xn)).';
    xn_all    = [xn_all;    xn];
    xnp1_all  = [xnp1_all;  xnp1];
    n_all     = [n_all;     n];
    trial_id  = [trial_id;  t*ones(numel(n),1)];
end

fprintf('Collected %d pairs from %d trials (target root ~ %.12g)\n', ...
    numel(xn_all), numel(unique(trial_id)), x_root);

err_n=abs(xn_all-x_root);
err_np1=abs(xnp1_all-x_root);
figure;
loglog(err_n,err_np1,'ro','markerfacecolor','r','markersize',1)

x_regression=[];
y_regression=[];
for n=1:length(trial_id)
    if err_n(n)>1e-15 && err_n(n)<1e-2 && ...
        err_np1(n)>1e-14 && err_np1(n)<1e-2 && ...
        trial_id(n)>2
        x_regression(end+1) = err_n(n);
        y_regression(end+1) = err_np1(n);
end
end

figure;
loglog(x_regression,y_regression,'ro','markerfacecolor','r','markersize',1)

Y = log(y_regression)';
X1 = log(x_regression)';
X2 = ones(length(X1),1);
%run the regression
coeff_vec = regress(Y,[X1,X2]);
%pull out the coefficients from the fit
p = coeff_vec(1);
k = exp(coeff_vec(2));

xx = logspace(log10(min(x_regression)), log10(max(x_regression)), 200);

yy = k .* xx.^p;

figure;                         
loglog(xx, yy, '-', 'LineWidth', 2);

function [dfdx,d2fdx2] = approximate_derivative(fun,x)
    delta_x = 1e-6;
    f_left = fun(x-delta_x);
    f_0 = fun(x);
    f_right = fun(x+delta_x);
    dfdx = (f_right-f_left)/(2*delta_x);
    d2fdx2 = (f_right-2*f_0+f_left)/(delta_x*delta_x);
end

f = test_func01;

[fp_fd, fpp_fd] = approximate_derivative(f, x_root);     % finite differences at the root
p_theory = 2;
k_theory = abs(0.5 * fpp_fd / fp_fd);

fprintf('Theory:   p = 2, k = %.6g\n', k_theory);
fprintf('Regression: p = %.4f, k = %.6g\n', p, k);


%STEP 8 TODO

%Example template for analysis function
%INPUTS:
%solver_flag: an integer from 1-4 indicating which solver to use
% 1->Bisection 2-> Newton 3->Secant 4->fzero
%fun: the mathematical function that we are using the
% solver to compute the root of
%x_guess0: the initial guess used to compute x_root
%guess_list1: a list of initial guesses for each trial
%guess_list2: a second list of initial guesses for each trial
% if guess_list2 is not needed, then set to zero in input
%filter_list: a list of constants used to filter the collected data

num_trials   = 1000;
guess_list_1 = linspace(-5, 0, num_trials);
guess_list_2 = linspace( 1, 5, num_trials);
newton_list  = linspace(-2, 2, num_trials);  % for Newton

% Bisection:
%convergence_analysis(1, test_func01, 1, guess_list_1, guess_list_2);
% Newton:
%convergence_analysis(2, test_func01, 0, newton_list, []);
% Secant:
%convergence_analysis(3, test_func01, 0, guess_list_1, guess_list_2);



function convergence_analysis(solver_flag, fun, x_guess0, guess_list1, guess_list2)
    my_recorder = input_recorder();
    f_record = my_recorder.generate_recorder_fun(fun);
    x_root = newton_solver(f_record,x_guess0);
    input_list = my_recorder.get_input_list();
    xn    = input_list(1:end-1);
    xnp1  = input_list(2:end);
    nVals = (1:numel(xn)).';
    xn_all = []; xnp1_all = []; n_all = []; trial_id = [];

    for t = 1:numel(guess_list1)
        rec = input_recorder();
        f_rec = rec.generate_recorder_fun(fun); 
        
        if solver_flag==1
            x_star = bisection_solver(f_rec,guess_list(t),guess_list2(t));
        elseif solver_flag==2
            x_star = newton_solver(f_rec, guess_list1(t)); 
        elseif solver_flag==3
            x_star= secant_solver(fun,guess_list(t),guess_list(t));
        elseif solver_flag==4
            %fzero
        end
                               
        xs = rec.get_input_list();                          
        if numel(xs) < 2, continue; end
    
        xn   = xs(1:end-1).';
        xnp1 = xs(2:end).';
        n    = (1:numel(xn)).';
        xn_all    = [xn_all;    xn];
        xnp1_all  = [xnp1_all;  xnp1];
        n_all     = [n_all;     n];
        trial_id  = [trial_id;  t*ones(numel(n),1)];
    end 

    fprintf('Collected %d pairs from %d trials (target root ~ %.12g)\n', ...
        numel(xn_all), numel(unique(trial_id)), x_root);
    
    err_n=abs(xn_all-x_root);
    err_np1=abs(xnp1_all-x_root);
    figure;
    loglog(err_n,err_np1,'ro','markerfacecolor','r','markersize',1)
    
    x_regression=[];
    y_regression=[];
    for n=1:length(trial_id)
        if err_n(n)>1e-15 && err_n(n)<1e-2 && ...
            err_np1(n)>1e-14 && err_np1(n)<1e-2 && ...
            trial_id(n)>2
            x_regression(end+1) = err_n(n);
            y_regression(end+1) = err_np1(n);
        end
    end
    
    figure;
    loglog(x_regression,y_regression,'ro','markerfacecolor','r','markersize',1)
    
    Y = log(y_regression)';
    X1 = log(x_regression)';
    X2 = ones(length(X1),1);
    coeff_vec = regress(Y,[X1,X2]);
    p = coeff_vec(1);
    k = exp(coeff_vec(2));
    xx = logspace(log10(min(x_regression)), log10(max(x_regression)), 200);
    yy = k .* xx.^p;
    figure;                         
    loglog(xx, yy, '-', 'LineWidth', 2);
end

%% 

%day 3
x_in=linspace(0,50,200);
[fvalclass, dfdxvalclass]=test_function03(x_in);
x_guess=27;

x_root_3=newton_solver(@test_function03, x_guess);
newton_success_list=[];
newton_fail_list=[];
secant_success_list=[];
secant_fail_list=[];

for n=1:length(x_in)
    x_guess=x_in(n);
    newton_root_attempt=newton_solver(@test_function03, x_guess);
    %secant_root_attempt=secant_solver(@test_function03, x_guess);
    if abs(x_root_3-newton_root_attempt)<.1
        newton_success_list(end+1)=x_guess;
    else
        newton_fail_list(end+1)=x_guess;
    end
end

[fsuccess_list,dfdxsuccess]=test_function03(newton_success_list);
[ffail_list,dfdxfail]=test_function03(newton_fail_list);

figure; hold on;
plot(x_in, fvalclass, 'k-', 'LineWidth', 1, 'DisplayName','Test Function 3');
yline(0, '--k', 'LineWidth', 1, 'DisplayName','F(x)=0');
scatter(newton_success_list, fsuccess_list, 26, 'g', 'filled', 'DisplayName','converges');
scatter(newton_fail_list, ffail_list, 26, 'r', 'filled', 'DisplayName','fails');
xlabel('initial guess x_0'); ylabel('f(x_0)');
title('Newton convergence by initial guess (sigmoid)');
legend('Location','best');

function [f_val,dfdx] = test_function03(x)
    %global input_list;
    %input_list(:,end+1) = x;
    a = 27.3; b = 2; c = 8.3; d = -3;
    H = exp((x-a)/b);
    dH = H/b;
    L = 1+H;
    dL = dH;
    f_val = c*H./L+d;
    dfdx = c*(L.*dH-H.*dL)./(L.^2);
end