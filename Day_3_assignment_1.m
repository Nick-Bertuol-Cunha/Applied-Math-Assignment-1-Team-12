x_in=linspace(0,50,200);
x_in_2=linspace(50,0,200);
[fvalclass, dfdxvalclass]=test_function03(x_in);
x_guess=27;
x_guess_2=20;

x_root_3=newton_solver(@test_function03, x_guess);
newton_success_list=[];
newton_fail_list=[];
fzero_success_list=[];
fzero_fail_list=[];
secant_success_list=[];
secant_success_list_2=[];
secant_fail_list=[];
secant_fail_list_2=[];

for n=1:length(x_in)
    x_guess=x_in(n);
    x_guess_2=x_in_2(n);
    newton_root_attempt=newton_solver(@test_function03, x_guess);
    fzero_root_attempt=fzero(@test_function03, x_guess);
    secant_root_attempt=secant_solver(@test_function03, x_guess, x_guess_2);
    if abs(x_root_3-newton_root_attempt)<.1
        newton_success_list(end+1)=x_guess;
    else
        newton_fail_list(end+1)=x_guess;
    end
    if abs(x_root_3-fzero_root_attempt)<.1
        fzero_success_list(end+1)=x_guess;
    else
        fzero_fail_list(end+1)=x_guess;
    end
    if abs(x_root_3-secant_root_attempt)<.1
        secant_success_list(end+1)=x_guess;
        secant_success_list_2(end+1)=x_guess_2;
    else
        secant_fail_list(end+1)=x_guess;
        secant_fail_list_2(end+1)=x_guess_2;
    end
end

%newton
[fsuccess_list,~]=test_function03(newton_success_list);
[ffail_list,~]=test_function03(newton_fail_list);

figure; hold on;
plot(x_in, fvalclass, 'k-', 'LineWidth', 1, 'DisplayName','Test Function 3');
yline(0, '--k', 'LineWidth', 1, 'DisplayName','F(x)=0');
scatter(newton_success_list, fsuccess_list, 26, 'g', 'filled', 'DisplayName','converges');
scatter(newton_fail_list, ffail_list, 26, 'r', 'filled', 'DisplayName','fails');
xlabel('initial guess x_0'); ylabel('f(x_0)');
title('Newton Convergence by Initial Guess');
legend('Location','best');

%fzero
[fsuccess_list,~]=test_function03(fzero_success_list);
[ffail_list,~]=test_function03(fzero_fail_list);

figure; hold on;
plot(x_in, fvalclass, 'k-', 'LineWidth', 1, 'DisplayName','Test Function 3');
yline(0, '--k', 'LineWidth', 1, 'DisplayName','F(x)=0');
scatter(fzero_success_list, fsuccess_list, 26, 'g', 'filled', 'DisplayName','converges');
scatter(fzero_fail_list, ffail_list, 26, 'r', 'filled', 'DisplayName','fails');
xlabel('initial guess x_0'); ylabel('f(x_0)');
title('fzero Convergence by Initial Guess');
legend('Location','best');

%secant
figure; hold on;
plot(x_in, x_in_2, 'k-', 'LineWidth', 1, 'DisplayName','Test Function 3');
yline(0, '--k', 'LineWidth', 1, 'DisplayName','F(x)=0');
scatter(secant_success_list, secant_success_list_2, 26, 'g', 'filled', 'DisplayName','converges');
scatter(secant_fail_list, secant_fail_list_2, 26, 'r', 'filled', 'DisplayName','fails');
xlabel('initial guess x_0'); ylabel('initial guess x_1');
title('Secant Convergence by Initial Guesses');
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