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
