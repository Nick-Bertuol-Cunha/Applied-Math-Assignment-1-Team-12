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