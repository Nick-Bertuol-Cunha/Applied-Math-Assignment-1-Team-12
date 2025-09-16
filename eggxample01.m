function eggxample01()
    egg_params = struct();
    egg_params.a = 3; egg_params.b = 2; egg_params.c = .15;
    x0 = 5; y0 = 5; theta = pi/6;
    hold on; axis equal; axis square
    axis([0,10,0,10])
    %plot the origin of the egg frame
    plot(x0,y0,'ro','markerfacecolor','r');
    %compute the perimeter of the egg
    [V_list, ~] = egg_func(linspace(0,1,100),x0,y0,theta,egg_params);
    %plot the perimeter of the egg
    plot(V_list(1,:),V_list(2,:),'k');
    %compute a single point along the egg (s=.8)
    %as well as the tangent vector at that point
    [V_single, G_single] = egg_func(.8,x0,y0,theta,egg_params);
    %plot this single point on the egg
    plot(V_single(1),V_single(2),'ro','markerfacecolor','r');
    %plot this tangent vector on the egg
    vector_scaling = .1;
    tan_vec_x = [V_single(1),V_single(1)+vector_scaling*G_single(1)];
    tan_vec_y = [V_single(2),V_single(2)+vector_scaling*G_single(2)];
    plot(tan_vec_x,tan_vec_y,'g')

    [xmin, xmax, ymin, ymax,~,~] = find_bounding_box(x0, y0, theta, egg_params);
    plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin])

    [tg, tw] = collision_func(@egg_trajectory01, egg_params, 0.0, 30.0);
    %fprintf('t_ground = %.4f, t_wall = %.4f\n', tg, tw);

    mypath1 = 'C:\Users\jvidaurrazaga\OneDrive - Olin College of Engineering\Documents\GitHub\Applied-Math-Assignment-1-Team-12\';
    fname='egg_animation.avi';
    input_fname = [mypath1,fname];
    writerObj = VideoWriter(input_fname);
    open(writerObj);
    fig1 = figure(1);
    %set up the plotting axis
    hold on; axis equal; axis square
    axis([0,30,0,30])
    %initialize the plot of the square
    egg_plot = plot(0,0,'k');
    xline(30,"LineWidth",2)
    yline(0,"LineWidth",2)
   
    if tg<tw
        t_final=tg;
    else
        t_final=tw;
    end

    for t=0:.001:t_final
        %compute the position of the square's center (travelling along ellipse)
        [position_x, position_y,new_theta]=egg_trajectory01(t);
        [V_list, ~] = egg_func(linspace(0,1,100),position_x,position_y,new_theta,egg_params);
        %update the coordinates of the egg plot
        set(egg_plot,'xdata',V_list(1,:),'ydata',V_list(2,:));
        %update the actual plotting window
        drawnow;
        %capture a frame (what is currently plotted)
        current_frame = getframe(fig1);
        %write the frame to the video
        writeVideo(writerObj,current_frame);
    end
    [~, x_max_final,~,y_max_final,x_max_y_final, y_min_x_final]=find_bounding_box(position_x,position_y,new_theta,egg_params);
    

    if t_final==tw
        plot(x_max_final,x_max_y_final,'ro','markerfacecolor','r')
    else
        plot(y_min_x_final,y_max_final,'ro','markerfacecolor','r')
    end
    drawnow;
    %capture a frame (what is currently plotted)
    current_frame = getframe(fig1);
    %write the frame to the video
    writeVideo(writerObj,current_frame);
    close(writerObj);
end


function [xmin, xmax, ymin, ymax, xmax_y, ymin_x] = find_bounding_box(x0, y0, theta, egg_params)
    egg_wrapper_func2_x= @(s_in)egg_wrapper_func1_x(s_in, x0, y0, theta, egg_params);
    egg_wrapper_func2_y= @(s_in)egg_wrapper_func1_y(s_in, x0, y0, theta, egg_params);
    s_guess_list = 0:0.2:1;
    dxtol=1e-14; ytol=1e-14; max_iter=200; dfdxmin = 1e-8;
    x_list=[];
    xy_list=[];
    y_list=[];
    yx_list=[];
    for s_guess = s_guess_list
        s_rootx=fzero(egg_wrapper_func2_x, s_guess);
        [V,~]=egg_func(s_rootx, x0, y0, theta, egg_params);
        x_list(end+1)=V(1);
        xy_list(end+1)=V(2);
        s_rooty=fzero(egg_wrapper_func2_y, s_guess);
        [V,~]=egg_func(s_rooty, x0, y0, theta, egg_params);
        y_list(end+1)=V(2);
        yx_list(end+1)=V(1);
    end

    xmin=100;
    xmax=0;
    ymin=100;
    ymax=0;
    for i = 1:length(x_list)
        if x_list(i)<xmin
            xmin=x_list(i);
        end
        if x_list(i)>xmax
            xmax=x_list(i);
            xmax_y=xy_list(i);
        end
    end
    for i = 1:length(y_list)
        if y_list(i)<ymin
            ymin=y_list(i);
            ymin_x=yx_list(i);
        end
        if y_list(i)>ymax
            ymax=y_list(i);
        end
    end
end

function Gx = egg_wrapper_func1_x(s, x0, y0, theta, egg_params)
    [~, G] = egg_func(s, x0, y0, theta, egg_params);
    Gx=G(1);
end
function Gy = egg_wrapper_func1_y(s, x0, y0, theta, egg_params)
    [~, G] = egg_func(s, x0, y0, theta, egg_params);
    Gy=G(2);
end

function [x0,y0,theta] = egg_trajectory01(t)
x0 = 7*t + 8;
y0 = -6*t.^2 + 20*t + 6;
theta = 5*t;
end

function [t_ground,t_wall] = collision_func(traj_fun, egg_params, y_ground, x_wall)
    g_ground = @(t) ymin_at_t(t) - y_ground;  
    g_wall   = @(t) xmax_at_t(t) - x_wall;    

    t_ground = time_finder(g_ground, 0);
    t_wall   = time_finder(g_wall,   0);

    function val = ymin_at_t(t)
        [x0,y0,theta] = traj_fun(t);
        [~, ~, ymin, ~] = find_bounding_box(x0, y0, theta, egg_params);
        val = ymin;
    end
    function val = xmax_at_t(t)
        [x0,y0,theta] = traj_fun(t);
        [~, xmax, ~, ~] = find_bounding_box(x0, y0, theta, egg_params);
        val = xmax;
    end

    function t_hit = time_finder(g, t0)
        tL = t0; gL = g(t0);
        dt = 0.1;   
        tR = tL + dt; gR = g(tR);     
        while sign(gL)*sign(gR) > 0
            tL = tR; gL = gR;
            dt = dt*2;
            tR = tR + dt;
            gR = g(tR);
        end
        t_hit = fzero(g, [tL tR]);
    end
end



function [V, G] = egg_func(s,x0,y0,theta,egg_params)
    a=egg_params.a; b=egg_params.b; c=egg_params.c;
    x = a*cos(2*pi*s); f = exp(-c*x/2); y = b*sin(2*pi*s).*f;
    dx = -2*pi*a*sin(2*pi*s); df = (-c/2)*f.*dx; dy = 2*pi*b*cos(2*pi*s).*f + b*sin(2*pi*s).*df;
    R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
    V = R*[x;y]+[x0*ones(1,length(theta));y0*ones(1,length(theta))];
    G = R*[dx;dy];
end

