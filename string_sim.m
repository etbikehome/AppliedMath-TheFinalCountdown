function string_sim()
    num_masses = 500;
    total_mass = 2;
    tension_force = 2;
    string_length = 3;
    damping_coeff = 0.005;

    dx = string_length/(num_masses+1);

    amplitude_Uf = 0.2;
    omega_Uf = pi;

    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);

%     Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
%     dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);

    Uf_func_in = @triangle_pulse;
    dUfdt_func_in = @triangle_pulse_derivative;

    width = 0.5;
    height = 5;

    Uf_func = @(t_in) Uf_func_in(t_in,width,height);
    dUfdt_func = @(t_in) dUfdt_func_in(t_in,width,height);

    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;

    %initial conditions
    U0 = zeros(num_masses,1);
    dUdt0 = zeros(num_masses,1);
    V0 = [U0;dUdt0];
    tspan = [0 10];

    %run the integration
    tlist_in = tspan(1):1/20:tspan(2);
    rate_func_wrapper = @(t,V) string_rate_func01(t,V,string_params);
    [tlist,Vlist] = ode45(rate_func_wrapper,tlist_in,V0);

    c = sqrt(tension_force/(total_mass/string_length));

    % Plot
    ymax_val = max(max(Vlist(:,1:num_masses)));
    ymin_val = min(min(Vlist(:,1:num_masses)));
    dy = ymax_val-ymin_val;
    figure()
    axis([0, string_length, ymin_val-.1*dy, ymax_val+.1*dy])
    hold on
    x_pos = linspace(0,string_length,num_masses+2);
    V_data = [0, Vlist(1,1:num_masses), string_params.Uf_func(tlist(1))];
    string_plot = plot(x_pos, V_data);
    line_plot = xline(string_length);
    frametimes = zeros(length(tlist),1);
    tic
    for i = 1:length(tlist)
        % Move vertical line
        x = string_length-c*tlist(i)+.5*width*c;
        x = mod(x,2*string_length);
        if x > string_length
            x = 2*string_length - x;
        end
        set(line_plot,'Value',x);

        V_data = [0, Vlist(i,1:num_masses), string_params.Uf_func(tlist(i))];
        set(string_plot,'YData',V_data)
        frametimes(i) = double(toc);
        pause(tlist(i) - frametimes(i));
        drawnow
    end
    close all
    figure()
    plot(frametimes - tlist)
    xlabel('Timestep')
    ylabel('Frame Lag (s)')
end