function string_simulation_template01()
    num_masses = 20;
    total_mass = 2;
    tension_force = 2;
    string_length = 3;
    damping_coeff = 0.05;

    dx = string_length/(num_masses+1);

    amplitude_Uf = 0.2;
    omega_Uf = pi;

    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);

    Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);

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
    tspan = [0 20];

    %run the integration
    rate_func_wrapper = @(t,V) string_rate_func01(t,V,string_params);
    [tlist,Vlist] = ode45(rate_func_wrapper,tspan,V0);

    figure()
    zoom = 1;
    axis([0 num_masses+1 -zoom zoom])
    hold on
    x_pos = 0:num_masses+1;
    V_data = [0, Vlist(1,1:num_masses), Uf_func(tlist(1))];
    string_plot = plot(x_pos, V_data, '-o');
    frametimes = zeros(length(tlist),1);
    tic
    for i = 1:length(tlist)
        V_data = [0, Vlist(i,1:num_masses), Uf_func(tlist(i))];
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