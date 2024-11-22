function mode_shape_analysis()
    num_masses = 20;
    total_mass = 2;
    tension_force = 2;
    string_length = 3;
    damping_coeff = 0.005;

    dx = string_length/(num_masses+1);

    amplitude_Uf = 0.2;
    omega_Uf = pi;

    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);

    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;

    %initial conditions
    U0 = zeros(num_masses,1);
    dUdt0 = zeros(num_masses,1);
    V0 = [U0;dUdt0];
    tspan = [0 10];

    [M_mat,K_mat] = construct_2nd_order_matrices(string_params);
    %Use MATLAB to solve the generalized eigenvalue problem
    [Ur_mat,lambda_mat] = eig(K_mat,M_mat);
    
    mode_num = 1;
    Ur = Ur_mat(:, mode_num);
    omega_Uf = sqrt(lambda_mat(mode_num, mode_num));
    
    string_params.Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    string_params.dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);

    %run the integration
    tlist_in = tspan(1):1/20:tspan(2);
    rate_func_wrapper = @(t,V) string_rate_func01(t,V,string_params);
    [tlist,Vlist] = ode45(rate_func_wrapper,tlist_in,V0);


    ymax_val = max(max(Vlist(:,1:num_masses)));
    ymin_val = min(min(Vlist(:,1:num_masses)));
    dy = ymax_val-ymin_val;
    figure()
    axis([0, num_masses+1, ymin_val-.1*dy, ymax_val+.1*dy])
    hold on
    x_pos = 0:num_masses+1;
    V_data = [0, Vlist(1,1:num_masses), string_params.Uf_func(tlist(1))];
    string_plot = plot(x_pos, V_data, '-o');
    frametimes = zeros(length(tlist),1);
    tic
    for i = 1:length(tlist)
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