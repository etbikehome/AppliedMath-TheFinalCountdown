clear all
close all
clc

total_mass = 2;
tension_force = 2;
string_length = 3;
damping_coeff = 0.005;

Uf_func_in = @triangle_pulse;
dUfdt_func_in = @triangle_pulse_derivative;

width = 0.5;
height = 5;

Uf_func = @(t_in) Uf_func_in(t_in,width,height);
dUfdt_func = @(t_in) dUfdt_func_in(t_in,width,height);

%generate the struct
string_params = struct();
string_params.M = total_mass;
string_params.Uf_func = Uf_func;
string_params.dUfdt_func = dUfdt_func;
string_params.Tf = tension_force;
string_params.L = string_length;
string_params.c = damping_coeff;

mode_num = 2;
c = sqrt(tension_force/(total_mass/string_length));
continuous_x = linspace(0,string_length,500);
[continuous_shape,continuous_freq] = modes(continuous_x,mode_num,string_length,c);
plot(continuous_x,continuous_shape); hold on

masses_list = [5 10 20 50 100];
for i = 1:length(masses_list)
    num_masses = masses_list(i);
    string_params.n = num_masses;
    string_params.dx = string_length/(num_masses+1);

    [M_mat,K_mat] = construct_2nd_order_matrices(string_params);
    %Use MATLAB to solve the generalized eigenvalue problem
    [discrete_shape_mat,discrete_freq_mat] = eig(K_mat,M_mat);
    
    discrete_x{i} = linspace(0,string_length,num_masses);
    discrete_shapes{i} = discrete_shape_mat(:, mode_num);
    discrete_shapes{i} = -1 * discrete_shapes{i}';
    
    plot(discrete_x{i},discrete_shapes{i})
    
    discrete_freq = sqrt(discrete_freq_mat(mode_num,mode_num));
    freq_errors(i) = abs(discrete_freq - continuous_freq);
end

legend(['Continous', string(masses_list)])
xlabel('Position')
ylabel('Displacement')

figure()
semilogy(masses_list,freq_errors,'o--')
title('Discrete Frequency Error')
xlabel('Number of Masses')
ylabel('Frequency Error (rad/s)')

num_masses = 100;
string_params.n = num_masses;
string_params.dx = string_length/(num_masses+1);

for i = 1:num_masses
    [M_mat,K_mat] = construct_2nd_order_matrices(string_params);
    %Use MATLAB to solve the generalized eigenvalue problem
    [discrete_shape_mat,discrete_freq_mat] = eig(K_mat,M_mat);
    
    discrete_harmonics(i) = sqrt(discrete_freq_mat(i,i));

    [~,continuous_freq] = modes(1,i,string_length,c);
    continuous_harmonics(i) = continuous_freq;
end

figure()
plot(continuous_harmonics,'o'); hold on
plot(discrete_harmonics,'o');
legend('Continuous','Discrete')
title('Harmonic Series, ' + string(num_masses) + ' Masses')
xlabel('Mode Index')
ylabel('Frequency (rad/s)')
    



