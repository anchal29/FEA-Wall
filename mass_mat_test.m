zeta_at_nodes = [-1, 1, 1, -1, -1, 1, 1, -1];
eta_at_nodes = [-1, -1, 1, 1, -1, -1, 1, 1];
nu_at_nodes = [-1, -1, -1, -1, 1, 1, 1, 1];
syms zeta eta nu;

shape_function_matrix = [];
for i = 1:8
    N(i) = 1/8*(1 + zeta*zeta_at_nodes(i))*(1 + eta*eta_at_nodes(i))*(1 + nu*nu_at_nodes(i));
    shape_function_matrix = [shape_function_matrix,N(i)*eye(3)];
end
%% Mass Matrix Calulation

integ = shape_function_matrix' * shape_function_matrix;
% Mass matrix is the integration of N*N' over the volume which we will be 
% doing by numerical integration.
[gaussian_points, weights] = gauss_quadrature(2);

%% 
% tic
% test_stiff = zeros(24, 24);
% for i = 1:length(weights)
%     temp = gaussian_points;
%     temp_pre_stiff = pre_stiff;
%     temp_pre_stiff = subs(temp_pre_stiff, [zeta, eta, nu], [temp(i, 1), temp(i, 2), temp(i, 3)]);
%     test_stiff = test_stiff + vpa(temp_pre_stiff);
% end
% toc
% Best performance from this instead of using subs. Using matlabFunction
% instead of syms subs as it is faster.
temp_fun = matlabFunction(integ);
el_mass = zeros(24, 24);
for i = 1:length(weights)
    temp = gaussian_points;
    temp_mass = temp_fun(temp(i, 2), temp(i, 3), temp(i, 1));
    el_mass = el_mass + double(temp_mass);
end
