function create_element_nodal_force_fun()
clear variables;
clear global;
clc;
tic
%% Shape Function
zeta_at_nodes = [-1, 1, 1, -1, -1, 1, 1, -1];
eta_at_nodes = [-1, -1, 1, 1, -1, -1, 1, 1];
nu_at_nodes = [-1, -1, -1, -1, 1, 1, 1, 1];
syms zeta eta nu;
shape_function_matrix = sym(zeros(3, 24));
for ii = 1:8
    N = 1/8*(1 + zeta*zeta_at_nodes(ii))*(1 + eta*eta_at_nodes(ii))*(1 + nu*nu_at_nodes(ii));
    shape_function_matrix(:, 3*(ii-1)+1:3*ii) = N*eye(3);
end
syms load_x load_y load_z;
load_vec = [load_x; load_y; load_z];
nodal_force = zeros(24, 1);
[gaussian_points, weights] = gauss_quadrature(2);
for jj = 1:length(weights)
    temp = gaussian_points;
    shape_fun = shape_function_matrix;
    ele_load = subs(shape_fun, [zeta, eta, nu], [temp(jj, 1), temp(jj, 2), temp(jj, 3)])'*load_vec;
%     ele_load
    nodal_force = nodal_force+ele_load;
end
toc
tic
matlabFunction(nodal_force, 'File', 'getEleNodalBodyForce', 'Optimize' ,true, 'vars', {'load_x', 'load_y', 'load_z'});
toc
disp('Done!!');
