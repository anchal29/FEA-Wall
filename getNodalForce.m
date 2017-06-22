function nodal_force = getNodalForce(force_type, no_elements, load_value)
%**************************************************************************
% Computes nodal force from provided constant force value.
%**************************************************************************
%
% Input parameters:
% force_type  - String storing force type i.e. body/surface/point load.
% no_elements - Total number of elements.
% load_value  - Load value vector consisting of constant force value in x,
%               y and z direction
%

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
%%
nodal_force = zeros(24, 1, no_elements);
for ii = 1:no_elements
    if(force_type == 'body')
        [gaussian_points, weights] = gauss_quadrature(2);
        for jj = 1:length(weights)
            temp = gaussian_points;
            ele_load = subs(shape_function_matrix, [zeta, eta, nu], [temp(jj, 1), temp(jj, 2), temp(jj, 3)])*load_value;
            nodal_force(:, :, ii) = nodal_force(:, :, ii) + double(ele_load);
        end

    end
    if(force_type == 'surface')
        [gaussian_points, weights] = gauss_quadrature(2);
        for jj = 1:length(weights)
            temp = gaussian_points;
            ele_load = subs(shape_function_matrix, [zeta, eta, nu], [temp(jj, 1), temp(jj, 2), temp(jj, 3)])*load_value;
            nodal_force(:, :, ii) = nodal_force(:, :, ii) + double(ele_load);
        end

    end
end