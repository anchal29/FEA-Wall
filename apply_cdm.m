function [nodal_disp, nodal_vel, nodal_acc] = apply_cdm(global_stiff, global_mass, force_vec, nodal_disp, time_step, time_index)
%**************************************************************************
% Apply Central Difference Method at time t.
%**************************************************************************
%
% Input parameters:
% global_stiff        - Global stiffness matrix for the given wall at t.
% global_mass - Global mass matrix for the wall element.
% force_vec   - Force vector at time t.
% time_step   - Choosen time step.
% time_index  - Index.
%
% Output:
% nodal_disp  - Resultant nodal displacement at time t+time_step.
% nodal_vec   - Resultant nodal velocity at time t.
% nodal_acc   - Resultant nodal acceleration at time t.

% For referece follow this initialization.
% total_dof = length(global_stiff);
% nodal_disp = zeros(total_dof, 1);
% nodal_vel = zeros(total_dof, 1);
% nodal_acc = zeros(total_dof, 1);
% time_index = t/time_step+1;

a_zero = 1/(time_step)^2;
a_one = 1/(time_step*2);
a_two = time_step*2;
a_three = 1/a_two;
u_minus_delta_t = u_zero - time_step*v_zero + a_three*a_zero;
damp_mat = 0;
eff_mass = a_zero*global_mass + a_one*damp_mat;
eff_load = force_vec - (global_stiff - a_two*global_mass)*nodal_disp(:, :, time_index) - (a_zero*global_mass - a_one*damp_mat)*nodal_disp(:, :, time_index-1);
nodal_disp(:, :, time_index+1) = eff_mass\eff_load;
nodal_vel(:, :, time_index) = a_zero(nodal_disp(:, :, time_index-1) - 2*nodal_disp(:, :, time_index)  + nodal_disp(:, :, time_index+1));
nodal_acc(:, :, time_index) = a_one(nodal_disp(:, :, time_index-1) + nodal_disp(:, :, time_index+1));
end