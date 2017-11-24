function [disp, vel, acc] = apply_newmarks(global_mass, force_vec, nodal_disp, nodal_vel, nodal_acc, time_step, time_index, da)
%**************************************************************************
% Apply Newmarks Method at time t.
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
% disp('Newmarks method preprocessing:');
% tic
alpha = 0.25;
delta = 0.5;

a = [1/(alpha*(time_step)^2);
     delta/(alpha*time_step);
     1/(time_step*alpha);
     (1/(2*alpha)) - 1;
     (delta/alpha) - 1;
     (time_step/2)*((delta/alpha) - 2);
     (time_step*(1 - delta));
     delta*time_step;
];

eff_load = force_vec + global_mass*(a(1)*(nodal_disp(:, :, time_index) - nodal_disp(:, :, time_index + 1))+a(3)*nodal_vel(:, :, time_index)+a(4)*nodal_acc(:, :, time_index));
% toc
% disp('Done!');
% disp('Solving equillibrium equation....');
% tic
delta_U = da\eff_load;
disp = nodal_disp(:, :, time_index+1) + delta_U;
acc = a(1)*(disp - nodal_disp(:, :, time_index)) - a(3)*nodal_vel(:, :, time_index) - a(4)*nodal_acc(:, :, time_index);
vel = nodal_vel(:, :, time_index) + a(7)*nodal_acc(:, :, time_index) + a(8)*acc;
% toc
% disp('Done!');
end