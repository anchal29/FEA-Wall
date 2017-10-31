function [nodal_disp, nodal_vel, nodal_acc] = apply_newmarks(eff_stiff, global_mass, force_vec, nodal_disp, nodal_vel, nodal_acc, time_step, time_index, L, D)
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
disp('Newmarks method preprocessing:');
tic
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

eff_load = force_vec + global_mass*(a(1)*nodal_disp(:, :, time_index)+a(3)*nodal_vel(:, :, time_index)+a(4)*nodal_acc(:, :, time_index));
toc
disp('Done!');
disp('Solving equillibrium equation....');
tic
nodal_disp(:, :, time_index+1) = L*D*L'\eff_load;
nodal_acc(:, :, time_index+1) = a(1)*(nodal_disp(:, :, time_index+1) - nodal_disp(:, :, time_index)) - a(3)*nodal_vel(:, :, time_index) - a(4)*nodal_acc(:, :, time_index);
nodal_vel(:, :, time_index+1)= nodal_disp(:, :, time_index) - a(7)*nodal_acc(:, :, time_index) - a(8)*nodal_acc(:, :, time_index+1);
toc
disp('Done!');

disp('Solving equillibrium equation....');
tic
nodal_disp(:, :, time_index+1) = eff_stiff\eff_load;
nodal_acc(:, :, time_index+1) = a(1)*(nodal_disp(:, :, time_index+1) - nodal_disp(:, :, time_index)) - a(3)*nodal_vel(:, :, time_index) - a(4)*nodal_acc(:, :, time_index);
nodal_vel(:, :, time_index+1)= nodal_disp(:, :, time_index) - a(7)*nodal_acc(:, :, time_index) - a(8)*nodal_acc(:, :, time_index+1);
toc
disp('Done!');

disp('Solving equillibrium equation....');
tic
LSopts.SYM = true;
LSopts.POSDEF = true;
nodal_disp(:, :, time_index+1) = linsolve(eff_stiff, eff_load, LSopts);
nodal_acc(:, :, time_index+1) = a(1)*(nodal_disp(:, :, time_index+1) - nodal_disp(:, :, time_index)) - a(3)*nodal_vel(:, :, time_index) - a(4)*nodal_acc(:, :, time_index);
nodal_vel(:, :, time_index+1)= nodal_disp(:, :, time_index) - a(7)*nodal_acc(:, :, time_index) - a(8)*nodal_acc(:, :, time_index+1);
toc
disp('Done!');

end