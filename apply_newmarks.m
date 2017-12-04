function [displacement, vel, acc, delta_U] = apply_newmarks(global_mass, force_vec, nodal_disp, nodal_vel, nodal_acc, time_step, time_index, da, eff_stiff, time_log)
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
% disp  - Resultant nodal displacement at time t+time_step for current iteration.
% vec   - Resultant nodal velocity at time t+time_step for current iteration.
% acc   - Resultant nodal acceleration at time t+time_step for current iteration.

if(nargin == 9)
    time_log = 0; %By default time_log stamps will be absent unless specified.
end

% For referece follow this initialization.
% total_dof = length(global_stiff);
% nodal_disp = zeros(total_dof, 1);
% nodal_vel = zeros(total_dof, 1);
% nodal_acc = zeros(total_dof, 1);
% time_index = t/time_step+1;
if(time_log == 1)
    disp('Newmarks method preprocessing:');
    tic
end

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
if(time_log == 1)
    toc
    disp('Done!');
    disp('Solving equillibrium equation....');
    tic
end
delta_U = eff_stiff\eff_load;
displacement = nodal_disp(:, :, time_index+1) + delta_U;
acc = a(1)*(displacement - nodal_disp(:, :, time_index)) - a(3)*nodal_vel(:, :, time_index) - a(4)*nodal_acc(:, :, time_index);
vel = nodal_vel(:, :, time_index) + a(7)*nodal_acc(:, :, time_index) + a(8)*acc;
if(time_log == 1)
    toc
    disp('Done!');
end
end