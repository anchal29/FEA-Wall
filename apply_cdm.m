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
a = [1/(time_step)^2;
     1/(time_step*2);
     2/time_step^2;
     ((time_step)^2)/2];
eff_mass = a(1)*global_mass;
% The created effective mass matrix should be close to symmetric. Make 
% it symmetric and handle the exception case.
if(max(max(abs(eff_mass - eff_mass.'))) < 1e-5)
    eff_mass= (eff_mass+ eff_mass.')/2;
else
    disp('Variations in effective mass matrix calculation crossed acceptable error. Aborting...');
    return;
end

eff_load = force_vec - (global_stiff - a(3)*global_mass)*nodal_disp(:, :, time_index) - (a(1)*global_mass)*nodal_disp(:, :, time_index-1);
nodal_disp(:, :, time_index+1) = eff_mass\eff_load;
nodal_vel = a(1)*(nodal_disp(:, :, time_index-1) - 2*nodal_disp(:, :, time_index)  + nodal_disp(:, :, time_index+1));
nodal_acc = a(2)*(nodal_disp(:, :, time_index-1) + nodal_disp(:, :, time_index+1));
end