alpha = 0.25;
delta = 0.5;
time_step = 0.28;
nodal_disp = zeros(2, 1, 14);
nodal_vel = zeros(2, 1, 14);
nodal_acc = zeros(2, 1, 14);
a = [1/(alpha*(time_step)^2);
     delta/(alpha*time_step);
     1/(time_step*alpha);
     (1/(2*alpha)) - 1;
     (delta/alpha) - 1;
     (time_step/2)*((delta/alpha) - 2);
     (time_step*(1 - delta));
     delta*time_step;
];

global_mass_bc = [2 0; 0 1;];
global_stiff_bc = [6 -2; -2 4;];
loa = [0; 10];
nodal_acc(:, :, 1) = [0; 10];

eff_stiff = global_stiff_bc + a(1)*global_mass_bc;
% The created effective stiffness matrix should be close to symmetric. Make 
% it symmetric and handle the exception case.
if(max(max(abs(eff_stiff - eff_stiff.'))) < 1e-5)
    eff_stiff= (eff_stiff+ eff_stiff.')/2;
else
    disp('Variations in effective stiffness matrix calculation crossed acceptable error. Aborting...');
    return;
end
da = decomposition(eff_stiff);
for i = 1:13
    disp('Applying Newmarks...')
    tic
    [nodal_disp, nodal_vel, nodal_acc] = apply_newmarks(eff_stiff, global_mass_bc, loa, nodal_disp, nodal_vel, nodal_acc, time_step, i, da);
    toc
    disp('Done!');
%     max_displ(end+1) = max(nodal_disp(:, :, i));
end
%% Plot
max_displ = [];
for i = 1:13
    max_displ(end+1) = max(nodal_disp(:, :, i));
end
plot(0:13-1, max_displ);
