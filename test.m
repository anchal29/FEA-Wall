%% Newmark's test
% Example taken from Bathe textbook to test the Newmark's method subroutine
% Example 9.4
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

%% CDM tests
% Example taken from Bathe textbook to test the Newmark's method subroutine
% Example 9.4
time_step = 0.28;
nodal_disp = zeros(2, 1, 14);
nodal_vel = zeros(2, 1, 14);
nodal_acc = zeros(2, 1, 14);
a = [1/(time_step)^2;
     1/(time_step*2);
     2/time_step^2;
     ((time_step)^2)/2];

global_mass_bc = [2 0; 0 1;];
global_stiff_bc = [6 -2; -2 4;];
loa = [0; 10];
nodal_acc(:, :, 2) = [0; 10];
nodal_disp(:, :, 1) = [0; 0.392;];
eff_mass = a(1)*global_mass_bc;
% The created effective mass matrix should be close to symmetric. Make 
% it symmetric and handle the exception case.
if(max(max(abs(eff_mass - eff_mass.'))) < 1e-5)
    eff_mass= (eff_mass+ eff_mass.')/2;
else
    disp('Variations in effective mass matrix calculation crossed acceptable error. Aborting...');
    return;
end
da = decompositio(neff_mass);
for i = 2:13
    disp('Applying CDM...')
    tic
    [nodal_disp, nodal_vel(:, :, i), nodal_acc(:, :, i)] = apply_cdm(global_stiff_bc, global_mass_bc, loa, nodal_disp, time_step, i, da);
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