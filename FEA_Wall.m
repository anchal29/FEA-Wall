%% Anchal Pandey 
% 201316114
% Conventions used - 
    
% Using hexahedron element Finite Element Analysis of a Wall.
% Working and tested on MATLAB R2017b
% Dependancy: Neural Network Toolbox.

clear variables;
clear global;
clc;
mkdir('../Logs');
diary(['../Logs/[',datestr(datetime, 'dd-mmm-yyyy, HH.MM AM'), ']Code_LogsNonLinearDynamic.txt'])


%% Input
choice = 0;
while(~(choice == 1 || choice == 2))
    choice = input('[1]. Provide new input values including wall dimensions and loads, or\n[2]. Use default values\nChoose [1/2]:    ');
end
if choice == 1
    %%% Wall Dimension 
    height = input('Enter the height of the Wall(in mm):    ');
    width = input('Enter the width of the Wall(in mm):    ');
    thickness = input('Enter the thickness of the Wall(in mm):    ');
    
    % For linear analysis take the following youngs modulus. % In Newton 
    % per milli-meter square.
    steel_E = input('Enter the Young''s modulus of steel:    ');
    pois_ratio = input('Enter the possion_ratio:    ');
    
    bar_dia = input('Enter the main steel bar diameter:    ');
    
%     div_x = input('Enter the division in x direction:    ');
%     div_y = input('Enter the division in y direction:    ');
%     div_z = input('Enter the division in z direction:    ');
    
    vertical_spacing = input('Enter the verticle reinforcement spacing:    ');
    horz_spacing = input('Enter the horizontal reinforcement spacing:    ');
    vertical_dia = input('Enter the diameter of verticle reinforcement bars:    ');
    horz_dia = input('Enter the diameter of horzontal reinforcement bars:    ');
    condition = 'all_fixed';

    conc_grade = input('Enter the grade of concrete:    ');
    steel_grade = input('Enter the grade of steel:    ');
    
    %%% Loads
    % One point load acting at centre.
else
    height = 3000;
    width = 5000;
    thickness = 170;
    steel_E = 2 * 10^5 * 10^3;
    pois_ratio = .3;
    bar_dia = 12; % 12mm diameter bars
    condition = 'one_fixed';
    vertical_spacing = 250; % 250 mm spacing of verticle reinforcement.
    horz_spacing = 300; % 300 mm spacing of horizontal reinforcement.
    vertical_dia = 6; % 6mm bars used for verticle reinforcement.
    horz_dia = 6;
    conc_grade = 25;
    steel_grade = 415;
end

conc_yield_strain = 0.002;
conc_E = 5000 * sqrt(conc_grade) * 10^3;

steel_Et = steel_E / 5;
conc_Et = conc_E / 10;
steel_yield_strain = steel_grade / steel_E;

% Dimensions matrix is in mm so to avoid float values.
dimension = [thickness, width, height];

%% Naive assumptions
% It seems that 15 divisions across y and z direction are enough. So
% creating the divisions accordingly.
div_y = 15;
div_z = 15;
div_x = 1; % This should vary according to the bars positioning.
divisions = [div_x, div_y, div_z];
vert_side_cover = 75;% Veticle bars side cover.
horz_side_cover = 75;
reinforcment_info = [vertical_dia, vertical_spacing, vert_side_cover;
                     horz_dia    , horz_spacing    , horz_side_cover];

%% Create mesh
disp('Creating Mesh...');
[nodal_connect, nodal_coordinate, faces, mesh_meta_data, bar_position] = createMesh(dimension, divisions, reinforcment_info);
disp('Done!');

%% Draw the mesh
disp('Plotting Mesh...');
% draw3DMesh(nodal_coordinate, faces);
disp('Done!');

%% Get element type
disp('Getting each element type...');
element_type_steel = getElementType(nodal_coordinate, nodal_connect, bar_position, thickness);
disp('Done!');
%%
% Total number of elements will be equal to the the size of the
% nodal_coordinate matrix.
no_elements = length(nodal_connect);

total_no_nodes = length(nodal_coordinate);
element_mapping = ElementMapping(nodal_connect, no_elements);
% Load vec must be in N, mm units.
% density = 2.5*10^(-6); % It is in kg/mm^3
conc_den = 2.5*10^(-6);
steel_den = 8.05*10^(-6);
g = 9.81;
% load_vec = [density*g; 0; 0];
% force_type = 'body';
% tic
% nodal_force = getNodalForce(force_type, no_elements, load_vec);
% toc
% tic
% global_force = global_force_vector(nodal_coordinate, nodal_connect, nodal_force, element_mapping);
% toc
% Initialize the E vector for each element with elastic modulus of
% elasticity.
element_mod_of_elas = steel_E.*element_type_steel(1:no_elements) + conc_E.*(~element_type_steel(1:no_elements));

%% Mass matrix calculation

tic
% density = 1;
density =  steel_den.*element_type_steel(1:no_elements) + conc_den.*(~element_type_steel(1:no_elements));
[mass, global_mass] = getMassMat(nodal_coordinate, nodal_connect, density, element_mapping);
% The created mass matrix should be close to symmetric. Make it symmetric
% and handle the exception case.
if(max(max(abs(global_mass - global_mass.'))) < 1e-5)
    global_mass = (global_mass + global_mass.')/2;
else
    disp('Variations in mass matrix calculation crossed acceptable error. Aborting...');
    return;
end
%%
% Since this is static analysis thus the global stiffness matrix will not
% vary as per time. Calculate stiffness matrix once.
global_stiff = getGlobalStiff(nodal_coordinate, nodal_connect, element_mod_of_elas);

disp('Assembling force from earthquake ground motion...');
tic
P = load('Elcentro.txt');
P(:, 2) = (P(:,2)*g*10^3);
force_time_history = P(:,2);
force_length = length(force_time_history);
dummy = repmat([1; 0; 0], total_no_nodes, 1);
mass_col = global_mass*dummy;
for i = 1:length(P)
    load_vec_time_history(:, :, i) = force_time_history(i)*mass_col;
end
time_step = 0.02;
toc

disp('Applying boundary conditions...');
tic
[global_stiff_bc, load_bc, global_mass_bc, boundary_pt_indices] = boundary_conditions(condition, global_stiff, mesh_meta_data, load_vec_time_history, global_mass);
all_indices = 1:total_no_nodes*3;
non_boundary_indices = setdiff(all_indices, boundary_pt_indices);
toc;
disp('Done!');

%% Initialization for applying newmarks's method
total_dof = length(global_stiff_bc);
nodal_disp = zeros(total_dof, 1, length(force_time_history));
nodal_vel = zeros(total_dof, 1, length(force_time_history));
nodal_acc = zeros(total_dof, 1, length(force_time_history));
nodal_acc(:, :, 1) = P(1, 2)*repmat([1; 0; 0], total_dof/3, 1);
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
% Damping matrix is zero otherwise add a damping matrix term here.
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
toc
%***************************************************************
% When yield strain are high. This should behave as linear one.
% Just for test purpose.
conc_yield_strain = 100;
steel_yield_strain = 100;
%***************************************************************
tic
%% Solution
overall_disp_now = zeros(total_dof, 1);
yielding_order = zeros(no_elements, force_length);
yield_counter = 1;
% Setting boundary point displacement to be zero
overall_disp_now(boundary_pt_indices) = 0; 
total_max_strain = zeros(1, no_elements);
each_ele_strain = zeros(8, 6, no_elements);
%*********************************************************************
% Above calculated informations are not needed to be calculated again.
%*********************************************************************
h = waitbar(0,'Nonlinear dynamic analysis...');
internal_force = zeros(length(global_mass_bc), 1);
prev_internal_force = internal_force;
delta_U_ini = zeros(length(global_mass_bc), 1);
delta_U = zeros(length(global_mass_bc), 1);
for i = 1:force_length
    h = waitbar(i/force_length);
    eff_force = load_bc(:, :, i) - internal_force;
    nodal_disp(:, :, i+1) = nodal_disp(:, :, i);
    nodal_vel(:, :, i+1) = nodal_vel(:, :, i);
    nodal_acc(:, :, i+1) = nodal_acc(:, :, i);
    iterat_count = 1;
    ini_internal_force = internal_force;
    % Force and energy convergence criterion
    conv_crit_one = norm(eff_force - global_mass_bc*nodal_acc(:,:,i+1)) > 10^(-4)*norm(load_bc(:, :, i) - ini_internal_force - global_mass_bc*nodal_acc(:,:,i));
    conv_crit_two = (delta_U'*(eff_force - global_mass_bc*nodal_acc(:,:,i+1)))/(delta_U_ini'*(load_bc(:, :, i) - ini_internal_force - global_mass_bc*nodal_acc(:,:,i))) > 10^(-3);
    disp([num2str(conv_crit_one), num2str(conv_crit_two)]);
    while(conv_crit_one || conv_crit_two)
        tic
        ['Iteration count:    ', num2str(iterat_count), '     and Residual force is:     ', num2str(norm(eff_force))]
        iterat_count = iterat_count + 1;
        disp('Applying newmark''s method to solve equillibrium equation...');
        [displacement, vel, acc, delta_U] = apply_newmarks(global_mass_bc, eff_force, nodal_disp, nodal_vel, nodal_acc, time_step, i, da, eff_stiff);
        toc
        disp('Done!');
        if(delta_U_ini == 0)
            delta_U_ini = delta_U;
        end
        nodal_disp(:, :, i+1) = displacement;
        nodal_vel(:, :, i+1) = vel;
        nodal_acc(:, :, i+1) = acc;
        % Calculate the stresses and strain.
        % Finding out stress and strain values for each element
        disp('Calculating element stresses and strains...');
        count = 0;
        overall_disp_now(non_boundary_indices) = delta_U;
        for ii = 1:no_elements
            ele_nodal_disp = overall_disp_now(element_mapping(ii, :));
        %     getStrainB(nodal_coordinate(nodal_connect(ii, :).', :), element_mod_of_elas(ii), zeta, eta, nu);
            [max_ele_strain, ele_strain] = ElementStressStrain(nodal_coordinate(nodal_connect(ii, :).', :), ele_nodal_disp);
            total_max_strain(ii) = total_max_strain(ii) + max_ele_strain;
            % Check if element is going into non lienar state or not. If it is then
            % update the element modulus of elasticity or the stress-strain slope.
            % Also do calcualtions here so that we can find out the force value
            % using stress. This will be used to find out the residual force.
            if(total_max_strain(ii) > conc_yield_strain && element_type_steel(ii) == 0)
                element_mod_of_elas(ii) = conc_Et;
                count = count + 1;
                yielding_order(ii, i) = iterat_count;
            elseif(total_max_strain(ii) > steel_yield_strain && element_type_steel(ii) == 1)
                yielding_order(ii, i) = iterat_count;
                element_mod_of_elas(ii) = steel_Et;
                count = count + 1;
            end
            each_ele_strain(:, :, ii) = each_ele_strain(:, :, ii) + ele_strain;
        end
        toc
        disp('Done!');
        fprintf('Number of elements that went into non-linear state in this iteration: %d\n', count);
%         global_stiff = getGlobalStiff(nodal_coordinate, nodal_connect, element_mod_of_elas);
%         global_stiff_bc = boundary_conditions(condition, global_stiff, mesh_meta_data, load_vec_time_history, global_mass);
%         eff_stiff = global_stiff_bc + a(1)*global_mass_bc;
        % The created effective stiffness matrix should be close to symmetric. Make 
        % it symmetric and handle the exception case.
        if(max(max(abs(eff_stiff - eff_stiff.'))) < 1e-5)
            eff_stiff= (eff_stiff+ eff_stiff.')/2;
        else
            disp('Variations in effective stiffness matrix calculation crossed acceptable error. Aborting...');
            return;
        end
        internal_force = internal_force + global_stiff_bc * delta_U;
        eff_force = load_bc(:, :, i) - internal_force;
        conv_crit_one = norm(eff_force - global_mass_bc*nodal_acc(:,:,i+1)) > 10^(-4)*norm(load_bc(:, :, i) - ini_internal_force - global_mass_bc*nodal_acc(:,:,i));
        conv_crit_two = (delta_U'*(eff_force - global_mass_bc*nodal_acc(:,:,i+1)))/(delta_U_ini'*(load_bc(:, :, i) - ini_internal_force - global_mass_bc*nodal_acc(:,:,i))) > 10^(-3);
        toc
    end
end
close(h);
toc
diary off
%% Plot
end_displ = [];
end_force = [];
for i = 1:length(force_time_history)
    end_force(end+1) = load_bc(end-2, :, i);
    end_displ(end+1) = (nodal_disp(end-2, :, i));
end
max_text = ['\leftarrow Max displacement = ',num2str((max(end_displ)/height)*100), '% of the wall height.'];
max_index = find(max_displ == max(max_displ));
figure;
plot(P(:,1), end_displ, 'LineWidth', 1.9);
text(P(max_index, 1), end_displ(max_index), max_text, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (sec)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Displacement (mm)', 'FontSize', 12, 'FontWeight', 'bold');
title('Displacement time history in x-direction for Elcentro GM(Applied in x-direction) - Last Node');
%% Plot
max_now = 0;
for j = 1:3:length(nodal_disp(:, :, 1))
    if(max_now <= max(abs(nodal_disp(j, :, i))))
        max_now = max(abs(nodal_disp(j, :, i)));
%         max_displ_new = [];
%         for i = 1:length(force_time_history)
%             max_displ_new(end+1) = nodal_disp(j, :, i);
%         end
%         plot(0:length(force_time_history)-1, max_displ_new);
%         hold on;
        max_index = j;
    end
end
%% Each level plot
%*********************************************************************
% Caculate first point's dof in x-direction on each of the levels while 
% moving upwards on wall. These points are will give displacement in 
% x-direction.
%*********************************************************************
first_points = zeros(mesh_meta_data(3), 1);
first_points(1) = 1;
for i = 2:mesh_meta_data(3)
    first_points(i) = first_points(i-1) + (mesh_meta_data(2)+1)*3;
end
max_acc_each_time_step = [];
max_disp_each_time_step = [];
for ii = 1:mesh_meta_data(3)
    max_displ = [];
    max_acc = [];
    for i = 1:length(force_time_history)
        max_displ(end+1) = (nodal_disp(first_points(ii), :, i));
        max_acc(end+1) = (nodal_acc(first_points(ii), :, i));
    end
    max_acc_each_time_step(end+1) = max(max_acc);
    max_disp_each_time_step(end+1) = max(max_displ);
end
% Distinct dz values
distinct_dz_coordinates = unique(nodal_coordinate(:, 3));
distinct_dz_coordinates = distinct_dz_coordinates(2:end);
%%
figure;
plot(distinct_dz_coordinates, max_acc_each_time_step/9.8/1000, '-p', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'white', 'LineWidth', 1.5);
xlabel('Height (mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Acceleration (g)', 'FontSize', 12, 'FontWeight', 'bold');
title('Max acceleration at each level');

figure;
plot(distinct_dz_coordinates, max_disp_each_time_step, '-s', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'white', 'LineWidth', 1.5);
xlabel('Height (mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Displacement (mm)', 'FontSize', 12, 'FontWeight', 'bold');
title('Max displacement at each level');
%%
% To find out overall displacement including the removed boundary points
final_disp = zeros(total_dof, 1);
% Setting boundary point displacement to be zero
final_disp(boundary_pt_indices) = 0; 
final_disp(non_boundary_indices) = nodal_disp(:, :, end);
nodal_delta_x = final_disp(1:3:length(final_disp));
nodal_delta_y = final_disp(2:3:length(final_disp));
nodal_delta_z = final_disp(3:3:length(final_disp));
% Here multiplying by a factor just to visualize the deformation.
new_nodal_coord = 100*[nodal_delta_x nodal_delta_y nodal_delta_z] + nodal_coordinate;
draw3DMesh(new_nodal_coord, faces);
counter_1 = 1;
displ_mesh = zeros(mesh_meta_data(3)+1, mesh_meta_data(2)+1);

for ii = 1:mesh_meta_data(3)+ 1
    for jj = 1:mesh_meta_data(2)+1
        displ_mesh(ii, jj) = final_disp(counter_1);
        counter_1 = counter_1 + 3;
    end
end
distinct_z = unique(nodal_coordinate(:, 3), 'rows');
distinct_y = unique(nodal_coordinate(:, 2), 'rows');

figure;
contourf(distinct_y, distinct_z, displ_mesh);
colorbar;
figure;
surf(distinct_y, distinct_z, displ_mesh);
colorbar;
%% Force-displacement curve
end_displ = [];
end_force = [];
for i = 1:length(force_time_history)
    end_force(end+1) = load_bc(1, :, i);
    end_displ(end+1) = (nodal_disp(1, :, i));
end
plot(end_displ, end_force);
%%
h = draw3DMesh(nodal_coordinate, faces);
filename = ['../Logs/[',datestr(datetime, 'dd-mmm-yyyy, HH.MM AM'), ']ResponseAnimatedNonLinear.gif'];
%%
animation_frames(length(force_time_history)) = struct('cdata',[],'colormap',[]);
for i = 1:length(force_time_history)
    % To find out overall displacement including the removed boundary points
    final_disp = zeros(total_dof, 1);
    % Setting boundary point displacement to be zero
    final_disp(boundary_pt_indices) = 0; 
    final_disp(non_boundary_indices) = nodal_disp(:, :, i);
    nodal_delta_x = final_disp(1:3:length(final_disp));
    nodal_delta_y = final_disp(2:3:length(final_disp));
    nodal_delta_z = final_disp(3:3:length(final_disp));
    % Here multiplying by a factor just to visualize the deformation.
    new_nodal_coord = 50*[nodal_delta_x nodal_delta_y nodal_delta_z] + nodal_coordinate;
%     draw3DMesh(new_nodal_coord, faces);
    set(h, 'Vertices',new_nodal_coord);
    drawnow limitrate;
          % Capture the plot as an image 
    animation_frames(i) = getframe(gca);
    frame = getframe(gca); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 

    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount', 0, 'Delay', 0.1); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end
end

%% Save workspace. Takes about 3GB space
save(['../Logs/[NonLinear_Dynamic',datestr(datetime, 'dd-mmm-yyyy, HH.MM AM'), ']Workspace_Variables.txt']);