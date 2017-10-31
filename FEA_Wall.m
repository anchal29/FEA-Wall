%% Anchal Pandey 
% 201316114
% Conventions used - 
    
% Using hexahedron element Finite Element Analysis of a Wall.

clear variables;
clear global;
clc;
    
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
    steel_E = 2 * 10^5;
    pois_ratio = .3;
    bar_dia = 12; % 12mm diameter bars
    condition = 'all_fixed';
    vertical_spacing = 250; % 250 mm spacing of verticle reinforcement.
    horz_spacing = 300; % 300 mm spacing of horizontal reinforcement.
    vertical_dia = 6; % 6mm bars used for verticle reinforcement.
    horz_dia = 6;
    conc_grade = 25;
    steel_grade = 415;
    % Point load acting at centre of wall.
end

conc_yield_strain = 0.002;
conc_E = 5000 * sqrt(conc_grade);

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
density = 2.5*10^(-6);
g = 9.81;
load_vec = [density*g; 0; 0];
force_type = "body";
tic
nodal_force = getNodalForce(force_type, no_elements, load_vec);
toc
tic
global_force = global_force_vector(nodal_coordinate, nodal_connect, nodal_force, element_mapping);
toc
% Initialize the E vector for each element with elastic modulus of
% elasticity.
element_mod_of_elas = steel_E.*element_type_steel(1:no_elements) + conc_E.*(~element_type_steel(1:no_elements));

%% Mass matrix calculation
tic
density = 1;
[mass, global_mass] = getMassMat(nodal_coordinate, nodal_connect, density, element_mapping);
% The created mass matrix should be close to symmetric. Make it symmetric
% and handle the exception case.
if(max(max(abs(global_mass - global_mass.'))) < 1e-5)
    global_mass = (global_mass + global_mass.')/2;
else
    disp('Variations in mass matrix calculation crossed acceptable error. Aborting...');
    return;
end

% Since this is static analysis thus the global stiffness matrix will not
% vary as per time. Calculate stiffness matrix once.
global_stiff = getGlobalStiff(nodal_coordinate, nodal_connect, element_mod_of_elas);

disp('Assembling force from earthquake ground motion...');
tic
P = load('Elcentro.txt');
P(:, 2) = (P(:,2)*g);
force_time_history = P(:,2);
for i = 1:length(P)
    load_vec_time_history(:, :, i) = force_time_history(i)*global_mass*ones(total_no_nodes*3,1);
end
time_step = 0.02;
toc

disp('Applying boundary conditions...');
tic
[global_stiff_bc, load_bc, global_mass_bc] = boundary_conditions(condition, global_stiff, mesh_meta_data, load_vec_time_history, global_mass);
toc;
disp('Done!');
total_max_strain = zeros(1, no_elements);
max_displ = [];


%% Initialization for applying newmarks's method
total_dof = length(global_stiff_bc);
nodal_disp = zeros(total_dof, 1, length(force_time_history));
nodal_vel = zeros(total_dof, 1, length(force_time_history));
nodal_acc = zeros(total_dof, 1, length(force_time_history));

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
[L, D, ~] = ldl(eff_stiff);
toc


%% Solution
%*********************************************************************
% Above calculated informations are not needed to be calculated again.
%*********************************************************************
for i = 1:length(force_time_history)
    disp('Applying Newmarks...')
    tic
    [nodal_disp, nodal_vel, nodal_acc] = apply_newmarks(eff_stiff, global_mass_bc, load_bc(:, :, i), nodal_disp, nodal_vel, nodal_acc, time_step, i, L, D);
    toc
    disp('Done!');
%     max_displ(end+1) = max(nodal_disp(:, :, i));
    break;
end
