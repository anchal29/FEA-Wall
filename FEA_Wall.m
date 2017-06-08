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
    mod_of_elas = input('Enter the Young''s modulus :    ');
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
    
    %%% Loads
    % Loads will be considered later.
    %%%% @TODO
else
    height = 3000;
    width = 5000;
    thickness = 170;
    steel_E = 2 * 10^5;
    pois_ratio = .3;
    bar_dia = 12; % 12mm diameter bars
    condition = 'all_fixed';
    vertical_spacing = 250; % 250 mm soacing of verticle reinforcement.
    horz_spacing = 300; % 300 mm soacing of horizontal reinforcement.
    vertical_dia = 6; % 6mm bars used for verticle reinforcement.
    horz_dia = 6;
    conc_grade = 25;
end


% % Converting all the input in SI unit 
% thickness = thickness * 10^(-3);
% mod_of_elas= mod_of_elas * 10^6;
% bar_dia = bar_dia * 10^-3;

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
draw3DMesh(nodal_coordinate, faces);
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
conc_E = 5000 * sqrt(conc_grade);
element_mod_of_elas = steel_E.*element_type_steel(1:no_elements) + conc_E.*(~element_type_steel(1:no_elements));

%*********************************************************************
% Above calculated informations are not needed to be calculated again.
%*********************************************************************
%% Stiffness Matrix Calculation
disp('Finding out local stiffness matrix for all the distinct elements...');
tic
[distinct_elements, distinct_coordinates] = getDistinctElements(nodal_coordinate, nodal_connect, element_mod_of_elas);

% global_stiff = sparse(total_no_nodes*3, total_no_nodes*3);
teemp11 = zeros(24, 24);

stiff = zeros(1, 24*24, length(distinct_elements));

% Calculating the stiffness matrix once for all the different types of
% element.
for ii = 1:length(distinct_elements)
    temp = nodal_coordinate(nodal_connect(distinct_elements(ii),:).', :);
%     temp = {temp(:).'};
    ele_stiff = getElementStiffness(temp, element_mod_of_elas(distinct_elements(ii)));
%     [ele_stiff, shape_function_matrix] = octa_element_stiff(element_mod_of_elas(distinct_elements(i)), nodal_coordinate(nodal_connect(distinct_elements(i),:).', :));
    stiff(:, :, ii) = ele_stiff(:).';
end
toc
disp('Done!');

%% Calculating the global stiffness matrix
disp('Assembling global stiffness matrix...')
tic
[global_stiff] = global_stiff_calculation(nodal_coordinate, nodal_connect, element_mod_of_elas, distinct_coordinates, distinct_elements, stiff);
toc
disp('Done!');
% Make the global stiffness matrix symmetric and handle any bug in
% calculation of global stiffness matrix calculation.
if(max(max(abs(global_stiff - global_stiff))) < 1e-5)
    global_stiff = (global_stiff + global_stiff.')/2;
else
    disp('Variations in global stiffness matrix crossed acceptable error. Aborting...');
    return;
end

%% Boundary conditions
disp('Applying boundary conditions...');
tic
%Fixed from all the sides
% displacement = sym('displacement', [total_no_nodes*3 1]);
displacement_index = (1:total_no_nodes*3).';
load = zeros(total_no_nodes*3, 1);
load(1:3:end) = 1000;
[displacement_index, global_stiff_, load_] = boundary_conditions(displacement_index, condition, global_stiff, mesh_meta_data, load);
toc;
disp('Done!');

%% Solving linear equation
% Not good to consider inverse so change it later with alternative sparse
% factorization implementaion.
disp('Solving for nodal displacement...');
tic
% load_ = sparse(load_);
reduced_displacement = global_stiff_\load_;
toc
disp('Done!');
nodal_displ = zeros(total_no_nodes*3, 1);
nodal_displ(displacement_index) = reduced_displacement;

%% Finding out stress and strain values for each elements
disp('Calculating element stresses and strains...')
tic
syms zeta eta nu;
for ii = 1:no_elements
    ele_nodal_disp = nodal_displ(element_mapping(ii, :));
%     getStrainB(nodal_coordinate(nodal_connect(ii, :).', :), element_mod_of_elas(ii), zeta, eta, nu);
    [max_ele_strain, ele_strain] = ElementStressStrain(nodal_coordinate(nodal_connect(ii, :).', :), element_mod_of_elas(ii), ele_nodal_disp);
    % Check if element is going into non lienar state or not. If it is then
    % update the element modulus of elasticity or the stress-strain slope.
    % Also do calcualtions here so that we can find out the force value
    % using stress. This will be used to find out the residual force.
end
toc
disp('Done!');
%% Displaying results
disp('Showing results...');
tic
nodal_delta_x = nodal_displ(1:3:length(nodal_displ));
nodal_delta_y = nodal_displ(2:3:length(nodal_displ));
nodal_delta_z = nodal_displ(3:3:length(nodal_displ));
new_nodal_coord = [nodal_delta_x nodal_delta_y nodal_delta_z] + nodal_coordinate;
draw3DMesh(new_nodal_coord, faces);
counter_new = 1;
counter_2 = 1;
displ_mesh = zeros(div_z+1, div_y+1);

for ii = 1:mesh_meta_data(3)+ 1
    for jj = 1:mesh_meta_data(2)+1
        if(displacement_index(counter_new) == counter_2)
            displ_mesh(ii, jj) = reduced_displacement(counter_new);
            counter_2 = counter_2 + 3;
            counter_new = counter_new + 3;
        else
            displ_mesh(ii, jj) = 0;
            counter_2 = counter_2 + 3;
        end
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
toc