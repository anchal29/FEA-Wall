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
    thickness = 230;
    mod_of_elas = 2 * 10^5;
    pois_ratio = .3;
    bar_dia = 12; % 12mm diameter bars
    div_x = 1;
    div_y = 15;
    div_z = 15;
    condition = 'two_fixed_opposite';
    vertical_spacing = 250; % 250 mm soacing of verticle reinforcement.
    horz_spacing = 300; % 300 mm soacing of horizontal reinforcement.
    vertical_dia = 6; % 6mm bars used for verticle reinforcement.
    horz_dia = 6;
end


% % Converting all the input in SI unit 
% thickness = thickness * 10^(-3);
% mod_of_elas= mod_of_elas * 10^6;
% bar_dia = bar_dia * 10^-3;

% Dimensions matrix is in mm so to avoid float values.
dimension = [thickness, width, height];
mod_of_elas  = mod_of_elas * 10^6;

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
[nodal_connect, nodal_coordinate, faces] = createMesh(dimension, divisions, reinforcment_info);

%% Draw the mesh
draw3DMesh(nodal_coordinate, faces);

%% Stiffness Matrix Calculation

% Total number of elements will be equal to the the size of the
% nodal_coordinate matrix.
no_elements = length(nodal_connect);
element_mod_of_elas = repmat(mod_of_elas, 1, no_elements);

[distinct_elements, distinct_coordinates] = getDistinctElements(nodal_coordinate, nodal_connect, element_mod_of_elas);

total_no_nodes = length(nodal_coordinate);
global_stiff = sparse(total_no_nodes*3, total_no_nodes*3);
teemp11 = zeros(24, 24);

tic
% Calculating the stiffness matrix once for all the different types of
% element.
for i = 1:length(distinct_elements)
    [stiff(:, :, i)] = octa_element_stiff(element_mod_of_elas(distinct_elements(i)), nodal_coordinate(nodal_connect(distinct_elements(i),:).', :));
end
toc

%%
% Calculating the global stiffness matrix
tic
[global_stiff] = global_stiff_calculation(nodal_coordinate, nodal_connect, element_mod_of_elas, distinct_coordinates, distinct_elements, stiff, global_stiff);
toc
%% Boundary conditions

%Fixed from all the sides
% displacement = sym('displacement', [total_no_nodes*3 1]);
displacement = [1:total_no_nodes*3].';
load = zeros(total_no_nodes*3, 1)*1000;
load(1:3:end) = 1000;

[displacement_, global_stiff_, load_] = boundary_conditions(displacement, condition, global_stiff, mesh_meta_data, load);

% Not good to consider inverse so change it later with alternative sparse
% factorization implementaion.
displacement = inv(global_stiff_)*load_;

displacemen = displacement;
a = size(displacemen);
counter_new = 1;
counter_2 = 1;
displ_mesh = zeros(div_z+1, div_y+1);

% displ_mesh(2:div_y, 2:div_z) = displacemen(1:(div_y+1)*(div_z+1));

for i = 1:div_z + 1
    for j = 1:div_y+1
        if(displacement_(counter_new) == counter_2)
            displ_mesh(i, j) = displacemen(counter_new);
            counter_2 = counter_2 + 3;
            counter_new = counter_new + 3;
        else
            displ_mesh(i, j) = 0;
            counter_2 = counter_2 + 3;
        end
    end
end

% 
% for i = 1:div_z + 1
%     if(i == 1 || i == div_z + 1)
%         for j = 1:div_y+1
%             displ_mesh(i, j) = displacemen(counter_new);
%             counter_new = counter_new + 3;
%         end
%     else
%         for j = 2:div_y
%             displ_mesh(i, j) = displacemen(counter_new);
%             counter_new = counter_new + 3;
%         end
%     end
% end    

figure;
contourf(1:div_y+1, 1:div_z+1, displ_mesh);
colorbar;
figure;
surf(1:div_y+1, 1:div_z+1, displ_mesh);
colorbar;