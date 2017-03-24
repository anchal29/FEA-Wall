%% Anchal Pandey
% 201316114
% Conventions used - 

% Using hexahedron element Finite Element Analysis of a Wall.
%      (z)                                   (5) ___________ (8)  
%       |                                       |\          |\
%       |                                       | \         | \
%       |                                       |  \________|__\    
%       |__________(y)                          |  |(6)     |  | (7)
%      /                                    (1) |__|________|  | 
%     /                                         \  |     (4)\  |   
%    /                                           \ |         \ |   
%  (x)                                         (2)\|__________\|(3)        

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
    
    div_x = input('Enter the division in x direction:    ');
    div_y = input('Enter the division in y direction:    ');
    div_z = input('Enter the division in z direction:    ');
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
    div_y = 24;
    div_z = 24;
end

% % Converting all the input in SI unit 
% thickness = thickness * 10^(-3);
% mod_of_elas= mod_of_elas * 10^6;
% bar_dia = bar_dia * 10^-3;

% Dimensions matrix is in mm so to avoid float values.
dimension = [thickness, width, height];
mod_of_elas  = mod_of_elas * 10^6;
a = mod_of_elas * (1-pois_ratio) / ((1- 2 * pois_ratio) * (1 + pois_ratio));
b = mod_of_elas * pois_ratio / ((1- 2 * pois_ratio) * (1 + pois_ratio));
G = mod_of_elas / (2 * (1 + pois_ratio));
D = [ a b b 0 0 0;
      b a b 0 0 0;
      b b a 0 0 0;
      0 0 0 G 0 0;
      0 0 0 0 G 0;
      0 0 0 0 0 G;
    ];


%% Mesh creation

% Our mesh size will be lower than or equal to size of diameter of the bars
% such that one bar could come inside the cube element. For ease of
% calculation considering the bar to be square of side equal to 

% Here simple assuming the cube side to be equal to the diameter of bar for
% ease.
cube_side_temp = bar_dia;

% Total number of elements will be equal to the multiplication of the
% divisions in all the directions.
no_elements = div_x * div_y * div_z;

% Mesh size contains the mesh sizes in the three directions i.e. x,y and z
% respectively.
mesh_size = [dimension(1)/div_x, dimension(2)/div_y, dimension(3)/div_z]*10^-3;

% mesh_meta_data contains no. of division in height, width and depth
% direction.
mesh_meta_data = [div_x, div_y, div_z];
total_no_nodes = (div_x + 1)*(div_y + 1)*(div_z + 1);
global_stiff = zeros(total_no_nodes*3, total_no_nodes*3);
teemp11 = zeros(24, 24);

% @Todo remove this part after confirming.
for ii = 1:no_elements
    element_no = ii;
    [jacobian, jacobian_testing, nodal_coordinates, pre_stiff, stiff, global_stiff] = octa_element_stiff(mesh_size, element_no, dimension, mesh_meta_data, D, global_stiff);
    if(jacobian == jacobian_testing)
        just_checking(ii) = 1;
    else
        just_checking(ii) = 0;
    end
    if(teemp11 == stiff)
        just_checking_stiff(ii) = 1;
    else
        just_checking_stiff(ii) = 0;
    end
    teemp11 = stiff;
    disp(ii);
end


%% Boundary conditions

%Fixed from all the sides
condition = 'all_fixed';
displacement = sym('displacement', [total_no_nodes*3 1]);
load = zeros(total_no_nodes*3, 1);
load(floor(37), 1) = 1000;

[displacement_, global_stiff_, load_] = boundary_conditions(displacement, condition, global_stiff, mesh_meta_data, load);

displacement = solve(global_stiff_*displacement_ == load_, displacement_);


displacemen	= vpa(struct2array(displacement));
a = size(displacemen);
counter_new = 1;
displ_mesh = zeros(div_y+1, div_z+1);
for i = 2:div_y
    for j = 2:div_z
        displ_mesh(i, j) = displacemen(counter_new);
        counter_new = counter_new + 3;
    end
end    % Loads will be considered later.
figure;
contourf(1:div_y+1, 1:div_z+1, displ_mesh);
colorbar;
figure;
surf(1:div_y+1, 1:div_z+1, displ_mesh);
colorbar;   