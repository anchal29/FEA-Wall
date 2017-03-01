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
    
    %%% Loads
    % Loads will be considered later.
    %%%% @TODO
else
    height = 3000;
    width = 5000;
    thickness = 230;
    mod_of_elas = 2 * 10^5;
    pois_ratio = -.3
    bar_dia = 12; % 12mm diameter bars
    div_
end

% % Converting all the input in SI unit 
% thickness = thickness * 10^(-3);
% mod_of_elas= mod_of_elas * 10^6;
% bar_dia = bar_dia * 10^-3;

% Dimensions matrix is in mm so to avoid float values.
dimension = [height, width, thickness];
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

% Mesh size will be equal to the GCD of thickness, height, width and bar
% diameter.
mesh_size = gcd(gcd(dimension(1), dimension(2)), gcd(dimension(3), cube_side_temp));
% mesh_size = gcd(gcd(dimension(1), dimension(2)), dimension(3));
mesh_size = 10;
fprintf('Mesh size: %d mm\n', mesh_size); 
% mesh_size = bar_dia;
no_elements = (height/mesh_size)*(width/mesh_size)*(thickness/mesh_size);
% mesh_meta_data contains no. of division in height, width and depth
% direction.
mesh_meta_data = [dimension(1)/mesh_size, dimension(2)/mesh_size, dimension(3)/mesh_size];
for ii = 1:no_elements
    element_no = ii;
    [ele_stiffness, jacobian, jacobian_testing, nodal_coordinates, pre_stiff, stiff] = octa_element_stiff(mesh_size, element_no, dimension, mesh_meta_data, D);
    if(jacobian == jacobian_testing)
        just_checking(ii) = 1;
    else
        just_checking(ii) = 0;
    end
%     break;
end