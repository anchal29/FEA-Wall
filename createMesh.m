% Function to create mesh for the provided wall considering the
% reinforcement detailings are also provided.
% Returns the nodal connectivity and nodal coordinate matrices.

% Input parameters:
% dimension = [thickness, width, height];
% divisions = [div_x, div_y, div_z];
% reinforcment_info = [vertical_dia, vertical_spacing, vert_side_cover;
%                      horz_dia    , horz_spacing    , horz_side_cover];

function [nodal_connect, nodal_coordinate] = createMesh(dimension, divisions, reinforcment_info)
% Total number of bars in y and z direction
num_bars = floor((dimension(2:3) - 2*reinforcment_info(:,3).')./reinforcment_info(:,2).');

% Verticle bars will be distributed in y direction while the horizontal
% bars are in z direction repectively.
% Bar position of say nth bar is simply the side cover + n-1 times the
% spacing.

vert_bars_poistion = (0:num_bars(1))*reinforcment_info(1,2) + reinforcment_info(1,3); % Along y-axis
horz_bars_poistion = (0:num_bars(2))*reinforcment_info(2,2) + reinforcment_info(2,3); % Along z-axis

% Gives x_dim, y_dim and z_dim i.e. the mesh size in all the directions.
mesh_size = 50 * floor((dimension./divisions)/50);

%% Create nodal connectivity and coordinate matrices

% Now we need to find out all the coordinate values for each elements.
y_coord = helper(mesh_size(2), vert_bars_poistion, reinforcment_info(1,1), dimension(2));
z_coord = helper(mesh_size(3), horz_bars_poistion, reinforcment_info(2,1), dimension(3));
x_coord = get_x(dimension(1), reinforcment_info, mesh_size(1));
% disp(y_coord);
% disp(z_coord);
% disp(mesh_size);
% disp(x_coord);
nodal_coordinate = combvec(y_coord, z_coord, x_coord);

no_elements = (length(x_coord) - 1) * (length(y_coord) - 1) * (length(z_coord) - 1);
total_no_nodes = length(nodal_coordinate);
% disp(no_elements);
% disp(total_no_nodes);

mesh_meta_data = [length(x_coord) - 1, length(y_coord) - 1, length(z_coord) - 1];
for ii = 1:no_elements
    element_no = ii;    
    layer = floor((element_no-1)/(mesh_meta_data(2)*mesh_meta_data(3)));
    temp = mod(element_no, mesh_meta_data(2)*mesh_meta_data(3));
    if temp
        element_no = temp;
    else
        element_no = mesh_meta_data(2)*mesh_meta_data(3);
    end
    row = floor((element_no-1)/mesh_meta_data(2));
    temp2 = mod(element_no, mesh_meta_data(2));
    if temp2
        index = temp2;
    else
        index = mesh_meta_data(2);
    end

    %% Mapping of local nodes to global nodes
    % First node of element. It is mapped to corresponding x node in global.
    node_one = layer*(mesh_meta_data(2)+1)*(mesh_meta_data(3)+1) + row*(mesh_meta_data(2)+1) + index;
    % disp(node_one);
    node_two = node_one + (mesh_meta_data(2)+1)*(mesh_meta_data(3)+1);
    mapping = [node_one, node_two, node_two + 1, node_one + 1, node_one + mesh_meta_data(2) + 1, node_two + mesh_meta_data(2) + 1, node_two + mesh_meta_data(2) + 2, node_one + mesh_meta_data(2) + 2];
    nodal_connect(ii, :) = mapping;
end
% nodal_connect = 0;
end

function [coordinates] = helper(mesh_size, bars_pos, bar_dia, dimension)
    coordinates = 0;
    index = 1;
    fine_meshing = 1;
    flag = 0;
    while(flag == 0)
        coord_val = min(coordinates(end) + mesh_size, dimension);
        if(index > length(bars_pos))
            coordinates(end + 1) = coord_val;
        elseif(coord_val >= bars_pos(index) - fine_meshing*bar_dia)
            temp = -fine_meshing:fine_meshing;
            % Considering fine_meshing number of elements around the bars.
            % So as to get finer mesh around it.
            coordinates = [coordinates, bars_pos(index) + bar_dia*temp];
            index = index + 1;
        else
            coordinates(end + 1) = coord_val;
        end
        flag = coordinates(end) >= dimension;
    end
end

function[x_coordinates] = get_x(dimension, reinforcment_info, mesh_size)

clear_cover = 40;
fine_meshing = 3;
% First and last case denotes that he horizontal steel comes before 
% verticle while moving along the thickness of wall and the middle one 
% shows the otherwise.
if(dimension(1) <= 170)
    bar_dia = [reinforcment_info(2,1), reinforcment_info(2,1)];
    bars_pos = [floor(dimension/2) - reinforcment_info(2,1), floor(dimension/2)];
elseif(dimension(1) > 220)
    bar_dia = [reinforcment_info(1,1), reinforcment_info(2,1), reinforcment_info(2,1), reinforcment_info(1,1)];
    bars_pos = [clear_cover, clear_cover + reinforcment_info(1,1), dimension - clear_cover - sum(reinforcment_info(:,1)), dimension - clear_cover - reinforcment_info(1,1)];
else
    bar_dia = [reinforcment_info(2,1), reinforcment_info(1,1), reinforcment_info(1,1), reinforcment_info(2,1)];
    bars_pos = [clear_cover, clear_cover + reinforcment_info(2,1), dimension - clear_cover - sum(reinforcment_info(:,1)), dimension - clear_cover - reinforcment_info(2,1)];
end
% disp(bars_pos)
% disp(bar_dia) 
x_coordinates  = helper2(mesh_size, bars_pos, bar_dia, dimension);
end


function [coordinates] = helper2(mesh_size, bars_pos, bar_dia, dimension)
    coordinates = 0;
    index = 1;
    fine_meshing = 1;
    flag = 0;
    choice = 1;
    temp_mat = [-fine_meshing:0;0:fine_meshing];
    while(flag == 0)
        coord_val = min(coordinates(end) + mesh_size, dimension);
%         disp(index)
%         disp(coord_val);
        if(index > length(bars_pos))
            coordinates(end + 1) = coord_val;
        elseif(coord_val >= bars_pos(index) - fine_meshing*bar_dia(index))
            temp = temp_mat(choice, :);
            if(choice == 1)
                choice = 2;
            else
                choice = 1;
            end
            % Considering fine_meshing number of elements around the bars.
            % So as to get finer mesh around it.
            coordinates = [coordinates, bars_pos(index) + bar_dia(index)*temp];
            index = index + 1;
        else
            coordinates(end + 1) = coord_val;
        end
        flag = coordinates(end) >= dimension;
    end
end
