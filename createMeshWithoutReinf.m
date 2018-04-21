function [nodal_connect, nodal_coordinate, faces, mesh_meta_data] = createMeshWithoutReinf(dimension, divisions)
%**************************************************************************
% Creates mesh for the provided wall.
% Returns the nodal connectivity and nodal coordinate matrices.
%**************************************************************************
%
% Input parameters:
% dimension         - [thickness, width, height];
% divisions         - [div_x, div_y, div_z];
%
% Output:
% nodal_connect     - Nodal connectivity matrix.
% nodal_coordinate  - Nodal coordinate matrix.
% faces             - Contains nodes to be connected to get a cube for each
%                     of the element using patch.
% mesh_meta_data    - Number of divisions in x, y and z direction
%                     respectively in the created mesh.


% Gives x_dim, y_dim and z_dim i.e. the mesh size in all the directions.
mesh_size = floor((dimension./divisions));
for index = 1:3
    for i = 0:mesh_size(index)
        if mod(dimension(index), (mesh_size(index)-i)) == 0
            mesh_size(index) = mesh_size(index)-i;
            break;
        end
    end
end
%% Create nodal connectivity and coordinate matrices
% % Now we need to find out all the coordinate values for each elements.
x_coord = 0:mesh_size(1):dimension(1);
y_coord = 0:mesh_size(2):dimension(2);
z_coord = 0:mesh_size(3):dimension(3);
nodal_coordinate = combvec(y_coord, z_coord, x_coord);

% Reorder the nodal coordinate matrix such that index 1 is x values, 2 for 
% y and 3 for z.
correct_nodal_coordinate = nodal_coordinate.';
correct_nodal_coordinate = correct_nodal_coordinate(:, [3 1 2]);
nodal_coordinate = correct_nodal_coordinate;

no_elements = (length(x_coord) - 1) * (length(y_coord) - 1) * (length(z_coord) - 1);

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
    face = [mapping(1) mapping(2) mapping(6) mapping(5);
            mapping(2) mapping(3) mapping(7) mapping(6);
            mapping(3) mapping(4) mapping(8) mapping(7);
            mapping(4) mapping(1) mapping(5) mapping(8);
            mapping(1) mapping(2) mapping(3) mapping(4);
            mapping(5) mapping(6) mapping(7) mapping(8)];
    nodal_connect(ii, :) = mapping;
    faces(6*(ii-1)+1:6*ii, :) = face;
end
% nodal_connect = 0;
end