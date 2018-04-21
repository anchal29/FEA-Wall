function [extended_nodal_connect, no_nodes, face_nodes, hex_equ_face_nodes, hex_equ_nodes] = getExtendedNodalConnect(mesh_meta_data)
%**************************************************************************
% Returns the nodal connectivity matrix for 20 noded brick element.
%**************************************************************************
%
% Input parameters:
% mesh_meta_data         - Number of divisions in x, y and z direction
%                          respectively in the created mesh.
% Output:
% extended_nodal_connect - Nodal connectivity matrix.
% no_nodes               - Total number of nodes present in 20 noded 
%                          brick element's mesh.

left_side_face_nodes = [];
right_side_face_nodes = [];
bottom_face_nodes = [];
top_face_nodes = [];
hex_equ_face_nodes = {[], [], [], []};
hex_equ_nodes = [];
no_elements = prod(mesh_meta_data);
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
    num_nodes_per_layer = 2*((mesh_meta_data(2)+1)*(mesh_meta_data(3)+1)) + (mesh_meta_data(2))*(mesh_meta_data(3)+1) + (mesh_meta_data(2)+1)*(mesh_meta_data(3));

    node_one = layer*num_nodes_per_layer + 2*row*(mesh_meta_data(2)+1) + row*mesh_meta_data(2) + 2*index - 1;
    % disp(node_one);    
    node_two = node_one + num_nodes_per_layer;
    node_five = node_one + 3*mesh_meta_data(2) + 2;
    node_six = node_two + 3*mesh_meta_data(2) + 2;
    node_twelve = layer*num_nodes_per_layer + num_nodes_per_layer - ((mesh_meta_data(2)+1)*(mesh_meta_data(3)+1)) + row*(mesh_meta_data(2)+ 1) + index;
    node_seventeen = node_two + (mesh_meta_data(2)+1) + (mesh_meta_data(2) - index) + 1;
    node_nineteen = node_one + (mesh_meta_data(2)+1) + (mesh_meta_data(2) - index) + 2; 
    node_fourteen = (node_twelve + 1) + (mesh_meta_data(2)+1);
    node_sixteen = node_fourteen - 1;
    mapping = [node_one,
               node_two,
               node_two + 2,
               node_one + 2,
               node_five,
               node_six,
               node_six + 2,
               node_five + 2,
               node_two + 1,
               node_twelve + 1,
               node_one + 1,
               node_twelve,
               node_six + 1,
               node_fourteen,
               node_five + 1,
               node_sixteen,
               node_seventeen,
               node_seventeen + 1,
               node_nineteen,
               node_nineteen - 1];
    extended_nodal_connect(ii, :) = mapping;
    mapping = mapping.';
    if row == 0
        bottom_face_nodes = [bottom_face_nodes, mapping(1:4), mapping(9:12)];
        hex_equ_face_nodes{1} = [hex_equ_face_nodes{1}, mapping(1:4)];
    end
    if row == mesh_meta_data(3)-1
        top_face_nodes = [top_face_nodes, mapping(5:8), mapping(13:16)];
        hex_equ_face_nodes{3} = [hex_equ_face_nodes{3}, mapping(5:8)];
    end
    if index == 1
        left_side_face_nodes = [left_side_face_nodes, mapping(1:2), mapping(5:6), mapping(12), mapping(16), mapping(17), mapping(20)];
        hex_equ_face_nodes{4} = [hex_equ_face_nodes{4}, mapping(1:2), mapping(5:6)];
    end
    if index == mesh_meta_data(2)
        right_side_face_nodes = [right_side_face_nodes, mapping(3:4), mapping(7:8), mapping(10), mapping(14), mapping(18), mapping(19)];
        hex_equ_face_nodes{2} = [hex_equ_face_nodes{2}, mapping(3:4), mapping(7:8)];
    end
    hex_equ_nodes = [hex_equ_nodes mapping(1:8)];
end

face_nodes{1} = unique(bottom_face_nodes);
face_nodes{2} = unique(right_side_face_nodes);
face_nodes{3} = unique(top_face_nodes);
face_nodes{4} = unique(left_side_face_nodes);

no_nodes = length(unique(extended_nodal_connect(:)));

hex_equ_nodes = unique(hex_equ_nodes);
hex_equ_nodes = [hex_equ_nodes*3, hex_equ_nodes*3-1, hex_equ_nodes*3-2];
hex_equ_nodes = sort(hex_equ_nodes);

hex_equ_face_nodes{1} = unique(hex_equ_face_nodes{1});
hex_equ_face_nodes{2} = unique(hex_equ_face_nodes{2});
hex_equ_face_nodes{3} = unique(hex_equ_face_nodes{3});
hex_equ_face_nodes{4} = unique(hex_equ_face_nodes{4});