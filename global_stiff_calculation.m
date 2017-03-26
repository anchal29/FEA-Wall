% Computes global stiffness matrix. Appends the element stiffness matrix
% into global stiffness matrix.

function [global_stiff] = global_stiff_calculation(mesh_meta_data, global_stiff, stiff, no_elements)

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

    %% Mapping of degree of freedom with global
    for i = 1:8
        final_mapping(3*(i-1)+1:3*i) = ((mapping(i) - 1) * 3) + 1 : mapping(i) * 3;
    end
%     disp(final_mapping);
    for i = 1:24
        for j = 1:24
            global_stiff(final_mapping(i), final_mapping(j)) = global_stiff(final_mapping(i), final_mapping(j)) + stiff(i, j);
        end
    end
%     disp(ii);
end

end