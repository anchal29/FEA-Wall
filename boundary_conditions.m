function [global_stiff, load, mass, index, reduced_index, hex_equ_index] =  boundary_conditions(condition, global_stiff, mesh_meta_data, load, mass, element_type, face_nodes, hex_equ_nodes)
%**************************************************************************
% This function applies the given boundary condition over load and 
% stiffness matrices.
%**************************************************************************
%
% Input parameters:
% condition      - Boundary condition to apply.
% global_stiff   - Original global stiffness matrix.
% mesh_meta_data - Mesh meta data consists of number division in all the
%                  directions
% load           - Load vector applied.
% mass           - Original mass matrix.
%
% Output:
% global_stiff   - Global Stiffness Matrix after application of boundary
%                  conditions.
% load           - Load vector after application of BC
% mass           - Mass matrix after application of BC.
% index          - Index of all the boundary points in the original frame.

switch condition
    % Figure for reference...
    %      _____________________
    %     |\ _______Face-3______|\
    %     | |                   | |
    %     | |                   | |
    %     |4|                   |2|
    %     | |                   | |
    %     | |                   | |
    %     |_|___________________| |
    %      \|_____ Face-1 _______\|
    % 
    % All fixed case means face 1,2,3,4 fixed and similarly are other
    % cases defined.
    %
    
    case 'all_fixed'
        face_fixed = [1, 1, 1, 1];
    case 'two_opposite_fixed_1'
        face_fixed = [1, 0, 1, 0];
    case 'two_opposite_fixed_2'
        face_fixed = [0, 1, 0, 1];
    case 'one_fixed'
        face_fixed = [1, 0, 0, 0];
    case 'one_free'
        face_fixed = [1, 1, 0, 1];
end
% face_fixed = [0, 0, 0, 0];
% Instead of removing the index from matrix. Substitute 0 at those index in
% stiffness with digonal values as 1 and in load 0 at each of the index.
% This means we are assuming fixidity at thode nodes.
index = [];
if strcmp(element_type, '8-noded')
    for ii = mesh_meta_data(1)+1:-1:1
        temp = (ii*mesh_meta_data(3)+ii-1)*(mesh_meta_data(2)+1)*3;
        temp2 = (mesh_meta_data(3)+1)*(mesh_meta_data(2)+1)*3;
        if(face_fixed(3))
            index = [index temp+1:ii*temp2];
        end
        for jj = mesh_meta_data(3):-1:2
            temp3 = jj*(mesh_meta_data(2)+1)*3 + temp2*(ii-1) - 3;
                if(face_fixed(2))
                    index = [index temp3+1:temp3+3];
                end
            temp3 = temp3 - (mesh_meta_data(2))*3;
            if(face_fixed(4))
                index = [index temp3+1:temp3+3];
            end
        end
        if(face_fixed(1))
            index = [index (ii-1)*temp2+1:(ii-1)*temp2 + (mesh_meta_data(2) +1)*3];
        end
    end
    reduced_index = [];
    hex_equ_index = [];
elseif strcmp(element_type, '20-noded')
    reduced_index = [];
    if(face_fixed(1))
        index = [index face_nodes{1}];
        reduced_index = [reduced_index hex_equ_nodes{1}];
    end
    if(face_fixed(2))
        index = [index face_nodes{2}];
        reduced_index = [reduced_index hex_equ_nodes{2}];
    end
    if(face_fixed(3))
        index = [index face_nodes{3}];
        reduced_index = [reduced_index hex_equ_nodes{3}];
    end
    if(face_fixed(4))
        index = [index face_nodes{4}];
        reduced_index = [reduced_index hex_equ_nodes{4}];
    end
    index = unique(index);
    index = [index*3-2, index*3-1, index*3];
    reduced_index = unique(reduced_index);
    reduced_index = [reduced_index*3 reduced_index*3-1 reduced_index*3-2];
    dummy = index;
    index = [];
        for ii = mesh_meta_data(1)+1:-1:1
        temp = (ii*mesh_meta_data(3)+ii-1)*(mesh_meta_data(2)+1)*3;
        temp2 = (mesh_meta_data(3)+1)*(mesh_meta_data(2)+1)*3;
        if(face_fixed(3))
            index = [index temp+1:ii*temp2];
        end
        for jj = mesh_meta_data(3):-1:2
            temp3 = jj*(mesh_meta_data(2)+1)*3 + temp2*(ii-1) - 3;
                if(face_fixed(2))
                    index = [index temp3+1:temp3+3];
                end
            temp3 = temp3 - (mesh_meta_data(2))*3;
            if(face_fixed(4))
                index = [index temp3+1:temp3+3];
            end
        end
        if(face_fixed(1))
            index = [index (ii-1)*temp2+1:(ii-1)*temp2 + (mesh_meta_data(2) +1)*3];
        end
    end
    hex_equ_index = index;
    index = dummy;
end
% Apply boundary condition at all the above calculated indices.
global_stiff(:, index) = [];
global_stiff(index, :) = [];
mass(:, index) = [];
mass(index, :) = [];
load(index, :, :) = [];
end