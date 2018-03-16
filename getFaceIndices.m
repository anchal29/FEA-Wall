function [face_indices, face_area] = getFaceIndices(mesh_meta_data, face_nun, dimension)
%**************************************************************************
% This function gets indeices of all the selected face.
%**************************************************************************
%
% Input parameters:
% mesh_meta_data - Mesh meta data consists of number division in all the
%                  directions
% face_num     - Face number from 1 to 6.
% Output:

% Figure for reference...
%      _____________________
%     |\ _______Face-3______|\
%     | |                   | |
%     | |                   | |
%     |4|    5 (Forward)    |2|
%     | |    6(Backward)    | |
%     | |                   | |
%(0,0)|_|___________________| |
%      \|_____ Face-1 _______\|
% 
indices = [];
temp2 = (mesh_meta_data(3)+1)*(mesh_meta_data(2)+1)*3; % Stores number of indices present in one verticle strip.
switch face_nun
    case 1
        for ii = 1:mesh_meta_data(1)+1
            indices = [indices (ii-1)*temp2+1:(ii-1)*temp2 + (mesh_meta_data(2) +1)*3];
        end
        face_area = dimension(1)*dimension(2);
    case 2
        for ii = 1:mesh_meta_data(1)+1
            for jj = 1:mesh_meta_data(3)
                temp3 = jj*(mesh_meta_data(2)+1)*3 + temp2*(ii-1) - 3;
                indices = [indices temp3+1:temp3+3];
            end
        end
        face_area = dimension(1)*dimension(3);
    case 3
        for ii = 1:mesh_meta_data(1)+1
            temp = (ii*mesh_meta_data(3)+ii-1)*(mesh_meta_data(2)+1)*3;
            indices = [indices (ii-1)*temp2+1:(ii-1)*temp2 + (mesh_meta_data(2)+1)*3];
        end
        face_area = dimension(1)*dimension(2);
    case 4
        for ii = 1:mesh_meta_data(1)+1
            for jj = 1:mesh_meta_data(3)
                temp3 = (jj-1)*(mesh_meta_data(2)+1)*3 + temp2*(ii-1) - 3;
                indices = [indices temp3+1:temp3+3];
            end
        end
        face_area = dimension(1)*dimension(3);
    case 5
        indices = mesh_meta_data(1)*temp2+1:(mesh_meta_data(1)+1)*temp2;
        face_area = dimension(2)*dimension(3);
    case 6
        indices = 1:temp2;
        face_area = dimension(2)*dimension(3);
end
face_indices = [indices(1:3:end); indices(2:3:end); indices(3:3:end)];