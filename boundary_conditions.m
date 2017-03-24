% This function returns applies the given boundary condition to
% displacement, load and stiffness matrices.

function [displacement, global_stiff, load] =  boundary_conditions(displacement, condition, global_stiff,mesh_meta_data, load)
    switch condition
        case 'all_fixed'
            for i = mesh_meta_data(1)+1:-1:1
                temp = (i*mesh_meta_data(3)+i-1)*(mesh_meta_data(2)+1)*3;
                temp2 = (mesh_meta_data(3)+1)*(mesh_meta_data(2)+1)*3;
                displacement = [displacement(1:temp); displacement(i*temp2 + 1:end)];
                for j = mesh_meta_data(3):-1:2
                    temp = j*(mesh_meta_data(2)+1)*3 + temp2*(i-1) - 3;
                    displacement = [displacement(1:temp); displacement(temp+4:end);];
                    
                    temp = temp - (mesh_meta_data(3)+1)*3;
                    displacement = [displacement(1:temp); displacement(temp+4:end);];
                end
                displacement = [displacement(1:(i-1)*temp2); displacement((i-1)*temp2 + (mesh_meta_data(2) +1)*3 + 1: end)];
            end
            disp('HOla');
            % For Global Stiffness Matrix
            for i = mesh_meta_data(1)+1:-1:1
                temp = (i*mesh_meta_data(3)+i-1)*(mesh_meta_data(2)+1)*3;
                temp2 = (mesh_meta_data(3)+1)*(mesh_meta_data(2)+1)*3;
                global_stiff = [global_stiff(:,1:temp) global_stiff(:,i*temp2 + 1:end)];
                for j = mesh_meta_data(3):-1:2
                    temp3 = j*(mesh_meta_data(2)+1)*3 + temp2*(i-1) - 3;
                    global_stiff = [global_stiff(:,1:temp3) global_stiff(:,temp3+4:end)];
                    
                    temp3 = temp3 - (mesh_meta_data(3)+1)*3;
                    global_stiff = [global_stiff(:,1:temp3) global_stiff(:,temp3+4:end)];
                end
                global_stiff = [global_stiff(:,1:(i-1)*temp2) global_stiff(:,(i-1)*temp2 + (mesh_meta_data(2) +1)*3 + 1: end)];
            end

            for i = mesh_meta_data(1)+1:-1:1
                temp = (i*mesh_meta_data(3)+i-1)*(mesh_meta_data(2)+1)*3;
                temp2 = (mesh_meta_data(3)+1)*(mesh_meta_data(2)+1)*3;
                global_stiff = [global_stiff(1:temp, :); global_stiff(i*temp2 + 1:end, :)];
                for j = mesh_meta_data(3):-1:2
                    temp = j*(mesh_meta_data(2)+1)*3 + temp2*(i-1) - 3;
                    global_stiff = [global_stiff(1:temp, :); global_stiff(temp+4:end, :);];
                    
                    temp = temp - (mesh_meta_data(3)+1)*3;
                    global_stiff = [global_stiff(1:temp, :); global_stiff(temp+4:end, :);];
                end
                global_stiff = [global_stiff(1:(i-1)*temp2, :); global_stiff((i-1)*temp2 + (mesh_meta_data(2) +1)*3 + 1:end , :)];
            end
            
            for i = mesh_meta_data(1)+1:-1:1
                temp = (i*mesh_meta_data(3)+i-1)*(mesh_meta_data(2)+1)*3;
                temp2 = (mesh_meta_data(3)+1)*(mesh_meta_data(2)+1)*3;
                load = [load(1:temp); load(i*temp2 + 1:end)];
                for j = mesh_meta_data(3):-1:2
                    temp = j*(mesh_meta_data(2)+1)*3 + temp2*(i-1) - 3;
                    load = [load(1:temp); load(temp+4:end);];
                    
                    temp = temp - (mesh_meta_data(3)+1)*3;
                    load = [load(1:temp); load(temp+4:end);];
                end
                load = [load(1:(i-1)*temp2); load((i-1)*temp2 + (mesh_meta_data(2) +1)*3 + 1: end)];
            end

    end
end

%%%%%%%%%%%%%%%%%%%%%%
% Consider the following example, mesh numbering is considered in the
% similar manner:
%
% This shows the numbering of first dof on each nodes( We will have three
% dof at each nodes so here just first dof is written in x-direction)
%
%    25 28 31 34
%    13 16 19 22
%    01 04 07 10  (1st Layer)
%
%    61 64 67 70
%    49 52 55 58
%    37 40 43 46  (2nd Layer)
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% For all sides fixed condition:
%   All the first dof which should become zero would be:
%       1, 4, 7, 10, 13, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 58, 61, 64
%       , 67, 70.
%   All the second dof which should become zero would be:
%       2, 5, 8, 11, 14, 23, 26, 27, 32, 35, 38, 41, 44, 47, 50, 59, 62, 65
%       , 68, 71.
%   So overall these numbers should have zero displacement.
%       01-12, 13-15, 22-24, 25-36 and 
%       37-48, 49-51, 58-60, 61-72.
% The displacements neglected by the following code will be 
% mesh_meta_data = [div_x, div_y, div_z] = [1, 3, 2];
%       
%           i = layer + 1 = 2;
%               temp = (2*2 + 2 - 1)(3+1)*3 = 60;
%               temp2 = 36;
%               [(1:59), 2*36+1 = 73:end]
% 
% 
%%%%%%%%%%%%%%%%%%%%%%
