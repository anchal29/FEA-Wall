function [displacement] =  boundary_conditions(displacement, condition, global_stiff,mesh_meta_data, mesh_size)
    switch condition
        case 'all_fixed'
            for i = mesh_meta_data(1)+1:1
                temp = (i*mesh_meta_data(3)+i-1)*(mesh_meta_data(2)+1)*3;
                temp2 = (mesh_meta_data(3)+1)*(mesh_meta_data(2)+1)*3;
                displacement = [displacement(1:temp), displacement(i*temp2 + 1:end)];
                for j = mesh_meta_data(3):2
                    temp = j*(mesh_meta_data(2)+1)*3 + temp2*(i-1) - 3;
                    displacement = [displacement(1:temp); displacement(temp+4:end);];
                    
                    temp = temp - (mesh_meta_data(3)+1)*3;
                    displacement = [displacement(1:temp); displacement(temp+4:end);];
                end
                displacement = displacement(1:(i-1)*temp2); displacement((i-1)*temp2 + (mesh_meta_data(2) +1)*3 + 1);
            end
        case 'top_free'
            
    end
            
end
