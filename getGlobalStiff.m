function [global_stiff, stiff] = getGlobalStiff(nodal_coordinate, nodal_connect, element_mod_of_elas, element_pois_ratio, element_type, element_mapping, total_no_nodes, no_elements)
%**************************************************************************
% Complete subroutine to get global stiffness matrix.
%**************************************************************************

disp('Finding out local stiffness matrix for all the distinct elements...');
tic
[distinct_elements, distinct_coordinates] = getDistinctElements(nodal_coordinate, nodal_connect, element_mod_of_elas);

if strcmp(element_type, '8-noded')
	stiff = zeros(1, 24*24, length(distinct_elements));
%     stiff = single(zeros(1, 24*24, length(distinct_elements)));
elseif strcmp(element_type, '20-noded')
	stiff = zeros(60, 60, length(distinct_elements));
%     stiff = single(zeros(60, 60, length(distinct_elements)));
end

% Calculating the stiffness matrix once for all the different types of
% element.
for ii = 1:length(distinct_elements)
	if strcmp(element_type, '8-noded')
	    temp = nodal_coordinate(nodal_connect(distinct_elements(ii),:).', :);
	    ele_stiff = getElementStiffness(temp, element_mod_of_elas(distinct_elements(ii)), element_pois_ratio(distinct_elements(ii)));
	    stiff(:, :, ii) = ele_stiff(:).';
%         stiff(:, :, ii) = single(ele_stiff(:).');
	elseif strcmp(element_type, '20-noded')
	    temp = nodal_coordinate(nodal_connect(distinct_elements(ii),:).', :);
        dx = nodal_coordinate(nodal_connect(distinct_elements(ii), 2), 1) - nodal_coordinate(nodal_connect(distinct_elements(ii), 1), 1);
        dy = nodal_coordinate(nodal_connect(distinct_elements(ii), 3), 2) - nodal_coordinate(nodal_connect(distinct_elements(ii), 2), 2);
        dz = nodal_coordinate(nodal_connect(distinct_elements(ii), 5), 3) - nodal_coordinate(nodal_connect(distinct_elements(ii), 1), 3);
	    ele_stiff = get20NodedElementStiffness([nodal_coordinate(nodal_connect(distinct_elements(ii), 1), :)], [dx dy dz], element_mod_of_elas(distinct_elements(ii)), element_pois_ratio(distinct_elements(ii)));
	    stiff(:, :, ii) = ele_stiff;
%         stiff(:, :, ii) = single(ele_stiff);
	end
end
% whos stiff;
toc
disp('Done!');

%% Calculating the global stiffness matrix
disp('Assembling global stiffness matrix...')
tic
[global_stiff] = global_stiff_calculation(nodal_coordinate, nodal_connect, element_mod_of_elas, distinct_coordinates, stiff, element_type, element_mapping, total_no_nodes, no_elements);
toc
disp('Done!');
% Make the global stiffness matrix symmetric and handle any bug in
% calculation of global stiffness matrix calculation.
if(max(max(abs(global_stiff - global_stiff.'))) < 1e-5)
    global_stiff = (global_stiff + global_stiff.')/2;
else
    disp('Variations in global stiffness matrix crossed acceptable error. Aborting...');
    return;
end

end