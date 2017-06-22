function global_stiff = getGlobalStiff(nodal_coordinate, nodal_connect, element_mod_of_elas)
%**************************************************************************
% Complete subroutine to get global stiffness matrix.
%**************************************************************************

disp('Finding out local stiffness matrix for all the distinct elements...');
tic
[distinct_elements, distinct_coordinates] = getDistinctElements(nodal_coordinate, nodal_connect, element_mod_of_elas);

stiff = zeros(1, 24*24, length(distinct_elements));

% Calculating the stiffness matrix once for all the different types of
% element.
for ii = 1:length(distinct_elements)
    temp = nodal_coordinate(nodal_connect(distinct_elements(ii),:).', :);
%     temp = {temp(:).'};
    ele_stiff = getElementStiffness(temp, element_mod_of_elas(distinct_elements(ii)));
%     [ele_stiff, shape_function_matrix] = octa_element_stiff(element_mod_of_elas(distinct_elements(i)), nodal_coordinate(nodal_connect(distinct_elements(i),:).', :));
    stiff(:, :, ii) = ele_stiff(:).';
end
whos stiff;
toc
disp('Done!');

%% Calculating the global stiffness matrix
disp('Assembling global stiffness matrix...')
tic
[global_stiff] = global_stiff_calculation(nodal_coordinate, nodal_connect, element_mod_of_elas, distinct_coordinates, stiff);
toc
disp('Done!');
% Make the global stiffness matrix symmetric and handle any bug in
% calculation of global stiffness matrix calculation.
if(max(max(abs(global_stiff - global_stiff))) < 1e-5)
    global_stiff = (global_stiff + global_stiff.')/2;
else
    disp('Variations in global stiffness matrix crossed acceptable error. Aborting...');
    return;
end

end