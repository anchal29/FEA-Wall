function [distinct_elements, distinct_coordinates] = getDistinctElements(nodal_coordinate, nodal_connect, element_mod_of_elas)
%**************************************************************************
% Function to get distinct elements calculated from the provided nodal
% connectivity and nodal coordinate matrices.
%**************************************************************************
% 
% Input parameters:
% nodal_connect       - Nodal connectivity matrix in which each row has the
%                       nodes present in that element.
%                       Row 1 => Element 1 => Terms in row 1 are the nodes
%                       present in element 1 in counter-clockwise order.
% nodal_coordinate    - Nodal coordinate matrix having each column as the
%                       x,y and z coordinate value of the repective node. 
%                       nodal_coordinate(1, 1): y coordinate value for node 1
%                       nodal_coordinate(2, 1): z coordinate value for node 1
%                       nodal_coordinate(3, 1): x coordinate value for node 1
% element_mod_of_elas - A vector containing modulus of elasticity value of
%                       all the elements.
% 
% Output:
% distinct_elements   - Stores distinct element numbers.
% distinct_coordinate - Stores respective distinct elements dimension plus
%                       its E value.
%
% Since we have cuboid elements thus if there dimensions are same we could 
% save some CPU cycle by not calculating stiffness matrix again for them.
% This function returns the distinct elements considering their dimensions
% and E value.

distinct_coordinates = [];
distinct_elements = [];
for i = 1:length(nodal_connect)
    dx = nodal_coordinate(nodal_connect(i, 2), 1) - nodal_coordinate(nodal_connect(i, 1), 1);
    dy = nodal_coordinate(nodal_connect(i, 3), 2) - nodal_coordinate(nodal_connect(i, 2), 2);
    dz = nodal_coordinate(nodal_connect(i, 5), 3) - nodal_coordinate(nodal_connect(i, 1), 3);
    if(i == 1)
        distinct_coordinates(end+1, :) = [dx, dy, dz, element_mod_of_elas(i)];
        distinct_elements(end+1) = i;
    elseif(~ismember([dx, dy, dz, element_mod_of_elas(i)], distinct_coordinates, 'rows'))
        distinct_coordinates(end+1, :) = [dx, dy, dz, element_mod_of_elas(i)];
        distinct_elements(end+1) = i;
    end
end
end