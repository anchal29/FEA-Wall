function [distinct_elements, distinct_coordinates] = getDistinctElements(nodal_coordinate, nodal_connect)
% Function to get distinct elements calculated from the provided nodal
% connectivity and nodal coordinate matrices.

%  Input parameters:
%  > nodal_connect - Nodal connectivity matrix in which each row has the nodes
%  present in that element. Row 1 => Element 1 => Terms in row 1 are the
%  nodes present in element 1 in counter-clockwise order.
%  > nodal_coordinate - Nodal coordinate matrix having each column as the x,y
%  and z coordinate value of the repective node. 
%     nodal_coordinate(1, 1): y coordinate value for node 1
%     nodal_coordinate(2, 1): z coordinate value for node 1
%     nodal_coordinate(3, 1): x coordinate value for node 1

distinct_coordinates = [];
distinct_elements = [];
for i = 1:length(nodal_connect)
    dx = nodal_coordinate(3, nodal_connect(i, 2)) - nodal_coordinate(3, nodal_connect(i, 1));
    dy = nodal_coordinate(1, nodal_connect(i, 3)) - nodal_coordinate(1, nodal_connect(i, 2));
    dz = nodal_coordinate(2, nodal_connect(i, 5)) - nodal_coordinate(2, nodal_connect(i, 1));
    if(i == 1)
        distinct_coordinates(end+1, :) = [dx, dy, dz];
        distinct_elements(end+1) = i;
    elseif(~ismember([dx, dy, dz], distinct_coordinates, 'rows'))
        distinct_coordinates(end+1, :) = [dx, dy, dz];
        distinct_elements(end+1) = i;
    end
end
end