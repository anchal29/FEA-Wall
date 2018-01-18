function element_type_steel = getElementType(nodal_coordinate, nodal_connect, bar_position, thickness)
%**************************************************************************
% Assigns element type to each elements. Returns a boolean vector storing
% whether the element is steel or not. 0 represents concrete element.
%**************************************************************************
%
% Input parameters:
% nodal_connect      - Nodal connectivity matrix in which each row has the
%                      nodes present in that element.
% nodal_coordinate   - Nodal coordinate matrix having each column as the
%                      x,y and z coordinate value of the repective node. 
% bar_position       - Reinforcement bars position in each directions.
% thickness          - Thickness of the wall used to determine how many
%                      layers of reinforcement is present and their 
%                      orientation.
% Output:
% element_type_steel - A boolean vector storing whether the element is 
%                      steel or not. 0 represents concrete element.

% First and last case denotes that the horizontal steel comes before 
% verticle while moving along the thickness of wall and the middle one 
% shows the otherwise.
bar_position_x = bar_position{1};
if(thickness <= 170)
    horz_bar_pos_x = bar_position_x(1);
    vert_bar_pos_x = bar_position_x(2);
elseif(thickness > 220)
    horz_bar_pos_x = [bar_position_x(2) bar_position_x(3)];
    vert_bar_pos_x = [bar_position_x(1) bar_position_x(4)];
else
    horz_bar_pos_x = [bar_position_x(1) bar_position_x(4)];
    vert_bar_pos_x = [bar_position_x(2) bar_position_x(3)];
end
no_elements = length(nodal_connect);
element_type_steel = zeros(no_elements, 1);
for ii = 1:no_elements
    first_point =  nodal_coordinate(nodal_connect(ii, 1), :);
    if(ismember(first_point(1), horz_bar_pos_x) && ismember(first_point(3), bar_position{3}))
        element_type_steel(ii) = 1;
    elseif(ismember(first_point(1), vert_bar_pos_x) && ismember(first_point(2), bar_position{2}))
        element_type_steel(ii) = 1;
    else
        element_type_steel(ii) = 0;
    end
end
end