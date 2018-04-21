function[element_mapping] = ElementMapping(nodal_connect, no_elements, element)
%**************************************************************************
% Computes map of nodal connectivity to nodal dofs connectivity matrix.
% Node 1 corresponds to 1st, 2nd and 3rd dof and node 2 corresponds to 4th,
% 5th and 6th dof. In similar manner this function calulated the mapping
% for each element.
%**************************************************************************
%
% Input parameters:
% nodal_connect   - Nodal connectivity matrix in which each row has the
%                   nodes present in that element.
%                   Row 1 => Element 1 => Terms in row 1 are the nodes
%                   present in element 1 in counter-clockwise order.
% Output:
% element_mapping - Nodal connectivity matrix mapping after considering
%                   three degree of freedom at each nodes.
%                   Row 1 => Element 1 => Terms in row 1 are the nodes
%                   present in element 1 in counter-clockwise order.

if strcmp(element, '8-noded')
	element_dof = 24;
	num_element_nodes = 8;
elseif strcmp(element, '20-noded')
	element_dof = 60;
	num_element_nodes = 20;
end
% element_dof = 60;

    element_mapping = zeros(no_elements, element_dof);
    for ii = 1:no_elements
        final_mapping = zeros(1, element_dof);
        for jj = 1:num_element_nodes
            final_mapping(3*(jj-1)+1:3*jj) = 3*(nodal_connect(ii, jj)-1)+1:3*nodal_connect(ii, jj);
        end
        element_mapping(ii, :) = final_mapping;
    end

end