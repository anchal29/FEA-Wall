function interpolated_nodal_coord = getInterpolatedNodalCoord(element_nodal_coordinate)
temp = 1:8;
pos_comb = nchoosek(temp, 2);
interpolated_nodal_coord = sym(zeros(20, 3));

interpolated_nodal_coord(1:8, :) = element_nodal_coordinate;

interpolated_nodal_coord(9, :)  = (element_nodal_coordinate(2, :) + element_nodal_coordinate(3, :))/2;
interpolated_nodal_coord(10, :) = (element_nodal_coordinate(3, :) + element_nodal_coordinate(4, :))/2;
interpolated_nodal_coord(11, :) = (element_nodal_coordinate(4, :) + element_nodal_coordinate(1, :))/2;
interpolated_nodal_coord(12, :) = (element_nodal_coordinate(1, :) + element_nodal_coordinate(2, :))/2;

interpolated_nodal_coord(13, :) = (element_nodal_coordinate(6, :) + element_nodal_coordinate(7, :))/2;
interpolated_nodal_coord(14, :) = (element_nodal_coordinate(7, :) + element_nodal_coordinate(8, :))/2;
interpolated_nodal_coord(15, :) = (element_nodal_coordinate(8, :) + element_nodal_coordinate(5, :))/2;
interpolated_nodal_coord(16, :) = (element_nodal_coordinate(5, :) + element_nodal_coordinate(6, :))/2;

interpolated_nodal_coord(17, :) = (element_nodal_coordinate(2, :) + element_nodal_coordinate(6, :))/2;
interpolated_nodal_coord(18, :) = (element_nodal_coordinate(3, :) + element_nodal_coordinate(7, :))/2;
interpolated_nodal_coord(19, :) = (element_nodal_coordinate(4, :) + element_nodal_coordinate(8, :))/2;
interpolated_nodal_coord(20, :) = (element_nodal_coordinate(1, :) + element_nodal_coordinate(5, :))/2;
