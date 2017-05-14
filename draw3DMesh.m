% Use patch to generate the 3D mesh for the given wall.
function draw3DMesh(nodal_coordinate, faces)
patch('Vertices',nodal_coordinate.','Faces',faces,...
    'FaceColor','c');
% material shiny;
% alpha('color');
% alphamap('rampdown');
view(55, 82);
axis off;
axis equal;
cameratoolbar('SetCoordSys','y');
cameratoolbar('setmode','orbit');
end