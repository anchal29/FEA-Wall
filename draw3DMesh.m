% Use patch to generate the 3D mesh for the given wall.
function draw3DMesh(nodal_coordinate, faces)
patch('Vertices',nodal_coordinate.','Faces',faces,...
    'FaceColor','c');
axis off;
axis equal;
cameratoolbar('SetCoordSys','y');
cameratoolbar('setmode','orbit');
ax = gca;
ax.CameraPosition = [-10016.766723143171 23693.23321378274 19389.13053720843];
ax.CameraUpVector = [0.378 0.719 -0.583];
% ax.
end