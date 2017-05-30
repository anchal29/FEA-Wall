function draw3DMesh(nodal_coordinate, faces)
%**************************************************************************
% Use patch to generate the 3D mesh for the given wall.
%**************************************************************************
%
% Input parameters:
% nodal_coordinate - @todo
% faces            - @todo

figure;
patch('Vertices',nodal_coordinate,'Faces',faces,...
    'FaceColor','c');
axis off;
axis equal;
cameratoolbar('SetCoordSys','z');
cameratoolbar('setmode','orbit');
ax = gca;
ax.CameraPosition = [-179534.00986074534 128291.59093640426 128119.23378678535];
ax.CameraUpVector = [0.40957602214449595 -0.28678821817552286 0.8660254037844387];
ax.CameraViewAngle = 1.4342442304781926;
end