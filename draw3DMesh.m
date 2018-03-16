function h = draw3DMesh(nodal_coordinate, faces)
%**************************************************************************
% Uses patch to draw 3D mesh for the given wall.
%**************************************************************************
%
% Input parameters:
% nodal_coordinate - Nodal coordinate matrix.
% faces            - Face matrix containing face connection for each
%                    elements together for plotting using patch.

figure;
h = patch('Vertices',nodal_coordinate,'Faces',faces,...
    'FaceColor','c');
axis off;
axis equal;
cameratoolbar('SetCoordSys','z');
cameratoolbar('setmode','orbit');
ax = gca;
ax.YAxis.Visible = 'off';
ax.ZAxis.Visible = 'off';
ax.XLim = [-3000 3000];
ax.CameraPosition = [-179534.00986074534 128291.59093640426 128119.23378678535];
ax.CameraUpVector = [0.40957602214449595 -0.28678821817552286 0.8660254037844387];
ax.CameraViewAngle = 1.4342442304781926;
end