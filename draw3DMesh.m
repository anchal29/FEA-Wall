% Use patch to generate the 3D mesh for the given wall.
function draw3DMesh(nodal_coordinate, faces)
figure;
patch('Vertices',nodal_coordinate,'Faces',faces,...
    'FaceColor','c');
axis off;
axis equal;
cameratoolbar('SetCoordSys','y');
cameratoolbar('setmode','orbit');
ax = gca;
ax.CameraPosition = [210734.89430597893 -124053.20007221 62763.93011415327];
ax.CameraTarget = [115 2500 1500];
ax.CameraUpVector = [-0.30161994759158284 0.21986843807447823 0.9277301747563999];
% ax.CameraViewAngle = 
% ax.
end