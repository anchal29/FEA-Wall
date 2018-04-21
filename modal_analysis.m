%% Anchal Pandey 
% 201316114
% Conventions used - 
    
% Using hexahedron element Finite Element Analysis of a Wall.
% Working and tested on MATLAB R2017b
% Dependancy: Neural Network Toolbox.
clear variables;
clc;
mkdir('../Logs');
folder_name = ['[',datestr(datetime, 'dd-mmm-yyyy, HH.MM AM'), '] Modal Analysis Logs'];
mkdir(['../Logs/',folder_name,'/']);
diary(['../Logs/',folder_name,'/Code_LogsLinearDynamic.txt'])
%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      INPUT/INIT SECTION      %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choice = 2;
diameter_now = 10;
damp_now = 0.05;
%     while(~(choice == 1 || choice == 2))
%         choice = input('[1]. Provide new input values including wall dimensions and loads, or\n[2]. Use default values\nChoose [1/2]:    ');
%     end
if choice == 1
    %%% Wall Dimension 
    height = input('Enter the height of the Wall(in mm):    ');
    width = input('Enter the width of the Wall(in mm):    ');
    thickness = input('Enter the thickness of the Wall(in mm):    ');

    % For linear analysis take the following youngs modulus. % In Newton 
    % per milli-meter square.
    steel_E = input('Enter the Young''s modulus of steel:    ');
    pois_ratio = input('Enter the possion_ratio:    ');

    bar_dia = input('Enter the main steel bar diameter:    ');

    %     div_x = input('Enter the division in x direction:    ');
    %     div_y = input('Enter the division in y direction:    ');
    %     div_z = input('Enter the division in z direction:    ');

    vertical_spacing = input('Enter the verticle reinforcement spacing:    ');
    horz_spacing = input('Enter the horizontal reinforcement spacing:    ');
    vertical_dia = input('Enter the diameter of verticle reinforcement bars:    ');
    horz_dia = input('Enter the diameter of horzontal reinforcement bars:    ');
    condition = 'all_fixed';

    conc_grade = input('Enter the grade of concrete:    ');
    steel_grade = input('Enter the grade of steel:    ');

    %%% Loads
    % One point load acting at centre.
else
    height = 3000;
    width = 5000;
    thickness = 170;
    reinf_present = true;
%   pois_ratio = .3;
    
    condition = 'one_fixed';
    
    conc_grade = 25;
    conc_pois_ratio = 0.18;
    conc_yield_strain = 0.002;
    face_num = 6;
    conc_den = 2.4*10^(-9); % In Mg/mm^3
    
    if reinf_present 
        bar_dia = diameter_now;
        vertical_spacing = 250; % 250 mm spacing of verticle reinforcement.
        horz_spacing = 300; % 300 mm spacing of horizontal reinforcement.
        vertical_dia = diameter_now;
        horz_dia = diameter_now;
    end
    steel_grade = 415;
    steel_pois_ratio = 0.3;
    steel_den = 8.05*10^(-9);
    steel_E = 2 * 10^5;
    g = 9.807*1000; % In mm/s^2
    element = '20-noded';
end

conc_E = 5000 * sqrt(conc_grade);
conc_Et = conc_E / 10;
steel_Et = steel_E / 5;
steel_yield_strain = steel_grade / steel_E;
num_eig_val = 10;

%% Setting height and width as per aspect ratio
% Comment it out later
%     if aspect_rat <= 1
%         height = 3000;
% %         width = 5000;
%         width = ceil(height/aspect_rat);
%     else
%         width = 3000;
%         height = ceil(width*aspect_rat);
%     end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Change steel to concrete %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Just to test concrete alone.
%     steel_den = conc_den;
%     steel_E = conc_E;
%     steel_Et = conc_Et;
%     steel_yield_strain = conc_yield_strain;
%     steel_pois_ratio = conc_pois_ratio;
    %%%% UPTO HERE ONLY

% Dimensions matrix is in mm so to avoid float values.
dimension = [thickness, width, height];
%%% Naive assumptions
% It seems that 15 divisions across y and z direction are enough. So
% creating the divisions accordingly.
div_y = 10;
div_z = 6;
div_x = 1; % This should vary according to the bars positioning.
divisions = [div_x, div_y, div_z];

if reinf_present
    vert_side_cover = 75;% Veticle bars side cover.
    horz_side_cover = 75;
    reinforcment_info = [vertical_dia, vertical_spacing, vert_side_cover;
                         horz_dia    , horz_spacing    , horz_side_cover];
end

%% Create mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%        MESH CREATION        %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Creating Mesh...');
if reinf_present
    [nodal_connect, nodal_coordinate, faces, mesh_meta_data, bar_position] = createMesh(dimension, divisions, reinforcment_info);
else
    [nodal_connect, nodal_coordinate, faces, mesh_meta_data] = createMeshWithoutReinf(dimension, divisions);
end
if strcmp(element, '20-noded')
    [extended_nodal_connect, total_no_nodes, face_nodes, hex_equ_face_nodes, hex_equ_nodes] = getExtendedNodalConnect(mesh_meta_data);
end
disp('Done!');

%% Print Input
fprintf('Wall dimesnsions:\n\tHeight: %d,\tWidth: %d,\tThickness: %d.\n', height, width, thickness);
fprintf('Concrete properties:\n\tGrade: %d,\n\tPoisson''s ratio: %0.2f,\n\tDensity: %0.3e,\n\tYield stress: %0.3e,\n\tYield Strain: %0.3e\n',conc_grade, conc_pois_ratio, conc_den, conc_E, conc_yield_strain);
if reinf_present
    fprintf('Steel properties:\n\tGrade: %d,\n\tPoisson''s ratio: %0.2f,\n\tDensity: %0.3e,\n\tYield stress: %0.3e,\n\tYield Strain: %0.3e\n',steel_grade, steel_pois_ratio, steel_den, steel_E, steel_yield_strain);
    fprintf('Verticle Bar specs:\n\tDiameter: %0.2f,\tSpacing: %0.2f,\t Side Cover: %0.2f.\n', vertical_dia, vertical_spacing, vert_side_cover);
    fprintf('Horizontal Bar specs:\n\tDiameter: %0.2f,\tSpacing: %0.2f,\t Side Cover: %0.2f.\n', horz_dia, horz_spacing, horz_side_cover);
end
fprintf('Boundary condition: %s\n', condition);
fprintf('Damping: %0.2f\n', damp_now);
fprintf('Element used: %s\n', element);
fprintf('Mesh Meta Data: \n\tNumber of divisions in x-dir: %d, \n\tNumber of divisions in y-dir: %d, \n\tNumber of divisions in z-dir: %d\n', mesh_meta_data(1), mesh_meta_data(2), mesh_meta_data(3));


%% Draw the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%        MESH 3D PLOT        %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     disp('Plotting Mesh...');
%     draw3DMesh(nodal_coordinate, faces, zeros(length(nodal_coordinate),
%     1));
%     disp('Done!');


%% Get element type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%          ASSIGN ELEMENTS TYPES          %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if reinf_present
    disp('Getting each element type...');
    element_type_steel = getElementType(nodal_coordinate, nodal_connect, bar_position, thickness);
    disp('Done!');
else
    no_elements = size(nodal_connect);
    no_elements = no_elements(1);
    element_type_steel = zeros(no_elements, 1);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%             HELPER INIT             %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total number of elements will be equal to the the size of the
% nodal_coordinate matrix.
no_elements = size(nodal_connect);
no_elements = no_elements(1);
if strcmp(element, '8-noded')
    total_no_nodes = length(nodal_coordinate);
    element_mapping = ElementMapping(nodal_connect, no_elements, element);
elseif strcmp(element, '20-noded')
    element_mapping = ElementMapping(extended_nodal_connect, no_elements, element);
end
% Initialize the E vector for each element with elastic modulus of
% elasticity.
element_mod_of_elas = steel_E.*element_type_steel(1:no_elements) + conc_E.*(~element_type_steel(1:no_elements));
element_pois_ratio = steel_pois_ratio.*element_type_steel(1:no_elements) + conc_pois_ratio.*(~element_type_steel(1:no_elements));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%            PERCENTAGE OF STEEL          %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_steel_vol = 0;
for i = 1:no_elements
    if element_type_steel(i)
        dx = nodal_coordinate(nodal_connect(i, 2), 1) - nodal_coordinate(nodal_connect(i, 1), 1);
        dy = nodal_coordinate(nodal_connect(i, 3), 2) - nodal_coordinate(nodal_connect(i, 2), 2);
        dz = nodal_coordinate(nodal_connect(i, 5), 3) - nodal_coordinate(nodal_connect(i, 1), 3);
        total_steel_vol = total_steel_vol + dx*dy*dz;
    end
end
total_vol = width*height*thickness;
percentage_steel = total_steel_vol/total_vol*100;
fprintf('Percentage of steel: %0.6f\n', percentage_steel);
%% Mass matrix calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%          MASS MATRIX CALCULATION          %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% density = 1;
density =  steel_den.*element_type_steel(1:no_elements) + conc_den.*(~element_type_steel(1:no_elements));
[mass, global_mass] = getMassMat(nodal_coordinate, nodal_connect, density, element_mapping, element, no_elements, total_no_nodes);

% The created mass matrix should be close to symmetric. Make it symmetric
% and handle the exception case.
if(max(max(abs(global_mass - global_mass.'))) < 1e-5)
    global_mass = (global_mass + global_mass.')/2;
else
    disp('Variations in mass matrix calculation crossed acceptable error. Aborting...');
    return;
end


%% Initial Stiffness matrix calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      INITIAL STIFFNESS MATRIX      %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Since this is static analysis thus the global stiffness matrix will not
% vary as per time. Calculate stiffness matrix once.
[global_stiff, stiff] = getGlobalStiff(nodal_coordinate, nodal_connect, element_mod_of_elas, element_pois_ratio, element, element_mapping, total_no_nodes, no_elements);

%% Force application
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%        EXTERNAL FORCES       %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Assembling force...');
tic
%% EARTHQUAKE FORCE SECTION
force_type = 'IMP';
if strcmp(force_type, 'EQ')
    EQ_name = 'Elcentro.txt';
    disp(['Using ', EQ_name]);
    P = load(EQ_name);
    P(:, 2) = (P(:,2)*g);
    time_vec = P(:, 1);
    force_time_history = P(:,2);
    num_time_steps = length(force_time_history);
    dummy = repmat([1; 0; 0], total_no_nodes, 1);
    mass_col = global_mass*dummy;
    for i = 1:length(P)
        load_vec_time_history(:, :, i) = force_time_history(i)*mass_col;
    end
    time_step = 0.02;
    toc
elseif strcmp(force_type, 'IMP') % Stands for impulsive surface pressure load
    %%
    pressure_time_history = [     0, 0.017345;% Time in sec.
                             0.04413,       0;]; % Pressure pulse for 125kg TNT kept at 20m.
    points = 0;
    temp = repmat([pressure_time_history(1, 2); 0], 1, points);
    for i = 1:points
        temp(1, i) = temp(1, i)*(i+1);
    end
    pressure_time_history = [pressure_time_history temp];
    time_vec = pressure_time_history(1, :);
    plot(pressure_time_history(1, :), pressure_time_history(2, :), 'LineWidth', 1.9);
    save(['../Logs/',folder_name,'/pressure_time_hist.mat'], 'pressure_time_history');
    close;
    [face_indices, face_area] = getFaceIndices(mesh_meta_data, face_num, dimension);
    num_time_steps = length(pressure_time_history(1, :));
    time_step = pressure_time_history(1, 2) - pressure_time_history(1, 1); % Assuming constant time step for now.
    load_vec = zeros(total_no_nodes*3, 1, num_time_steps);
    for ii = 1:num_time_steps
        net_force = face_area*pressure_time_history(2, ii);
        load_vec(face_indices(1, :), :, ii) = net_force/length(face_indices(1, :));
    end
    load_vec_time_history = load_vec;
    %%
end
% load_vec = [density*g; 0; 0];
% force_type = 'body';
% tic
% nodal_force = getNodalForce(force_type, no_elements, load_vec);
% toc
% tic
% global_force = global_force_vector(nodal_coordinate, nodal_connect, nodal_force, element_mapping);
% toc

%%
disp('Applying boundary conditions...');
tic;
if strcmp(element, '8-noded')
    [global_stiff_bc, load_bc, global_mass_bc, boundary_pt_index, ~, ~] = boundary_conditions(condition, global_stiff, mesh_meta_data, load_vec_time_history, global_mass, element, {}, {});
elseif strcmp(element, '20-noded')
    [global_stiff_bc, load_bc, global_mass_bc, boundary_pt_index, reduced_index, hex_equ_boundary_index] = boundary_conditions(condition, global_stiff, mesh_meta_data, load_vec_time_history, global_mass, element, face_nodes, hex_equ_face_nodes);
end
all_indices = 1:total_no_nodes*3;
non_boundary_indices = setdiff(all_indices, boundary_pt_index);
toc;
disp('Done!');

%% Rayleigh damping matrix calculation
disp('Damping matrix calculations...');
tic
zeta = damp_now;
return;
% smallest_eig_value = eigs(global_stiff_bc, global_mass_bc, 1, 'smallestabs');
[eig_vecs, smallest_eig_vals] = eigs(global_stiff_bc, global_mass_bc, num_eig_val, 'smallestabs', 'Display', true);


%%
% 
% counter_1 = 1;
% displ_mesh = zeros(mesh_meta_data(3)+1, mesh_meta_data(2)+1);
% final_disp = zeros(total_dof, 1);
% % Setting boundary point displacement to be zero
% final_disp(boundary_pt_index) = 0; 
% final_disp(non_boundary_indices) = eig_vecs(:, 2);
% for ii = 1:mesh_meta_data(3)+ 1
%     for jj = 1:mesh_meta_data(2)+1
%         displ_mesh(ii, jj) = final_disp(counter_1);
%         counter_1 = counter_1 + 3;
%     end
% end
% distinct_z = unique(nodal_coordinate(:, 3), 'rows');
% distinct_y = unique(nodal_coordinate(:, 2), 'rows');
% 
% figure;
% contourf(distinct_y, distinct_z, displ_mesh);
% colorbar;
% savefig(['../Logs/',folder_name,'/2D displacement contour.fig']);
% % @TODO Comment/Remove
% % close;
% figure;
% surf(distinct_y, distinct_z, displ_mesh);
% colorbar;
% savefig(['../Logs/',folder_name,'/3D displacement contour.fig']);
% % @TODO Comment/Remove
% % close;

%% Participation factor
norm_eig_vecs = eig_vecs;
part_fac = [];
for mode = 1:num_eig_val
    phi = norm_eig_vecs(1:3:end, mode);
%     phi = phi/norm(phi);
    part_fac(end+1) = (phi.'*global_mass_bc(1:3:end, 1:3:end)*ones(length(eig_vecs)/3, 1)) / ((phi.')*global_mass_bc(1:3:end, 1:3:end)*phi);
end
n = vpa(abs(part_fac(1:num_eig_val))./sum(abs(part_fac(1:num_eig_val)))*100);
%%
for mode = 1:num_eig_val
    % To find out overall displacement including the removed boundary points
    % Setting boundary point displacement to be zero
    if strcmp(element, '8-noded')
        final_disp = zeros(total_no_nodes*3, 1);
        final_disp(boundary_pt_index) = 0;
        final_disp(non_boundary_indices) = eig_vecs(:, mode);
    elseif strcmp(element, '20-noded')
        final_disp = zeros(length(nodal_coordinate)*3, 1);
        final_disp(setdiff(1:length(nodal_coordinate)*3, hex_equ_boundary_index)) = eig_vecs(arrayfun(@(x)find(non_boundary_indices==x,1),setdiff(hex_equ_nodes, reduced_index)), mode);
        % final_disp = final_disp(setdiff(hex_equ_nodes, reduced_index));
    end

    nodal_delta_x = final_disp(1:3:length(final_disp));
    nodal_delta_y = final_disp(2:3:length(final_disp));
    nodal_delta_z = final_disp(3:3:length(final_disp));
    % Here multiplying by a factor just to visualize the deformation.
    new_nodal_coord = 1000*[nodal_delta_x nodal_delta_y nodal_delta_z] + nodal_coordinate;
    draw3DMesh(new_nodal_coord, faces, final_disp(1:3:end));
    savefig(['../Logs/',folder_name,'/Final deflected 3D mode shape for mode number- ', num2str(mode),'.fig']);
    % @TODO Comment/Remove
    close;
end

%% Saving video
v = VideoWriter(['../Logs/',folder_name,'/ModalAnalysisDeformation.avi']);
v.FrameRate = 0.5;
open(v);
h = draw3DMesh(nodal_coordinate, faces, zeros(length(nodal_coordinate), 1));
filename = ['../Logs/',folder_name,'/ModalAnalysisDeformation.gif'];
for mode = 1:num_eig_val
    % To find out overall displacement including the removed boundary points
    % Setting boundary point displacement to be zero
    if strcmp(element, '8-noded')
        final_disp = zeros(total_no_nodes*3, 1);
        final_disp(boundary_pt_index) = 0;
        final_disp(non_boundary_indices) = eig_vecs(:, mode);
    elseif strcmp(element, '20-noded')
        final_disp = zeros(length(nodal_coordinate)*3, 1);
        final_disp(setdiff(1:length(nodal_coordinate)*3, hex_equ_boundary_index)) = eig_vecs(arrayfun(@(x)find(non_boundary_indices==x,1),setdiff(hex_equ_nodes, reduced_index)), mode);
        % final_disp = final_disp(setdiff(hex_equ_nodes, reduced_index));
    end

    nodal_delta_x = final_disp(1:3:length(final_disp));
    nodal_delta_y = final_disp(2:3:length(final_disp));
    nodal_delta_z = final_disp(3:3:length(final_disp));
    % Here multiplying by a factor just to visualize the deformation.
    new_nodal_coord = 1000*[nodal_delta_x nodal_delta_y nodal_delta_z] + nodal_coordinate;
%     draw3DMesh(new_nodal_coord, faces, final_disp(1:3:end));
%     savefig(['../Logs/',folder_name,'/Final deflected 3D mode shape for mode number- ', num2str(mode),'.fig']);
%     % @TODO Comment/Remove
%     close;
    set(h, 'Vertices',new_nodal_coord, 'FaceVertexCData', final_disp(1:3:end));
    drawnow;
%     savefig(['../Logs/',folder_name,'/Final deflected 3D mode shape for mode number- ', num2str(mode),'.fig']);

    % Capture the plot as an image 
%     animation_frames(mode) = getframe(gca);
    frame = getframe(gca); 
    animation_frames(mode) = frame;
    writeVideo(v,frame);
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 

    % Write to the GIF File 
    if mode == 1
      imwrite(imind,cm,filename,'gif', 'Loopcount', inf, 'Delay', 1); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append', 'Delay', 1); 
    end
end

close(v);
diary off;