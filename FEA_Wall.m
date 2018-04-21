%% Anchal Pandey 
% 201316114
% Conventions used - 
    
% Using hexahedron element Finite Element Analysis of a Wall.
% Working and tested on MATLAB R2017b
% Dependancy: Neural Network Toolbox.
var_name = 'Damping ratio';
var_start = 0.0;
var_diff = 0.01;
var_end = 0.3;
disp_th_per_var = [];
max_disp_per_var = [];
var_vec = [];
now_index = 1;
% damp_now = 0.05;
diameter_now = 10;
for var_val = var_start:var_diff:var_end
    damp_now = var_val;
%     clear variables;
%     clear global;
    clearvars -except diameter_now damp_now var_name var_start var_diff var_end disp_th_per_var max_disp_per_var var_vec now_index var_val
    clc;
    mkdir('../Logs');
    folder_name = ['[',datestr(datetime, 'dd-mmm-yyyy, HH.MM AM'), '] Logs'];
    mkdir(['../Logs/',folder_name,'/']);
    diary(['../Logs/',folder_name,'/Code_LogsLinearDynamic.txt'])

    %% Input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%      INPUT/INIT SECTION      %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    choice = 2;
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
        thickness = 220;
        steel_E = 2 * 10^5;
        pois_ratio = .3;
        bar_dia = diameter_now; % 12mm diameter bars
        condition = 'one_fixed';
        vertical_spacing = 250; % 250 mm spacing of verticle reinforcement.
        horz_spacing = 300; % 300 mm spacing of horizontal reinforcement.
        vertical_dia = diameter_now;
        horz_dia = diameter_now;
        conc_grade = 25;
        steel_grade = 415;
        conc_pois_ratio = 0.18;
        steel_pois_ratio = 0.3;
        conc_yield_strain = 0.002;
        face_num = 6;
        conc_den = 2.49*10^(-9); % In Mg/mm^3
        steel_den = 8.05*10^(-9);
        g = 9.81*1000; % In mm/s^2
%         damp_now = damp_now; % 
    end
    conc_E = 5000 * sqrt(conc_grade);
    steel_Et = steel_E / 5;
    conc_Et = conc_E / 10;
    steel_yield_strain = steel_grade / steel_E;
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
    div_y = 15;
    div_z = 15;
    div_x = 1; % This should vary according to the bars positioning.
    divisions = [div_x, div_y, div_z];
    vert_side_cover = 75;% Veticle bars side cover.
    horz_side_cover = 75;
    reinforcment_info = [vertical_dia, vertical_spacing, vert_side_cover;
                         horz_dia    , horz_spacing    , horz_side_cover];


    %% Print Input
    fprintf('Wall dimesnsions:\n\tHeight: %d,\tWidth: %d,\tThickness: %d.\n', height, width, thickness);
    fprintf('Concrete properties:\n\tGrade: %d,\n\tPoisson''s ratio: %0.2f,\n\tDensity: %0.3e,\n\tYield stress: %0.3e,\n\tYield Strain: %0.3e\n',conc_grade, conc_pois_ratio, conc_den, conc_E, conc_yield_strain);
    fprintf('Steel properties:\n\tGrade: %d,\n\tPoisson''s ratio: %0.2f,\n\tDensity: %0.3e,\n\tYield stress: %0.3e,\n\tYield Strain: %0.3e\n',steel_grade, steel_pois_ratio, steel_den, steel_E, steel_yield_strain);
    fprintf('Boundary condition: %s\n', condition);
    fprintf('Verticle Bar specs:\n\tDiameter: %0.2f,\tSpacing: %0.2f,\t Side Cover: %0.2f.\n', vertical_dia, vertical_spacing, vert_side_cover);
    fprintf('Horizontal Bar specs:\n\tDiameter: %0.2f,\tSpacing: %0.2f,\t Side Cover: %0.2f.\n', horz_dia, horz_spacing, horz_side_cover);
    fprintf('Damping: %0.2f\n', damp_now);

    %% Create mesh
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%        MESH CREATION        %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Creating Mesh...');
    [nodal_connect, nodal_coordinate, faces, mesh_meta_data, bar_position] = createMesh(dimension, divisions, reinforcment_info);
    disp('Done!');


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
    disp('Getting each element type...');
    element_type_steel = getElementType(nodal_coordinate, nodal_connect, bar_position, thickness);
    disp('Done!');

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%             HELPER INIT             %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Total number of elements will be equal to the the size of the
    % nodal_coordinate matrix.
    no_elements = length(nodal_connect);
    total_no_nodes = length(nodal_coordinate);
    element_mapping = ElementMapping(nodal_connect, no_elements);
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
    [mass, global_mass] = getMassMat(nodal_coordinate, nodal_connect, density, element_mapping);
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
    global_stiff = getGlobalStiff(nodal_coordinate, nodal_connect, element_mod_of_elas, element_pois_ratio);


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
    [global_stiff_bc, load_bc, global_mass_bc, boundary_pt_index] = boundary_conditions(condition, global_stiff, mesh_meta_data, load_vec_time_history, global_mass);
    all_indices = 1:total_no_nodes*3;
    non_boundary_indices = setdiff(all_indices, boundary_pt_index);
    toc;
    disp('Done!');

    %% Rayleigh damping matrix calculation
    disp('Damping matrix calculations...');
    tic
    zeta = damp_now;
    smallest_eig_values = eigs(global_stiff_bc, global_mass_bc, 20, 'smallestabs');
    larg_eig_value = smallest_eig_values(end);
    smallest_eig_value = smallest_eig_values(1);
%     [eig_vecs, smallest_eig_vals] = eigs(global_stiff_bc, global_mass_bc, 50, 'smallestabs', 'Display', true, 'SubspaceDimension', 100);
%     larg_eig_value = eigs(global_stiff_bc, global_mass_bc, 1, 'largestabs', 'Display', true, 'SubspaceDimension', 100);
%     larg_eig_value = larg_eig_value(1);
    omega_low = sqrt(smallest_eig_value);
    omega_high = sqrt(larg_eig_value);
    % Above mentioned values are calculated for wall - [200 5000 3000]
    beta = 2*zeta/(omega_low + omega_high);
    alpha = omega_low*omega_high*beta;
    global_damp_bc = alpha*global_mass_bc + beta*global_stiff_bc;
    toc
    disp('Done!');
    fprintf('Alpha: %0.14f and Beta: %0.14f\n', alpha, beta);
    %% Initialization for applying newmarks's method
    disp('Initialization for applying newmarks''s method...');
    total_dof = length(global_stiff_bc);
    nodal_disp = zeros(total_dof, 1, num_time_steps);
    nodal_vel = zeros(total_dof, 1, num_time_steps);
    nodal_acc = zeros(total_dof, 1, num_time_steps);
    if strcmp(force_type, 'EQ')
        nodal_acc(:, :, 1) = P(1, 2)*repmat([1; 0; 0], total_dof/3, 1);
    elseif strcmp(force_type, 'IMP')
        nodal_acc(:, :, 1) = global_mass_bc\load_bc(:, :, 1);
    end
    disp('Newmarks method preprocessing:');
    tic
    alpha = 0.25;
    delta = 0.5;

    a = [1/(alpha*(time_step)^2);
         delta/(alpha*time_step);
         1/(time_step*alpha);
         (1/(2*alpha)) - 1;
         (delta/alpha) - 1;
         (time_step/2)*((delta/alpha) - 2);
         (time_step*(1 - delta));
         delta*time_step;
    ];
    % Damping matrix is zero otherwise add a damping matrix term here.[ADDED]
    eff_stiff = global_stiff_bc + a(1)*global_mass_bc + a(2)*global_damp_bc;
    % The created effective stiffness matrix should be close to symmetric. Make 
    % it symmetric and handle the exception case.
    if(max(max(abs(eff_stiff - eff_stiff.'))) < 1e-5)
        eff_stiff= (eff_stiff+ eff_stiff.')/2;
    else
        disp('Variations in effective stiffness matrix calculation crossed acceptable error. Aborting...');
        return;
    end
    da = decomposition(eff_stiff);
    toc
    disp('Done!');
    tic
    %% Solution
    %*********************************************************************
    % Above calculated informations are not needed to be calculated again.
    %*********************************************************************
    h = waitbar(0,'Linear dynamic analysis going on...');
    for i = 1:num_time_steps
        waitbar(i/num_time_steps);
        [nodal_disp, nodal_vel, nodal_acc] = apply_newmarks(eff_stiff, global_mass_bc, load_bc(:, :, i), nodal_disp, nodal_vel, nodal_acc, time_step, i, da, global_damp_bc);
    end
    close(h);
    toc
    %% Plot
    max_displ = [];
    for i = 1:num_time_steps
        max_displ(end+1) = (nodal_disp(end-2, :, i));
    end
    disp_th_per_var(:, :, now_index) = max_displ;
    max_text = ['\leftarrow Max displacement = ',num2str((max(max_displ)/height)*100), '% of the wall height.'];
    max_index = find(max_displ == max(max_displ));
    max_disp_per_var(now_index) = max(max_displ);
    var_vec(now_index) = var_val;
    figure;
    plot(time_vec, max_displ, 'LineWidth', 1.9, 'DisplayName',[var_name, ': ', num2str(var_val)]);
    text(time_vec(max_index), max_displ(max_index), max_text, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Time (sec)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Displacement (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Displacement time history in x-direction for Elcentro GM(Applied in x-direction) - Last Node');
    savefig(['../Logs/',folder_name,'/Displ time history of top node.fig']);
    close;
    %% Plot
    max_now = 0;
    for j = 1:3:length(nodal_disp(:, :, 1))
        if(max_now <= max(abs(nodal_disp(j, :, i))))
            max_now = max(abs(nodal_disp(j, :, i)));
    %         max_displ_new = [];
    %         for i = 1:num_time_steps
    %             max_displ_new(end+1) = nodal_disp(j, :, i);
    %         end
    %         plot(0:num_time_steps-1, max_displ_new);
    %         hold on;
            max_index = j;
        end
    end
    %% Each level plot
    %*********************************************************************
    % Caculate first point's dof in x-direction on each of the levels while 
    % moving upwards on wall. These points are will give displacement in 
    % x-direction.
    %*********************************************************************
    first_points = zeros(mesh_meta_data(3), 1);
    first_points(1) = 1;
    for i = 2:mesh_meta_data(3)
        first_points(i) = first_points(i-1) + (mesh_meta_data(2)+1)*3;
    end
    max_acc_each_time_step = [];
    max_disp_each_time_step = [];
    for ii = 1:mesh_meta_data(3)
        max_displ = [];
        max_acc = [];
        for i = 1:num_time_steps+1
            max_displ(end+1) = (nodal_disp(first_points(ii), :, i));
            max_acc(end+1) = (nodal_acc(first_points(ii), :, i));
        end
        max_acc_each_time_step(end+1) = max(max_acc);
        max_disp_each_time_step(end+1) = max(max_displ);
    end
    % Distinct dz values
    distinct_dz_coordinates = unique(nodal_coordinate(:, 3));
    distinct_dz_coordinates = distinct_dz_coordinates(2:end);

    figure;
    plot(distinct_dz_coordinates, max_acc_each_time_step/9.8/1000, '-p', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'white', 'LineWidth', 1.5);
    xlabel('Height (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Acceleration (g)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Max acceleration at each level');
    savefig(['../Logs/',folder_name,'/Max acc at level.fig']);
    % @TODO Comment/Remove
    close;
    figure;
    plot(distinct_dz_coordinates, max_disp_each_time_step, '-s', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'white', 'LineWidth', 1.5);
    xlabel('Height (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Displacement (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Max displacement at each level');
    savefig(['../Logs/',folder_name,'/Max disp at level.fig']);
    % @TODO Comment/Remove
    close;
%     %%
%     h = draw3DMesh(nodal_coordinate, faces, zeros(length(nodal_coordinate), 1));
%     filename = ['../Logs/',folder_name,'/ResponseAnimatedLinear.gif'];
% 
%     animation_frames(num_time_steps) = struct('cdata',[],'colormap',[]);
%     for i = 1:num_time_steps+1
%         % To find out overall displacement including the removed boundary points
%         final_disp = zeros(total_dof, 1);
%         % Setting boundary point displacement to be zero
%         final_disp(boundary_pt_index) = 0; 
%         final_disp(non_boundary_indices) = nodal_disp(:, :, i);
%         nodal_delta_x = final_disp(1:3:length(final_disp));
%         nodal_delta_y = final_disp(2:3:length(final_disp));
%         nodal_delta_z = final_disp(3:3:length(final_disp));
%         % Here multiplying by a factor just to visualize the deformation[Fixed]
%         new_nodal_coord = [nodal_delta_x nodal_delta_y nodal_delta_z] + nodal_coordinate;
%     %     draw3DMesh(new_nodal_coord, faces, final_disp(1:3:end));
%         set(h, 'Vertices',new_nodal_coord);
%         if strcmp(force_type, 'EQ')
%             drawnow limitrate;
%         elseif strcmp(force_type, 'IMP')
%             drawnow; % Use for Impulsive loads as then we will have lesser frame
%     %     rate.
%         end
%               % Capture the plot as an image 
%         animation_frames(i) = getframe(gca);
%         frame = getframe(gca); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256); 
% 
%         % Write to the GIF File 
%         if i == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount', 0, 'Delay', 0.1); 
%         else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%         end
%     end
%     close;
    %%
    counter_1 = 1;
    displ_mesh = zeros(mesh_meta_data(3)+1, mesh_meta_data(2)+1);
    final_disp = zeros(total_dof, 1);
    % Setting boundary point displacement to be zero
    final_disp(boundary_pt_index) = 0; 
    final_disp(non_boundary_indices) = nodal_disp(:, :, end);
    for ii = 1:mesh_meta_data(3)+ 1
        for jj = 1:mesh_meta_data(2)+1
            displ_mesh(ii, jj) = final_disp(counter_1);
            counter_1 = counter_1 + 3;
        end
    end
    distinct_z = unique(nodal_coordinate(:, 3), 'rows');
    distinct_y = unique(nodal_coordinate(:, 2), 'rows');

    figure;
    contourf(distinct_y, distinct_z, displ_mesh);
    colorbar;
    savefig(['../Logs/',folder_name,'/2D displacement contour.fig']);
    % @TODO Comment/Remove
    close;
    figure;
    surf(distinct_y, distinct_z, displ_mesh);
    colorbar;
    savefig(['../Logs/',folder_name,'/3D displacement contour.fig']);
    % @TODO Comment/Remove
    close;
    %%
    % clear global_stiff_bc global_stiff global_mass global_mass_bc da eff_stiff
    diary off;
    % save(['../Logs/',folder_name,'/Workspace_Variables.mat']);


    %%
    % To find out overall displacement including the removed boundary points
    final_disp = zeros(total_dof, 1);
    % Setting boundary point displacement to be zero
    final_disp(boundary_pt_index) = 0; 
    final_disp(non_boundary_indices) = nodal_disp(:, :, end);
    nodal_delta_x = final_disp(1:3:length(final_disp));
    nodal_delta_y = final_disp(2:3:length(final_disp));
    nodal_delta_z = final_disp(3:3:length(final_disp));
    % Here multiplying by a factor just to visualize the deformation.
    new_nodal_coord = 10*[nodal_delta_x nodal_delta_y nodal_delta_z] + nodal_coordinate;
    draw3DMesh(new_nodal_coord, faces, final_disp(1:3:end));
    savefig(['../Logs/',folder_name,'/Final deflected 3D mesh.fig']);
    % @TODO Comment/Remove
    close;
    now_index = now_index + 1;
end
%%
figure(1001);
% var_name = 'Diameter';
% var_vec = 4:1:22;
plot(var_vec, max_disp_per_var, '-d','MarkerFaceColor', 'red', 'LineWidth', 1.5);
xlabel(var_name, 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Max Displacement (mm)', 'FontSize', 12, 'FontWeight', 'bold');
title(['Max displacement vs ', var_name]);
savefig(['../Logs/',folder_name,'/', var_name, '_vs_max_disp_plot.fig']);
% close;
a = [var_vec', max_disp_per_var'];
save(['../Logs/',folder_name,'/', var_name, '_vs_max_disp.mat'], 'a');
% How to get the matrix
% example = matfile(['../Logs/',folder_name,'/damp_vs_max_disp.mat'])
% example.a