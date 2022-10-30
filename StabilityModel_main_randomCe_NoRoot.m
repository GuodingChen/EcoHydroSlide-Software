clear;
clc;
tic
% read the root Biomass
 
H5fileName = "./LandslideBasics/daily&monthly&annual_data.h5";
H5infor = h5info(H5fileName);

% load the modelling moment time
data = load("./LandslideBasics/LandslideModelTime.mat");
ModelTime = data.final_maxtirx(:,2);
TimeIndex = data.final_maxtirx(:,3);

Vegetation_style = ["Bare","Normal","NoVegDynamics"]; % Bare, Normal, NoVegDynamics
% set the parallel pool
Parallelism_set = parpool(25);

for i_vege = 1 : 3%length(Vegetation_style)
Vegetation_name = Vegetation_style(i_vege);

    
switch Vegetation_name
    case "Bare"
        SMFile_path = "./Results/SM_bare/OutDT_SM_";
        OutputH5name = "./Results/FS&PF_percentage_bare";
        flag = 0;
    case "Normal"
        SMFile_path = "./Results/SM_DVeg/OutDT_SM_";
        OutputH5name = "./Results/FS&PF_percentage_DVeg";
        flag = 0;
        % read the Monthly_Broot
        Broot = h5read(H5fileName,'/MonthlyData/Monthly_Broot');
    case "NoVegDynamics"
        SMFile_path = "./Results/SM_SVeg/OutDT_SM_";
        OutputH5name = "./Results/FS&PF_percentage_SVeg";
        flag = 0;
        % read the MOYforMonthly_Broot
        Broot = h5read(H5fileName,'/DOY&MOYData/MOYforMonthly_Broot');
end

if flag == 1 || flag == 2
    
    Initial_Broot = mean(Broot);
end

% get the cohesion for vegetation
if flag == 1 || flag == 2
    [root_cohesion_initial] = RootCohesion_calculate();
else
    root_cohesion_initial = 0;
end
% read the basical data: hydro + land + downscaling
[Hydrodata, LandslideData, DownscalingData] = ReadBasicData();

% read the soil moisture file
ncols = 598; 
nrows = 650;
FormatString=repmat('%f',1,ncols);
Nodata_value = -9999;
  





% initial FS and PF statistical output
FS_PF_CountPercen = [];


Out_i = 1;
for i_SM = 1 : length(ModelTime)
    TimeMoment = ModelTime(i_SM);
    TimeMomentIndex = TimeIndex(i_SM);
    Monitor_variable = strcat(Vegetation_name,": ", TimeMoment);
    disp(Monitor_variable)
    % Number of SM files is 1212 (one per month)
    if flag == 1 || flag == 2
        Broot_moment = Broot(i_SM);
        root_cohesion = Broot_moment / Initial_Broot * root_cohesion_initial;
    else
        root_cohesion = 0;
    end
    
    
    
    SM_FileName = strcat(SMFile_path, TimeMoment, ".asc");
    fid_SM = fopen(SM_FileName,'r');
    SM_coarse = cell2mat(textscan(fid_SM, FormatString, nrows, 'HeaderLines', 6));
    fclose(fid_SM);
    % soil downscaling
    SM_fine = Soil_downscaling(SM_coarse, Hydrodata, LandslideData, DownscalingData);
    
    % set the soil parameters
    % case "sandy loam"
    SoilCohesion = 6;  % kPa
    Friction_angle = 32;  % angle: degree
    K_sat = 50;  % mm/h
    Theta_r = 0.041;
    Thera_s = 0.55;
    Dry_Unit_Weight = 15;  % kN/m^3
    VG1980_param_alpha = 0.5;
    VG1980_param_n = 1.35;
    Porosity = 0.4;
    
    % define the soil data for 3D model
    SoilPara.Code_cohesion = SoilCohesion;
    SoilPara.Porosity = Porosity;
    SoilPara.Theta_r = Theta_r;
    SoilPara.Thera_s = Thera_s;
    SoilPara.Dry_Unit_Weight = Dry_Unit_Weight;
    SoilPara.VG1980_param_alpha = VG1980_param_alpha;
    SoilPara.VG1980_param_n = VG1980_param_n;
    SoilPara.Friction_angle = Friction_angle;
    
    % landslide data need
        Landslide_geoinfo = LandslideData.geoinfo;
        % only part area of yuehe basin is considered
        % cell with soil moisture > 0 is valid 
        mask_fine = SM_fine;
        mask_fine(mask_fine>0) = 1;
        z_matrix = LandslideData.DEM;
        slope_matrix = LandslideData.slope;
        aspect_matrix = LandslideData.aspect;
        
        
        % landslide model basic setting
        % user-defined density and tile number
        ellipse_dendity = 20;
        total_tile_number = 120;
        
        % user-defined random ellipse geo restriction
        min_ae = 50; 
        max_ae = 100;
        min_be = 50; 
        max_be = 100;
        min_ce = 0.5;
        max_ce = 5;
        a_e_range = (min_ae : 5 : max_ae)';
        b_e_range = (min_be : 5 : max_be)';
        c_e_range = (min_ce : 0.1 : max_ce)';
        % data adjustment----------------
        z_matrix( z_matrix<0 ) = 0;
        aspect_matrix( aspect_matrix < 0 ) = 1;
        slope_matrix(slope_matrix == 0) = 1;
        all_slope = (slope_matrix / 180) * pi;
        all_aspect = (aspect_matrix / 180) * pi;
        
        % coordinate information (12.5 m)
        row_12_5 = Landslide_geoinfo.RasterSize(1);
        column_12_5 = Landslide_geoinfo.RasterSize(2);
        
        %-----------------generate the user-defined coordinate (resolution = 12.5 m)
        % let the coordinate transformation occurs within a Cartesian coordinate system
        % it will be More convenient
        x_cell = 12.5;
        y_cell = 12.5;
        x_min_set = 0;
        x_max_set = x_cell * (column_12_5-1);
        y_min_set = 0; 
        y_max_set = y_cell * (row_12_5-1);
        x = x_min_set : x_cell : x_max_set;
        y = y_max_set : -y_cell : y_min_set; % y coordinate should be carefully defined
        [x_all, y_all] = meshgrid(x,y);
        
        %-----% split the study area into sevaral tiles------
        [i_mask, j_mask] = find(mask_fine>0);
        mask_vector = [i_mask, j_mask];
        
        % calculate the cell number for each tile
        CellNumber_InEachTile = floor( length(mask_vector) / total_tile_number );
        % calculate the area ofe each tile (m^2)
        A_s = CellNumber_InEachTile * 12.5 * 12.5; % 
        
        % initial the whole FS_3D map
        % it is no matter what cell value in FS3D_whole_map due to that the it will
        % be replaced by each tile in final step.
        FS3D_whole_map = mask_fine;
        Volume_whole_map = mask_fine;
        
        Area_whole_map = mask_fine;
        PF_UnstableCount_map = mask_fine;
        PF_TotalCount_map = mask_fine;
        
        FS3D_whole_map(mask_fine > 0) = 2;
        PF_UnstableCount_map(mask_fine > 0) = 1;
        PF_TotalCount_map(mask_fine > 0) = 1;

        
        % transfer double matrix to single
        Volume_whole_map = single(Volume_whole_map);
        Area_whole_map = single(Area_whole_map);
        FS3D_whole_map = single(FS3D_whole_map);
        PF_UnstableCount_map = single(PF_UnstableCount_map);
        PF_TotalCount_map = single(PF_TotalCount_map);
        mask_vector = single(mask_vector);
        mask_fine = single(mask_fine);
        x_all = single(x_all); 
        y_all = single(y_all);
        z_all = single(z_matrix);
        SM_fine = single(SM_fine);
        
        
        % initial the cell type matrix
        FS3D_cell = cell(total_tile_number, 1);
        Volume_cell = cell(total_tile_number, 1);
        Area_cell = cell(total_tile_number, 1);
        PF_count_cell = cell(total_tile_number, 1);
        PF_unstable_cell = cell(total_tile_number, 1);
        RandomIndex_tile = cell(total_tile_number, 1);
        
        % here use parfor
        parfor tile_class = 1 : total_tile_number
            
            
            if tile_class == 1
                random_matrix = mask_vector(1:CellNumber_InEachTile,:);
            
            else
                random_matrix = mask_vector( ((tile_class-1) * CellNumber_InEachTile) : (tile_class*CellNumber_InEachTile),: );
            end
            
            FS3D_tile_matrix = ones(length(random_matrix), 1) * 1000;
            Volume_tile_matrix = zeros(length(random_matrix), 1);
            Area_tile_matrix = zeros(length(random_matrix), 1);
            PF_unstable_matrix = zeros(length(random_matrix), 1);
            PF_count_matrix = zeros(length(random_matrix), 1);        
            FS3D_tile_final = FS3D_tile_matrix;
            Volume_tile_final = Volume_tile_matrix;
            Area_tile_final = Area_tile_matrix;
            random_matrix_MinColumn = min(random_matrix(:,2));
            random_matrix_MaxColumn = max(random_matrix(:,2));
            % each tile corresponding to a random_matrix_LineIndex, and there
            % is impossible to get the intersection between any two matrix
            random_matrix_LineIndex = sub2ind(size(mask_fine), random_matrix(:,1), random_matrix(:,2));
            %disp(strcat('tile = ', num2str(tile_class)))
            ellipsoid_number = ellipse_dendity * 16 * A_s / (pi * (min_ae + max_ae) * ...
                (min_be + max_be)); 
            ellipsoid_number = ceil(ellipsoid_number);
            % begin to generate the random ellipse surface
            for num_ellipsoid = 1 : ellipsoid_number

                    %disp(['tile : ', num2str(tile_class), 'elli_N = ', num2str(num_ellipsoid)]);
                    %disp(strcat('ellipsoid number = ', num2str(num_ellipsoid)))

                    center_index = randi(length(random_matrix));
                    % (i, j) is the index in whole study region 
                    i = random_matrix(center_index, 1);
                    j = random_matrix(center_index, 2);

                %---------------------------------------------------  
 
                    % determine slope and aspect in 12m map 
                    slope = all_slope(i, j);
                    landslide_aspect = all_aspect(i, j);


                    %---determine the ellipse center in 12.5m map
                    x_center = x_all(i, j);
                    y_center = y_all(i, j);
                    z_center = z_all(i, j);
                    
                    % make a smaller matrix which contain a single ellipse
                    x_single_all = x_all((i-30) : (i+30), (j-30) : (j+30));
                    y_single_all = y_all((i-30) : (i+30), (j-30) : (j+30));
                    z_single_all = z_all((i-30) : (i+30), (j-30) : (j+30));
                    slope_single_all = all_slope((i-30) : (i+30), (j-30) : (j+30));
                    aspect_single_all = all_aspect((i-30) : (i+30), (j-30) : (j+30));
                    SM_fine_single = SM_fine((i-30) : (i+30), (j-30) : (j+30));
                    all_x_transition = (x_single_all-x_center) * cos(landslide_aspect) ...
                                           + (y_single_all-y_center) * sin(landslide_aspect);
                    all_y_transition = (y_single_all-y_center) * cos(landslide_aspect) ...
                                           - (x_single_all-x_center) * sin(landslide_aspect); 

                
                   
                    % random the length and width of ellipsoid
                    ae_random_index = randi([1,length(a_e_range)]);
                    be_random_index = randi([1,length(b_e_range)]);
                    ce_random_index = randi([1,length(c_e_range)]);
                    a_e_pre = a_e_range(ae_random_index);
                    b_e = b_e_range(be_random_index);
                    c_e = c_e_range(ce_random_index);
                    project_length = a_e_pre * cos(slope);
                    
                    %------find the involved ellipse---------
                    % jugge_location is the index in small single region which can 
                    % include a single ellipse
                    judge_location = (all_x_transition).^2/project_length^2+(all_y_transition).^2/b_e^2;
                    ellipse_location = find(judge_location<1);  % value<1 represent that the points are included in project ellipse
                    [ellipse_i, ellipse_j] = find(judge_location<1);
                    ellipse_matrix = [ellipse_i, ellipse_j];
                    ellipse_matrix = single(ellipse_matrix);
                  %---------------------------get Ellipsoidal region----------------
                    z_grid = z_single_all(ellipse_location);
                    x_transition = all_x_transition(ellipse_location);
                    y_transition = all_y_transition(ellipse_location);
                    slope_c = slope_single_all(ellipse_location);
                    aspect_c = aspect_single_all(ellipse_location);
                    % find soil moisture
                    SM_InEllipse = SM_fine_single(ellipse_location);
                    

                    main_slope = mean(slope_c,1);
                    main_aspect = mean(aspect_c,1);
                    a_e = project_length/cos(main_slope);
                    % ensure there are no nodata cell in SM matrix
                    
                    if isempty(find(SM_InEllipse<0, 1))

                    
                        [FSR_3D, Failure_volumn, Failure_area] = StabilityEQ_3D(a_e, b_e, c_e , z_center, ...
                        x_cell, y_cell, SoilPara, landslide_aspect, ...
                           x_transition, y_transition, main_slope, main_aspect, ...
                            z_single_all, ellipse_location, ellipse_matrix, z_grid, judge_location, ...
                            ellipse_i, ellipse_j, SM_InEllipse, root_cohesion);

                        % transform the ellipse_location index to (i,j) index map
                        [MovingWindow_center_i, ~] = find(y_single_all == y_center);
                        [~, MovingWindow_center_j] = find(x_single_all == x_center);
                        MovingWindow_center_i = unique(MovingWindow_center_i);
                        MovingWindow_center_j = unique(MovingWindow_center_j);

                        % update the ellipse_matrix into (i,j) index map (whole region)
                        ellipse_matrix(:,1) = ellipse_matrix(:,1) + (i - MovingWindow_center_i);
                        ellipse_matrix(:,2) = ellipse_matrix(:,2) + (j - MovingWindow_center_j);

                        % exclude the index which is out of the random_maxtrix,
                        % FS_ellipse_matrix represents the coordinate corresponding to
                        % the mask_fine (the total area)
                        %FS_ellipse_matrix = ellipse_matrix()
                        FS_column_index = find((ellipse_matrix(:,2) >= random_matrix_MinColumn)...
                            & (ellipse_matrix(:,2) <= random_matrix_MaxColumn));

                        FS_ellipse_matrix = ellipse_matrix(FS_column_index, :);
                        %FS_ellipse_matrix = intersect(random_matrix, ellipse_matrix, 'rows');

                        % FS_tile_index represents the normal index (begin with 1) in
                        % each tile of FS3D_tile_maxtrix
                        FS_ellipse_matrix_LineIndex = sub2ind(size(mask_fine), FS_ellipse_matrix(:,1), FS_ellipse_matrix(:,2));
                        %FS_tile_index = ismember(random_matrix, FS_ellipse_matrix, 'rows');

                        FS_tile_index = ismember(random_matrix_LineIndex, FS_ellipse_matrix_LineIndex);
                        FS3D_tile_matrix(FS_tile_index) = FSR_3D;



                        % count the cell considered into random test
                        PF_count_matrix(FS_tile_index) = PF_count_matrix(FS_tile_index) + 1;

                        % count the test which computed as unstable cell
                        if FSR_3D < 1
                            Volume_tile_matrix(FS_tile_index) = Failure_volumn;
                            Area_tile_matrix(FS_tile_index) = Failure_area;
                            PF_unstable_matrix(FS_tile_index) = PF_unstable_matrix(FS_tile_index) + 1;
                        end
                    end

                    % update the FS value

                    FS3D_tile_final = min(FS3D_tile_final, FS3D_tile_matrix);
                    Volume_tile_final = max(Volume_tile_final, Volume_tile_matrix);
                    Area_tile_final = max(Area_tile_final, Area_tile_matrix);
  
                    
            end
            % adjust the nodata in raster

            FS3D_tile_final(FS3D_tile_final == 1000) = 10;
            Volume_tile_final(Volume_tile_final == 0) = 0;
            Area_tile_final(Area_tile_final == 0) = 0;
            % storage each tile in tile class
            FS3D_cell{tile_class, 1} = FS3D_tile_final;
            Volume_cell{tile_class, 1} = Volume_tile_final;
            Area_cell{tile_class, 1} = Area_tile_final;
            RandomIndex_tile{tile_class, 1} = random_matrix_LineIndex;
            PF_count_cell{tile_class, 1} = PF_count_matrix;
            PF_unstable_cell{tile_class, 1} = PF_unstable_matrix;

        end
        
        disp('# # # # # # End of current time calculation # # # # # # ')
        % reconstruct the tiles to the whole FS/PF map
        % In theory, the cell number of FS ~= nan will equal to PF ~= nan, but
        % this program initial the PF_unstable_matrix as zero matrix, and dose
        % not consider the cell which is not randomly selected. Therefor, PF ~= nan 
        % will have more cell number than that of FS ~= nan
        for i = 1 : total_tile_number
            EachTile_index = RandomIndex_tile{i};
            FS3D_whole_map(EachTile_index) = FS3D_cell{i};
            PF_TotalCount_map(EachTile_index) = PF_count_cell{i};
            PF_UnstableCount_map(EachTile_index) = PF_unstable_cell{i};
            Volume_whole_map(EachTile_index) = Volume_cell{i};
            Area_whole_map(EachTile_index) = Area_cell{i};

        end
        PF_whole_map = PF_UnstableCount_map ./ PF_TotalCount_map;
        PF_whole_map(PF_TotalCount_map == 0) = 0;
        PF_whole_map(mask_fine ~= 1) = Nodata_value;
         % output the file in each observed time step
%         path_output = "../Results/";
%         FS_FileName = strcat(path_output, 'FS_3D_', ModelTime(i_SM));
%         PF_FileName = strcat(path_output, 'PF_', ModelTime(i_SM));
%         Volume_Filename = strcat(path_output, 'Volume_', ModelTime(i_SM));
%         Area_Filename = strcat(path_output, 'Area_', ModelTime(i_SM));
%         geotiffwrite(FS_FileName, FS3D_whole_map, Landslide_geoinfo);  
%         geotiffwrite(PF_FileName, PF_whole_map, Landslide_geoinfo);
%         geotiffwrite(Volume_Filename, Volume_whole_map, Landslide_geoinfo);
%         geotiffwrite(Area_Filename, Area_whole_map, Landslide_geoinfo);
        

        % count the necessary index of moment
        TOT_CellNumber = sum(sum(mask_fine == 1));
        
        
        NUMA_UnstableCell = sum(sum(FS3D_whole_map < 1 & FS3D_whole_map > 0));
        
        NUMA_PF10_cell = sum(sum(PF_whole_map > 0.1));
        NUMA_PF20_cell = sum(sum(PF_whole_map > 0.2));
        NUMA_PF30_cell = sum(sum(PF_whole_map > 0.3));
        NUMA_PF40_cell = sum(sum(PF_whole_map > 0.4));
        NUMA_PF50_cell = sum(sum(PF_whole_map > 0.5));
        
        Percen_Unstable = NUMA_UnstableCell / TOT_CellNumber;
        
        Percen_PF10 = NUMA_PF10_cell / TOT_CellNumber;
        Percen_PF20 = NUMA_PF20_cell / TOT_CellNumber;
        Percen_PF30 = NUMA_PF30_cell / TOT_CellNumber;
        Percen_PF40 = NUMA_PF40_cell / TOT_CellNumber;
        Percen_PF50 = NUMA_PF50_cell / TOT_CellNumber;
        
        FS_PF_CountPercen(Out_i, 1) = TimeMoment;
        FS_PF_CountPercen(Out_i, 2) = TimeMomentIndex;
        FS_PF_CountPercen(Out_i, 3) = Percen_Unstable;
        FS_PF_CountPercen(Out_i, 4) = Percen_PF10;
        FS_PF_CountPercen(Out_i, 5) = Percen_PF20;
        FS_PF_CountPercen(Out_i, 6) = Percen_PF30;
        FS_PF_CountPercen(Out_i, 7) = Percen_PF40;
        FS_PF_CountPercen(Out_i, 8) = Percen_PF50;
        
        % count the volume 
        FS_PF_CountPercen(Out_i, 9) = sum(Volume_whole_map(Volume_whole_map>0));
        FS_PF_CountPercen(Out_i, 10) = sum(Volume_whole_map(PF_whole_map>0.1));
        FS_PF_CountPercen(Out_i, 11) = sum(Volume_whole_map(PF_whole_map>0.2));
        FS_PF_CountPercen(Out_i, 12) = sum(Volume_whole_map(PF_whole_map>0.3));
        FS_PF_CountPercen(Out_i, 13) = sum(Volume_whole_map(PF_whole_map>0.4));
        FS_PF_CountPercen(Out_i, 14) = sum(Volume_whole_map(PF_whole_map>0.5));
        
        % count the area
        FS_PF_CountPercen(Out_i, 15) = sum(Area_whole_map(Area_whole_map>0));
        FS_PF_CountPercen(Out_i, 16) = sum(Area_whole_map(PF_whole_map>0.1));
        FS_PF_CountPercen(Out_i, 17) = sum(Area_whole_map(PF_whole_map>0.2));
        FS_PF_CountPercen(Out_i, 18) = sum(Area_whole_map(PF_whole_map>0.3));
        FS_PF_CountPercen(Out_i, 19) = sum(Area_whole_map(PF_whole_map>0.4));
        FS_PF_CountPercen(Out_i, 20) = sum(Area_whole_map(PF_whole_map>0.5));
        
        % hazard metric
        FS_PF_CountPercen(Out_i, 21) = sum(PF_whole_map(PF_whole_map>=0).*...
                                            Volume_whole_map(Volume_whole_map>=0));
        
        
        % Check if the database already exists
        OutputH5name_moment = strcat(OutputH5name, "RandomCe-100Year-NoRoot", ".h5");
        if ~exist(OutputH5name_moment,'file')

            h5create(OutputH5name_moment, '/FS_PF_CountPercen', size(FS_PF_CountPercen))
            h5write(OutputH5name_moment, '/FS_PF_CountPercen', FS_PF_CountPercen)
        else
            delete(OutputH5name_moment)
            h5create(OutputH5name_moment, '/FS_PF_CountPercen', size(FS_PF_CountPercen))
            h5write(OutputH5name_moment, '/FS_PF_CountPercen', FS_PF_CountPercen)
        end
        
        
        Out_i = Out_i + 1;
end



end



delete(Parallelism_set)    
disp('# # # # # # All time calculation ends # # # # # # ')
toc
