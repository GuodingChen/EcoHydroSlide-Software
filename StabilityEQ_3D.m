function [FS_3D, Failure_volumn, Failure_area] = StabilityEQ_3D(a_e, b_e, c_e , z_center, x_cell, y_cell, SoilPara, ...
                landslide_aspect, x_transition, y_transition, main_slope, main_aspect, ...
                z_single_all, ellipse_location, ellipse_matrix, z_grid, judge_location,...
                ellipse_i, ellipse_j, SM_InEllipse, root_cohesion)
    
    % get the basic soil parameter
    
    SoilCohesion = SoilPara.Code_cohesion;
    Friction_angle = (SoilPara.Friction_angle / 180) * pi;
    Porosity = SoilPara.Porosity;
    
    Theta_r = SoilPara.Theta_r;
    Thera_s = SoilPara.Thera_s;

    
    Dry_Unit_Weight = SoilPara.Dry_Unit_Weight;
    VG1980_param_alpha = SoilPara.VG1980_param_alpha;
    VG1980_param_n = SoilPara.VG1980_param_n;               
    gamma_w = 9.81;                
    % relationship between root depth and landslide depth
    root_depth = 1.5;  % m
    if c_e > root_depth
        
        c3 = max(root_cohesion) * (root_depth / c_e);
    else
        
        c3 = max(root_cohesion);
    end
 
    %-----------------get the coordinate in ellipsoid system---------------
    x_ellipsoid = x_transition / cos(main_slope);
    y_ellipsoid = y_transition;
    if max(x_ellipsoid.^2/a_e^2+y_ellipsoid.^2/b_e^2) - 1 > 0.01
        error('check the a_e and delta_x')
    end
    delta_x=(-2*x_ellipsoid/(a_e^2)-...
        sqrt((2*x_ellipsoid/(a_e^2)).^2-4*(1/a_e^2+1/(c_e^2*(tan(main_slope))^2))*...
            (x_ellipsoid.^2/a_e^2+y_ellipsoid.^2/b_e^2-1)))/(2*(1/a_e^2+1/(c_e^2*(tan(main_slope))^2)));
    z_ellipsoid = delta_x / tan(main_slope);
    cell_number = length(ellipse_matrix);
    
    
    z_grid_InSlip = z_center + (z_ellipsoid-x_transition * sin(main_slope)) / cos(main_slope);
    z_single_all(ellipse_location) = z_grid_InSlip;
    
    
    D_raster = z_grid - z_grid_InSlip; 

%     if sum(sum(D_raster < 0)) > ceil(length(D_raster) / 2)
%         FS_3D = -9999;
%         Failure_volumn = -9999;
%         Failure_area = -9999;
%         return       
%     end
    
    % get the apparent dip and dip for each colume cell at slip surface 
    
      % z_grid_InSlip is the GIS coordinates at slip surface
    SlipSurface_aspect = single(zeros(cell_number,1));
%     SlipSurface_dip = single(zeros(cell_number,1));
%     aspect_result = single(zeros(cell_number,1));
  
    % SlipSurface_aspect = gpuArray.zeros(1,cell_number);
    % SlipSurface_dip = gpuArray.zeros(1,cell_number);
    z_a = z_single_all(sub2ind(size(judge_location), ellipse_i-1, ellipse_j-1));
    z_b = z_single_all(sub2ind(size(judge_location), ellipse_i-1, ellipse_j));
    z_c = z_single_all(sub2ind(size(judge_location), ellipse_i-1, ellipse_j+1));
    z_d = z_single_all(sub2ind(size(judge_location), ellipse_i, ellipse_j-1));
    z_f = z_single_all(sub2ind(size(judge_location), ellipse_i, ellipse_j+1));
    z_g = z_single_all(sub2ind(size(judge_location), ellipse_i+1, ellipse_j-1));
    z_h = z_single_all(sub2ind(size(judge_location), ellipse_i+1, ellipse_j));
    z_i = z_single_all(sub2ind(size(judge_location), ellipse_i+1, ellipse_j+1));
    dz_dx = ((z_c+2.*z_f+z_i)-(z_a+2.*z_d+z_g))/8;
    dz_dy = ((z_g+2.*z_h+z_i)-(z_a+2.*z_b+z_c))/8;
    aspect_result = atan2(dz_dy, -dz_dx);
    SlipSurface_aspect(aspect_result < 0) = pi / 2 - aspect_result(aspect_result < 0);
    SlipSurface_aspect(aspect_result > pi/2) = 2 * pi - aspect_result(aspect_result > pi/2) + pi / 2;
    SlipSurface_aspect(SlipSurface_aspect == 0) = pi / 2 - aspect_result(SlipSurface_aspect == 0);
    rise_run = sqrt((dz_dx ./ x_cell).^2 + (dz_dy ./ y_cell).^2);
    SlipSurface_dip = atan(rise_run);

    %Main_SlipAspect = mean(SlipSurface_aspect);
    Apparent_Dip_yz = atan(tan(SlipSurface_dip).* cos(SlipSurface_aspect));
    Apparent_Dip_xz = atan(tan(SlipSurface_dip).*sin(SlipSurface_aspect));
    SlipApparent_dip = atan(tan(SlipSurface_dip).*abs(cos(SlipSurface_aspect - main_aspect)));
    slipe_surface = x_cell * y_cell * sqrt(1 - (sin(Apparent_Dip_xz)).^2.*(sin(Apparent_Dip_yz)).^2)./...
    (cos(Apparent_Dip_xz).*cos(Apparent_Dip_yz));
    
    
    

    % 3D FS calculation 
    A = 80;          % ori = 20 KPa
    lambda = 0.7;        % ori = 0.4
    
    alpha = 3.4;      % ori = 3.4;
    %G_s = gamma_d / gamma_w;
    GroundWater_aspect = landslide_aspect;
    % compute the FS with finner sr
    
    Sr = SM_InEllipse( D_raster>0 ) / 100;
    %Sr = Sr - 0.1;
    m_t = Sr;
    %m_t = W_t ./ (c_e * n * (1-Sr));
    m_t(m_t>1) = 1;                  
    m_t(Sr==1) = 1;
    %V_raster = D_raster(D_raster>0) * x_cell * y_cell;
    c1 = SoilCohesion; % note cohesion used N/m^2 in readsoil.m
    %c2 = A .* Sr .* (1-Sr).^lambda .* (1-m_t).^alpha;
    c2 = A .*(1-lambda .* m_t.^alpha);
    c_final = c1 + c2 + c3;
    
    z_w = (abs(z_ellipsoid(D_raster>0)) - (c_e - m_t .* c_e .* cos(main_slope)))/cos(main_slope);


    z_w(z_w<0) = 0;

    G_soil = x_cell .* y_cell .* (Dry_Unit_Weight .* D_raster(D_raster>0) + gamma_w .* z_w .* (Porosity-1)+...
                gamma_w .* (D_raster(D_raster>0) - z_w) .* Porosity .* Sr);
    
    % prepare for seepage force---N_s & T_s
%     N_s = zeros(length(D_raster(D_raster>0)),1);
%     T_s = zeros(length(D_raster(D_raster>0)),1);
    %d_sat = m_t .* c_e;
    
    % S = gamma_w .*  d_sat(d_sat>0) .* x_cell .* y_cell .* sin(main_slope);
    S = gamma_w .*  z_w .* x_cell .* y_cell .* sin(main_slope);
    S_h = S .* cos(main_slope);
    S_v = S .* sin(main_slope);
    S_ch = S_h .* cos(GroundWater_aspect - SlipSurface_aspect(D_raster>0));
    S_mh = S_h .* cos(GroundWater_aspect - landslide_aspect);
    S_c = sqrt(S_v.^2 + S_ch.^2);
    S_m = sqrt(S_v.^2 + S_mh.^2);
    beta_Sc = acos(S_ch ./ S_c);
    beta_Sm = acos(S_mh ./ S_m);
    
    % compute the seepage force in two directions
    N_s = S_c .* sin(beta_Sc - SlipSurface_dip(D_raster>0));
    T_s = S_m .* cos(beta_Sm - SlipApparent_dip(D_raster>0));
    N_s(z_w == 0) = 0;
    T_s(z_w == 0) = 0;
    resist_force = (c_final .* slipe_surface(D_raster>0) + (G_soil .* cos(SlipSurface_dip(D_raster>0)) + N_s) .* tan(Friction_angle))...
                .* cos(SlipApparent_dip(D_raster>0));
    drive_force = (G_soil .* sin(SlipApparent_dip(D_raster>0)) + T_s) .* cos(SlipApparent_dip(D_raster>0));
    FS_3D = sum(resist_force) / sum(drive_force);
    % chech the FS value
    if FS_3D <= 0
        error('Invalid value of safety factor! Please check the program.')
    end
    Failure_volumn = sum( D_raster(D_raster > 0) * x_cell * y_cell);
    Failure_area = length( D_raster(D_raster > 0) ) * x_cell * y_cell;

end