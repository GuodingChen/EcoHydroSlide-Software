function SM_12_5 = Soil_downscaling(SM_coarse, rundata, LandslideData, DownscalingData)
% read 90 m file
    % hydrological output data
    Nodata = -9999;
        mask = rundata.mask;
        Hydro_geoinfo = rundata.geoinfo;
        slope_90 = rundata.slopeDeg;  % degrees
        aspect_90 = DownscalingData.aspect;
        curvature_90 = DownscalingData.curvature_coarse;
        TWI_90 = DownscalingData.TWI_coarse;
      
        SM_90 = SM_coarse;
       
 
        % landslide data need
        Landslide_geoinfo = LandslideData.geoinfo;
        mask_fine = LandslideData.mask;
        %z_matrix = LandslideData.dem;
        slope_12_5 = LandslideData.slope;   % degrees
        aspect_12_5 = LandslideData.aspect;
        curvature_12_5 = DownscalingData.curvature_fine;
        TWI_12_5 = DownscalingData.TWI_fine;
        
        


    slope_90_mean = mean(slope_90(mask >0));
    
    [i_90, j_90] = find(mask>0);
    matrix_90 = [i_90, j_90];
    C_90 = mask - 1;
    C_90(mask<0) = Nodata;
    % read the 12.5 map
    
    slope_125_mean = mean(slope_12_5(mask_fine>0));
    [i_12_5, j_12_5] = find(mask_fine>0);
    Kw_12_5 = mask_fine - 1;
    SM_12_5 = mask_fine - 1;
    Kw_12_5(mask_fine<0) = Nodata;
    SM_12_5(mask_fine<0) = Nodata;
    matrix_12_5 = [i_12_5, j_12_5];
    
    SM_90_mean = mean(SM_90(mask>0));
    Kw_90 = SM_90 / SM_90_mean;
    Kw_90(mask<0) = Nodata;

        % prepare the regress parameter p and q
        for index_90 = 1 : length(matrix_90)
            i = matrix_90(index_90, 1);
            j = matrix_90(index_90, 2);
            A_90 = aspect_90(i,j);
            if (A_90 > 45) && (A_90 <= 135)
                k_e = 0.002;

            elseif (A_90 > 135) && (A_90 <= 225)
                k_e = 0.005;

            elseif (A_90 > 225) && (A_90 <= 315)
                k_e = -0.003;
            else
                k_e = -0.01;
            end
            ka_90 = (1 - k_e * slope_90(i,j)) / (1- k_e * slope_90_mean);

            % compute the c parameters
            if curvature_90(i,j) <= 0
                C_90(i,j) = Kw_90(i,j) / ka_90 - 0.16 * cosd(A_90) - 0.09 * sind(A_90);
            else
                C_90(i,j) = Kw_90(i,j) / ka_90 - 0.14 * cosd(A_90) - 0.10 * sind(A_90)...
                    + 0.02 * cosd(2 * A_90);
            end


        end

        fit_c_90 = C_90( C_90>0 );
        fit_TWI_90 = double(TWI_90( mask>0 ));

        line_regress = fit(fit_TWI_90, fit_c_90,  'Poly1', 'Robust', 'Bisquare');
        p = line_regress.p1;
        q = line_regress.p2;

        % compute the downscaled soil moisture

        for index_12_5 = 1 : length(matrix_12_5)
            i = matrix_12_5(index_12_5, 1);
            j = matrix_12_5(index_12_5, 2);
            A_12_5 = aspect_12_5(i,j);
            if (A_12_5 > 45) && (A_12_5 <= 135)
                k_e = 0.002;

            elseif (A_12_5 > 135) && (A_12_5 <= 225)
                k_e = 0.005;

            elseif (A_12_5 > 225) && (A_12_5 <= 315)
                k_e = -0.003;
            else
                k_e = -0.01;
            end
            ka_12_5 = (1 - k_e * slope_12_5(i,j)) / (1- k_e * slope_125_mean);

            % compute the c parameters
            if curvature_12_5(i,j) <= 0
                Kw_12_5(i,j) = ( p * TWI_12_5(i,j) + q + ...
                             0.16 * cosd(A_12_5) + 0.09 * sind(A_12_5) ) * ka_12_5;

            else
                Kw_12_5(i,j) = ( p * TWI_12_5(i,j) + q + ...
                             0.14 * cosd(A_12_5) + 0.10 * sind(A_12_5) -...
                             0.02 * cosd(2*A_12_5) ) * ka_12_5;

            end


        end
        % compute the soil moisture in finer map
        % coordinate information (90 m)
        x_min_90 = Hydro_geoinfo.LongitudeLimits(1);  
        x_max_90 = Hydro_geoinfo.LongitudeLimits(2);  
        y_min_90 = Hydro_geoinfo.LatitudeLimits(1);
        y_max_90 = Hydro_geoinfo.LatitudeLimits(2);
        cell_size_90 = Hydro_geoinfo.CellExtentInLatitude;
        % coordinate information (12.5 m)
        x_min_12_5 = Landslide_geoinfo.LongitudeLimits(1);  
        x_max_12_5 = Landslide_geoinfo.LongitudeLimits(2);  
        y_min_12_5 = Landslide_geoinfo.LatitudeLimits(1);
        y_max_12_5 = Landslide_geoinfo.LatitudeLimits(2);
        cell_size_12_5 = Landslide_geoinfo.CellExtentInLatitude;

        %-----------------generate big coordinate (resolution = 12.5 m)
        x_big = (x_min_12_5 + cell_size_12_5/2) : cell_size_12_5 : (x_max_12_5 - cell_size_12_5/2);
        y_big = (y_max_12_5 - cell_size_12_5/2) : -cell_size_12_5 : (y_min_12_5 + cell_size_12_5/2); % 
        [x_big_all, y_big_all] = meshgrid(x_big, y_big);



        i = matrix_12_5(1, 1);
        j = matrix_12_5(1, 2);
        x_WGS = x_big_all(i,j);
        y_WGS = y_big_all(i,j);
            %-----------find the slope location: i,j, in 90m map
        L_location_row = ceil( (y_max_90 - y_WGS) / cell_size_90 );
        L_location_column = ceil( (x_WGS - x_min_90) / cell_size_90 );


        row_initial = L_location_row;
        col_initial = L_location_column;
        Sum_SM = 0;
        n = 0;
        SM_compare = SM_90(L_location_row, L_location_column);

        for index_12_5 = 1 : length(matrix_12_5)

            i = matrix_12_5(index_12_5, 1);
            j = matrix_12_5(index_12_5, 2);
            x_WGS = x_big_all(i,j);
            y_WGS = y_big_all(i,j);
                %-----------find the slope location: i,j, in 90m map
            L_location_row = ceil( (y_max_90 - y_WGS) / cell_size_90 );
            L_location_column = ceil( (x_WGS - x_min_90) / cell_size_90 );

            SM_12_5(i,j) = Kw_12_5(i,j) * SM_90(L_location_row, L_location_column);


            if (row_initial ~= L_location_row) || (col_initial ~= L_location_column)
                delta = Sum_SM / n - SM_compare;

                row_initial = L_location_row;
                col_initial = L_location_column;
                Sum_SM = 0;
                n = 0;
                SM_compare = SM_90(L_location_row, L_location_column);
            end
            Sum_SM = Sum_SM + SM_12_5(i, j);
            n = n + 1;
        end
        % output the downscaling SM file with the CRS = R_fine

        % Output_FileName = strcat('SM_12_5_', dateStr_array(i_FileName));
        SM_12_5(SM_12_5<0) = -9999;
        
        %geotiffwrite("SM_12_5", SM_12_5, Landslide_geoinfo);  
        %disp(dateStr_array(i_FileName));
        
end















