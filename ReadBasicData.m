function  [Hydrodata, LandslideData, DownscalingData] = ReadBasicData()

        %% .asc file info
        ncols = 598; 
        nrows = 650;
        FormatString=repmat('%f',1,ncols);
%         NoData = -9999;
%         cellSize = 0.000833;
%         XLLCorner = 108.335815;
%         YLLCorner = 32.654663;


        %%  read the basical hydrological part file
        [HydroMask, R_hydro] = geotiffread("./Basics/DCatchment_Mask.tif");
        fid_slope = fopen("./Basics/DSlope.asc",'r');
        slope_hydro = cell2mat(textscan(fid_slope, FormatString, nrows, 'HeaderLines', 6));
        slope_hydro(HydroMask > 0) = slope_hydro(HydroMask > 0) / pi * 180;
        [aspect_hydro, ~] = geotiffread("./DownscalingBasicData/aspect_90.tif");
        [curvature_hydro, ~] = geotiffread("./DownscalingBasicData/curvature_coarse.tif");
        [TWI_hydro, ~] = geotiffread("./DownscalingBasicData/TWI_coarse.tif");

        %% read the landslide part file
        [LandMask, R_land] = geotiffread("./LandslideBasics/Dmask_fine.tif");
        [slope_land, ~] = geotiffread("./LandslideBasics/slope_fine.tif");
        [aspect_land, ~] = geotiffread("./LandslideBasics/aspect_fine.tif");
        [DEM_land, ~] = geotiffread("./LandslideBasics/DEM_fine.tif");
        [TWI_land, ~] = geotiffread("./DownscalingBasicData/TWI_fine.tif");
        [curvature_land, ~] = geotiffread("./DownscalingBasicData/curvature_fine.tif");
            

        fclose(fid_slope);
        Hydrodata.mask = HydroMask;
        Hydrodata.geoinfo = R_hydro;
        Hydrodata.slopeDeg = slope_hydro;
        DownscalingData.aspect = aspect_hydro;
        DownscalingData.curvature_coarse = curvature_hydro;
        DownscalingData.TWI_coarse = TWI_hydro;
        
        LandslideData.geoinfo = R_land;
        LandslideData.mask = LandMask;
        LandslideData.slope = slope_land;
        LandslideData.DEM = DEM_land;
        LandslideData.aspect = aspect_land;
        DownscalingData.curvature_fine = curvature_land;
        DownscalingData.TWI_fine = TWI_land;
        


end




