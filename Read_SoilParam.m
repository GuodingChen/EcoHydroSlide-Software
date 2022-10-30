% get the soil information
% please read me !!!!!!
% -----------------------------------------------------
%                       unit define
% Soil parameters unit when read from the soil USDA code
% Cohesion: (kpa)
% Friction_angle: (degree)
% Porosity: Non
% K_s: (cm/day)
% Theta_r: (cm^3/cm^3)
% Thera_s: (cm^3/cm^3)
% Dry_Unit_Weight: kN/m^3
% water_weight = 9.8KN/m^3
% VG1980_param_alpha = kPa^(-1)
% VG1980_param_n = unitless
% ------------------------------------------------------

function [SoilPara]= Read_SoilParam(Soil_USDA_code)
      
      
            No_data = -9999;
            soil = Soil_USDA_code;
            switch soil 
                % Look-up Table for HWSD Soil Texture, soil is coded by 
                % USDA Soil Textural Classes code
                case 0 % No_Soil
                    Code_cohesion = No_data; 
                    Friction_angle = 0;
                    Porosity = 0; 
                    K_s = 0;
                    Theta_r = 0;
                    Thera_s = 0;
                    Theta_fc = 0;
                    Theta_pw = 0;
                    Dry_Unit_Weight = 0;
                    VG1980_param_alpha = 0;
                    VG1980_param_n = 0;
                case 1  %  clay (heavy)                                                       
                    Code_cohesion = 50;
                    Friction_angle = 35;
                    Porosity = 0.4;
                    K_s = 14.750;
                    Theta_r = 0.2;
                    Thera_s = 0.4;
                    Theta_fc = 0.36;
                    Theta_pw = 0.21;
                    Dry_Unit_Weight = 25;
                    VG1980_param_alpha = 0.013 * 100 / 9.8;
                    VG1980_param_n = 1.15;
                case 2  %  silty clay                                                     
                    Code_cohesion = 30;
                    Friction_angle = 18.5;
                    Porosity = 0.49;
                    K_s = 9.614;
                    Theta_r = 0.198;
                    Thera_s = 0.481;
                    Theta_fc = 0.36;
                    Theta_pw = 0.21;
                    Dry_Unit_Weight = 18;
                    VG1980_param_alpha = 0.012 * 100 / 9.8;
                    VG1980_param_n = 2;
                case 3 % clay                                                       
                    Code_cohesion = 40;
                    Friction_angle = 16.5;
                    Porosity = 0.47;
                    K_s = 14.750;
                    Theta_r = 0.226;
                    Thera_s = 0.459;
                    Theta_fc = 0.36;
                    Theta_pw = 0.21;
                    Dry_Unit_Weight = 19.5;
                    VG1980_param_alpha = 0.013 * 100 / 9.8;
                    VG1980_param_n = 1.15;
                case 4  % silty clay loam
                    Code_cohesion = 50;
                    Friction_angle = 16.5;
                    Porosity = 0.48;
                    K_s = 11.108;
                    Theta_r = 0.172;
                    Thera_s = 0.481;
                    Theta_fc = 0.34;
                    Theta_pw = 0.19;
                    Dry_Unit_Weight = 14;
                    VG1980_param_alpha = 0.010 * 100 / 9.8;
                    VG1980_param_n = 1.92;
                case 5 % clay loam
                    Code_cohesion = 35;
                    Friction_angle = 20;
                    Porosity = 0.46;
                    K_s = 7.064;
                    Theta_r = 0.171;
                    Thera_s = 0.438;
                    Theta_fc = 0.34;
                    Theta_pw = 0.21;
                    Dry_Unit_Weight = 14;
                    VG1980_param_alpha = 0.017 * 100 / 9.8;
                    VG1980_param_n = 1.35;
                case 6 % silt                                                        
                    Code_cohesion = 9;
                    Friction_angle = 26.5;
                    Porosity = 0.52;
                    K_s = 43.747;
                    Theta_r = 0.126;
                    Thera_s = 0.489;
                    Theta_fc = 0.32;
                    Theta_pw = 0.165;
                    Dry_Unit_Weight = 16.5;
                    VG1980_param_alpha = 0.017 * 100 / 9.8;
                    VG1980_param_n = 1.33;
                case 7  % silt loam 
                    Code_cohesion = 9;
                    Friction_angle = 24;
                    Porosity = 0.46;
                    K_s = 18.471;
                    Theta_r = 0.115;
                    Thera_s = 0.439;
                    Theta_fc = 0.3;
                    Theta_pw = 0.15;
                    Dry_Unit_Weight = 14;
                    VG1980_param_alpha = 0.017 * 100 / 9.8;
                    VG1980_param_n = 1.41;
                case 8 % sandy caly
                    Code_cohesion = 24.5;
                    Friction_angle = 22.5;
                    Porosity = 0.41;
                    K_s = 11.353;
                    Theta_r = 0.237;
                    Thera_s = 0.385;
                    Theta_fc = 0.31;
                    Theta_pw = 0.23;
                    Dry_Unit_Weight = 18.5;
                    VG1980_param_alpha = 0.028 * 100 / 9.8;
                    VG1980_param_n = 1.35;
                case 9 % loam
                    Code_cohesion = 10;
                    Friction_angle = 22.5;
                    Porosity = 0.43;
                    K_s = 13.339;
                    Theta_r = 0.126;
                    Thera_s = 0.403;
                    Theta_fc = 0.26;
                    Theta_pw = 0.12;
                    Dry_Unit_Weight = 13;
                    VG1980_param_alpha = 0.022 * 100 / 9.8;
                    VG1980_param_n = 1.52;
                case 10 % sandy clay loam
                    Code_cohesion = 29;
                    Friction_angle = 20;
                    Porosity = 0.39;
                    K_s = 13.231;
                    Theta_r = 0.169;
                    Thera_s = 0.382;
                    Theta_fc = 0.33;
                    Theta_pw = 0.175;
                    Dry_Unit_Weight = 15;
                    VG1980_param_alpha = 0.028 * 100 / 9.8;
                    VG1980_param_n = 1.64;
                case 11 % sandy loam
                    Code_cohesion = 6;
                    Friction_angle = 32;
                    Porosity = 0.4;
                    K_s = 37.450;
                    Theta_r = 0.087;
                    Thera_s = 0.387;
                    Theta_fc = 0.23;
                    Theta_pw = 0.1;
                    Dry_Unit_Weight = 15;
                    VG1980_param_alpha = 0.037 * 100 / 9.8;
                    
                    %VG1980_param_n = 3.70;
                    % modified
                    VG1980_param_n = 1.50;
                case 12  % loamy sand
                    Code_cohesion = 7.5;
                    Friction_angle = 28.5;
                    Porosity = 0.42;
                    K_s = 108.199;
                    Theta_r = 0.077;
                    Thera_s = 0.390;
                    Theta_fc = 0.14;
                    Theta_pw = 0.06;
                    Dry_Unit_Weight = 20.5;
                    VG1980_param_alpha = 0.042 * 100 / 9.8;
                    %VG1980_param_n = 3.13;
                    % modified
                    VG1980_param_n = 1.50;
                 case 13 % sand
                    Code_cohesion = 5;
                    Friction_angle = 40;
                    Porosity = 0.43;
                    K_s = 642.954;
                    Theta_r = 0.061;
                    Thera_s = 0.375;
                    Theta_fc = 0.12;
                    Theta_pw = 0.04;
                    Dry_Unit_Weight = 21;
                    VG1980_param_alpha = 0.040 * 100 / 9.8;
                    %VG1980_param_n = 6.67;
                    % modified
                    VG1980_param_n = 1.50;
                otherwise
                    Code_cohesion = No_data; 
                    Friction_angle = 0;
                    Porosity = 0; 
                    K_s = 0;
                    Theta_r = 0;
                    Thera_s = 0;
                    Theta_fc = 0;
                    Theta_pw = 0;
                    Dry_Unit_Weight = 0;
                    VG1980_param_alpha = 0;
                    VG1980_param_n = 0;
            end
            
            
            
            
            SoilPara.Code_cohesion = Code_cohesion;
            SoilPara.Friction_angle = Friction_angle;
            SoilPara.Porosity = Porosity;
            SoilPara.K_s = K_s;
            SoilPara.Theta_r = Theta_r;
            SoilPara.Thera_s = Thera_s;
            
            SoilPara.Theta_fc = Theta_fc;
            SoilPara.Theta_pw = Theta_pw;
            SoilPara.Dry_Unit_Weight = Dry_Unit_Weight;
            SoilPara.VG1980_param_alpha = VG1980_param_alpha;
            SoilPara.VG1980_param_n = VG1980_param_n;
            
end













