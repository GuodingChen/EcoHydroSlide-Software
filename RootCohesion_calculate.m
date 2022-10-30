% Main Parameters of Root System, from arnone 2016 (WRR)
% clear;
% clc;
function [root_cohesion] = RootCohesion_calculate()
    


            % Root depth (m)
            Zr = 1.5;

            % Root Collar Diameter (mm)
            Dmax = 21.5;

            % Total number of branches in the root system
            Ntot = 123;

            % Exponent of root tensile strength power function
            a = -0.72; 

            % Coeff. of root tensile strength power function (Mpa mm)
            A0 = 28.1;

            % Exponent of Young module power function
            b = -1;

            % Scaling factor of Young module power function (Mpa mm)
            E0 = 696;


            % Root tortuosity parameter 
            R = 0.5;

            % Exponent of Weibull survivor function
            W = 11;

            % Single root length (mm)
            L0 = 100;

            % Soil matric potential at which the stomatal closure (Mpa)

            Psi = -0.5;

            % Soil matric potential at which the stomatal closure (Mpa)
            Psi_w = -0.3;

            % parameters of distribution profile

            r = 11.015;

            p = 0.631382;

            % proportionality factor
            alpha_d = 1.15;

            % rooted soil area (m^2)
            RootSoil_area = 0.697695501;

            % Define the landslide failure depth (m)
            SoilDepth = (0.1 : 0.1 : 1.5)';

            % define the constant root number part
            ConstantRootNumber = 10;

            % PDF pamameters: Exponential distribution
            PDF_para =  0.2390;

        



    % calculate the total topological length, TL

    TL = Zr * 1000 / L0;

    % in ConstantDepth, the root structure will not change
    ConstantDepth = 200; % mm
    ConstantLevelNumber = ConstantDepth / L0;

    % discrete probability function: f(i|r,p)


    % 
    % f_i = ( gamma( RL + r ) ./ ( gamma(r) .* factorial(RL) ) ) ...
    %        * (p.^r) .* ((1 - p).^RL);


    % calculate the number of roots at each topological level, Ni, round to a
    % integer number
    Ntot_PDF = Ntot - ConstantRootNumber * ConstantLevelNumber;


    % calculate the diamter for each root level

    % calculate the diamter for each root level
    % calculate the diamter for each root level

    Diameter_RootLevel = zeros(TL, 1); %mm
    Ni = zeros(TL, 1);
    Diameter_RootLevel(1) = Dmax;  % unit: mm 
    % compute root number
    for i = 1 : TL 
        if i <= ConstantLevelNumber

            Ni(i) = ConstantRootNumber;  

        else

            Ni(i) = round((1 - PDF_para)^(i-1-ConstantLevelNumber) * PDF_para * Ntot_PDF);

        end
    end
    % compute root diameter

    for i = 2 : TL
        Diameter_RootLevel(i) = Diameter_RootLevel(i - 1) * sqrt( Ni(i - 1) / (alpha_d * Ni(i)));
    end   

    Diameter_RootLevel(Ni == 0) = 0;
    % compute the maximum root tension
    Tr = A0 * Diameter_RootLevel.^a;  % in Mpa
    Tr(Ni == 0) = 0;
    % Young's modulus 

    E = E0 * Diameter_RootLevel.^b; % in Mpa
    E(Ni == 0) = 0;
    % maximum displacement
    X_max = Tr * L0 ./ (R * E);  % in mm
    X_max(Ni == 0) = 0;
    % tensile force for a single root

    Fmax_SingleRoot = zeros(TL, 1); %mm
    
    for i_force = 1 : TL
        D_i = Diameter_RootLevel(i_force);
        
        Xmax_i = X_max(i_force);
        X = (0 : 0.1 : 30);
        Ei = E(i_force);
        Fcurve_SingleRoot = R*Ei*pi*(D_i^2) * X ./ (4*L0) .* exp(-(X/Xmax_i).^W) / 1000; % KN
        Fmax_SingleRoot(i_force) = max(Fcurve_SingleRoot);

    end

    % tensile force for bundle of roots
    Ftot = Fmax_SingleRoot .* Ni;  % kN
    root_cohesion = Ftot ./ RootSoil_area;   % kPa
    
    % calculate the root cohesion at specific depth
%     if Landslide_depth <= Zr
%         DepthLocation_matrix = abs(SoilDepth - Landslide_depth);
%         LandDepth_index = find(DepthLocation_matrix == min(DepthLocation_matrix));
%         if length(LandDepth_index) > 1
%             LandDepth_index = LandDepth_index(1);
%         end
%         Veg_cohesion = root_cohesion(LandDepth_index);
%     else
%         Veg_cohesion = 0;
%     end


end












