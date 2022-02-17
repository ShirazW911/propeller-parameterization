
function [xU, yU, xL, yL] = AerofoilViewFunction(PropData)

%This code only generates the aerofoil cross-section at a the percentage along the
%blade 

Percentage = PropData.CrossSectionDisplayed/100;

PropData.numPoints = (PropData.NoDataPointsBlade + PropData.NoDataPointsGap)*PropData.numBlades;     % Prop.numPoints = number of points that define circumference of the blades
PropData.r = PropData.BladeDiameter + PropData.HubRadius;                                                  % PropData.r = radius of blades + hub
PropData.c = PropData.r*2*pi;                                                                              % PropData.c = outer circumference of propeller                                   
PropData.MaxArcAngle = (2*pi/PropData.numBlades);

%Output [Max Camber , Max Camber Pos, Max Thickness] NACA Matrix,
%ChordLength Matrix and Pitch Angle Matrix
[CurrentValues] = GeometryVariation(PropData, Percentage);

% Identify the x and y values that represent the outer values of my blade
[xU, yU, xL, yL] = CreateAerofoil(CurrentValues(1),CurrentValues(2),CurrentValues(3),PropData.NoDataPointsBlade);

[xU, yU] = Rotate(CurrentValues(5), xU, yU);
[xL, yL] = Rotate(CurrentValues(5), xL, yL);

%Move to positive axes 
x = [xU; xL];
xTranslate = 0 - min(x); 
xU = xU + xTranslate; 
xL = xL + xTranslate; 

end


%% Function to apply pitch angle to cross-sections

function [X_Rot, Y_Rot] = Rotate(Angle, X, Y)  

X_Rot = X.*(cos(degtorad(Angle))) - Y.*(sin(degtorad(Angle)));              % Equation to rotate x coordinate values
Y_Rot = Y.*(cos(degtorad(Angle))) + X.*(sin(degtorad(Angle)));              % Equation to rotate y coordinate values

end 


function [CurrentValues] = GeometryVariation(PropData, Percentage)

%Build NACA Matrix
No_MaxCamber = size(PropData.MaxCamber,2)-1;
CurrentValues(1) = interp1(0:1:No_MaxCamber,PropData.MaxCamber,No_MaxCamber*Percentage,'spline'); %Use the interpolation to just calculate the value required at the percentage along blade

No_PosMaxCamber = size(PropData.PosMaxCamber,2)-1;
CurrentValues(2) = interp1(0:1:No_PosMaxCamber,PropData.PosMaxCamber,No_PosMaxCamber*Percentage,'spline');

No_MaxThickness = size(PropData.MaxThickness,2)-1;
CurrentValues(3) = interp1([0:1:No_MaxThickness,No_MaxThickness+0.05],[PropData.MaxThickness, 0],No_MaxThickness*Percentage,'spline');

%Build Chord Length Matrix
No_ChordLength = size(PropData.ChordLength,2)-1;
CurrentValues(4) = interp1(0:1:No_ChordLength,PropData.ChordLength,No_ChordLength*Percentage,'spline');

%Build Pitch Angle Matrix
No_PitchAngle = size(PropData.PitchAngle,2)-1;
CurrentValues(5) = interp1(0:1:No_PitchAngle,PropData.PitchAngle,No_PitchAngle*Percentage,'spline');

end

 function [xU, yU, xL, yL] = CreateAerofoil(MaxCamber, PosMaxCamber, MaxThickness, num_points)
            
            % Constants used for thickness distribution
            t0 = 0.2969;
            t1 = -0.1260;
            t2 = -0.3516;
            t3 = 0.2843;
            t4 = -0.1036;
            
            % Normalised Values
            MaxCamber_Norm = MaxCamber/100;                          % M% of chordlength
            PosMaxCamber_Norm = PosMaxCamber/10;                   % 10P% of chordlength
            MaxThickness_Norm = MaxThickness/100;                  % XX% of chordlength
            
            % Aerofoil grid
            x = linspace(0,1,num_points)';          % divide 1 into equal lenghth section
            
            % Camber and Gradient
            yc      = zeros (num_points,1);
            dyc_dx  = zeros (num_points,1);
            theta   = zeros (num_points,1);
            
            condition1 = (x >= 0).*(x < PosMaxCamber_Norm);
            yc1       = (MaxCamber_Norm./PosMaxCamber_Norm.^2).*((2.*PosMaxCamber_Norm.*x)-x.^2).*condition1;
            dy_dx1   = ((2.*MaxCamber_Norm)./(PosMaxCamber_Norm.^2)).*(PosMaxCamber_Norm-x).*condition1;
            
            condition2 = (x >= PosMaxCamber_Norm) .* (x <= 1);
            yc2       = (MaxCamber_Norm./(1-PosMaxCamber_Norm).^2).*(1-(2.*PosMaxCamber_Norm)+(2.*PosMaxCamber_Norm.*x)-(x.^2)).*condition2;
            dy_dx2   = ((2.*MaxCamber_Norm)./((1-PosMaxCamber_Norm).^2)).*(PosMaxCamber_Norm-x).*condition2;
            
            yc = yc1+yc2;
            dy_dx = dy_dx1+dy_dx2;
            
            theta = atan(dy_dx);            

            % Thickness distribution
            
            yt = 5*MaxThickness_Norm.*((t0.*sqrt(x)) + (t1.*x) + (t2.*x.^2) + (t3.*x.^3) + (t4.*x.^4));
            
            % Upper surface points
            xU = x(:) - yt(:).*sin(theta);
            yU = yc(:) + yt(:).*cos(theta);
            
            % Lower surface points
            xL = x(:) + yt(:).*sin(theta);
            yL = yc(:) - yt(:).*cos(theta);
                        
        end