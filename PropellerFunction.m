function OutputCoordinates = PropellerFunction(PropData)

%% Derived variables

PropData.numPoints = (PropData.NoDataPointsBlade + PropData.NoDataPointsGap)*PropData.numBlades;     % Prop.numPoints = number of points that define circumference of the blades
PropData.r = PropData.BladeDiameter + PropData.HubRadius;                                                  % PropData.r = radius of blades + hub
PropData.c = PropData.r*2*pi;                                                                              % PropData.c = outer circumference of propeller                                   
PropData.MaxArcAngle = (2*pi/PropData.numBlades);                                                        % PropData.MaxArcAngle = the maximum angle allowed for each blade and gap    

%% Output Matrices

OutputCoordinates.X = [];       % Create empty output matrix for X
OutputCoordinates.Y = [];       % Create empty output matrix for Y
OutputCoordinates.Z = [];       % Create empty output matrix for Z

PropData.endPos=0;              % Z coordinate for bottom layer of cylinder

%% Create Cross-Section Geometry

r_CrossSections = linspace(PropData.HubRadius, PropData.r, PropData.NoCrossSections);             % Linearly space the cross-sections from the hub to the tip of the blade

%Output [Max Camber , Max Camber Pos, Max Thickness] NACA Matrix, ChordLength Matrix and Pitch Angle Matrix
[NACA, ChordLength, PitchAngle] = GeometryVariation(PropData);


%% Create Aerofoil Crossections

for i = 1:PropData.NoCrossSections %Starting at hub we create all aerofoil cross-sections data points up to the blade tip
    
    [XX, YY, ZZ, AngleStore] = CreateAerofoilDataPoints(r_CrossSections(i), PitchAngle(i), ChordLength(i), PropData, NACA(i,:)); 
    
    %Separate Upper and Lower edges of Aerofoil Cross-section to rearrange
    %order in final matrix
    XX_L(i,:) = XX(1,:);
    XX_U(i,:) = XX(2,:);
    YY_L(i,:) = YY(1,:);
    YY_U(i,:) = YY(2,:);
    ZZ_L(i,:) = ZZ(1,:);
    ZZ_U(i,:) = ZZ(2,:);
    
    if i == 1 %Store angles for which x,y ,z were evaluated at the hub cross-section to create cylinder 
        CylinderAngles = AngleStore;
    end
end

%Reverse order of Upper edge points for correct construction of blade (the
%lower edge will be created first, and the upper edge then in reverse
%order)
XX_U = flipud(XX_U);
YY_U = flipud(YY_U);
ZZ_U = flipud(ZZ_U);

%Concatenate Lower Edge Values above Top Edge Values
OutputCoordinates = ConcatenateResults(OutputCoordinates, XX_L, YY_L, ZZ_L);
OutputCoordinates = ConcatenateResults(OutputCoordinates, XX_U, YY_U, ZZ_U);


%% Create Outer Cylinder
HubMax_L = min(ZZ_L(1,:));                      % HubMax_L is the lowest z value on the hub where the blade attaches
HubMax_U = max(ZZ_U(end,:));                    % HubMax_U is the highest z value one the hub where the blade attaches
HubMax_Height = HubMax_U-HubMax_L;              % HubMax_Height is the height of the hub from the minimum and maximum attachment coordinates points

[xx,yy,zz] = CreateCylinderCoordinates(PropData.HubRadius, CylinderAngles);
zz(1,:) =(PropData.endPos + HubMax_L - PropData.HubHeight); %The z value at the bottom of the cyinder = (the lowest z value for the hub blade cross-section - the additional hub height) 
zz(2,:) =(PropData.endPos + HubMax_U + PropData.HubHeight); %The z value at the top of the cyinder = (the highest z value for the hub blade cross-section + the additional hub height) 
OutputCoordinates = ConcatenateCylinderResults(OutputCoordinates, xx, yy, zz, 1);

%% Create Inner Cylinder
[xx,yy,zz] = CreateCylinderCoordinates(PropData.HubRadius*0.25, CylinderAngles);
zz(1,:) =(PropData.endPos + HubMax_L - PropData.HubHeight);
OutputCoordinates = ConcatenateCylinderResults(OutputCoordinates, xx, yy, zz, 0 );

%% Create Nose
No_Discretizations = 100; %Number of z points to evaluate at to create nose
 
% Create nose arc using equation of circle
z_curv = linspace(0,PropData.HubRadius,No_Discretizations); 
NoseR = sqrt(PropData.HubRadius.^2 - z_curv.^2); %Calculate radius for each z value
x_Nose = NoseR' * cos(CylinderAngles);
y_Nose = NoseR' * sin(CylinderAngles);

CylinderHeight = HubMax_U + PropData.HubHeight; %z starting height of nose

z_Nose = linspace(CylinderHeight,CylinderHeight+PropData.NoseHeight,No_Discretizations)';
z_Nose = repmat(z_Nose,1,PropData.numPoints); %replicate z values for number of r values

OutputCoordinates = ConcatenateResults(OutputCoordinates,x_Nose, y_Nose, z_Nose); %Concatenate nose points to the end of matrix


%% Plot Figure
mesh(OutputCoordinates.X,OutputCoordinates.Y,OutputCoordinates.Z)
axis equal

end 


%% END

%% Function to create thickness variation along the cross-sections

function [NACA, ChordLength, PitchAngle] = GeometryVariation(PropData)
% This function takes the input matrices given for the MaxCamber,
% Max Camber Position, Max Thickness, Chord Lengths and Pitch angles and
% interpolates using a spline in order to calculate the values of each of these
% variables at the radii the blade aerofoil cross-sections are being created

%Build NACA Matrix
No_MaxCamber = size(PropData.MaxCamber,2)-1;                                                                            %The number of inputs in the MaxCamber matrix minus one, since the x values given to the interp1 function start at 0 and increment in 1. The x values can be set this way as we assume the Camber Values are linearly spaced across the blade diameter.
MaxCamber = interp1(0:1:No_MaxCamber,PropData.MaxCamber,linspace(0,No_MaxCamber,PropData.NoCrossSections),'spline');  % Spline interpolation.  Linspace is used to create the evaluation points spaced equally to the number of cross-sections. 

No_PosMaxCamber = size(PropData.PosMaxCamber,2)-1;
PosMaxCamber = interp1(0:1:No_PosMaxCamber,PropData.PosMaxCamber,linspace(0,No_PosMaxCamber,PropData.NoCrossSections),'spline');

No_MaxThickness = size(PropData.MaxThickness,2)-1;
MaxThickness = interp1([0:1:No_MaxThickness,No_MaxThickness+0.05],[PropData.MaxThickness, 0],linspace(0,No_MaxThickness,PropData.NoCrossSections),'spline');

NACA =  [MaxCamber', PosMaxCamber', MaxThickness']; %Create NACA matrix containing only aerofoil geometry inputs  

%Build Chord Length Matrix
No_ChordLength = size(PropData.ChordLength,2)-1;
ChordLength = interp1(0:1:No_ChordLength,PropData.ChordLength,linspace(0,No_ChordLength,PropData.NoCrossSections),'spline');

%Build Pitch Angle Matrix
No_PitchAngle = size(PropData.PitchAngle,2)-1;
PitchAngle = interp1(0:1:No_PitchAngle,PropData.PitchAngle,linspace(0,No_PitchAngle,PropData.NoCrossSections),'spline');

end

%% Main Function to create data points

function [XX, YY, ZZ, AngleStore] = CreateAerofoilDataPoints(r, PitchAngle, ChordLength, PropData, NACA) 

%This function creates the aerofoil data points given a set radius from
%origin, pitch angle, chord length and NACA matrix 

arcLength = PropData.MaxArcAngle;                   % arclength is the allocated angle to each blade
L1 = ChordLength/r;                                 % L1 is the angle that will create a cross-section with required ChordLength

ic = 0;                                             % starting value for the counter
BladePosition = cell(1, PropData.numBlades);        % cell array is used to store only data points associated with blades, each cell entry entry will be matrix of the indices related to each blade
AngleStore = [];                                    % Matrix that will store all the angles used to create blade/gap arcs, used later to create the same number of hub and nose circle data points

for i = 1:PropData.numBlades %For each blade, we define the arc points in the (X,Y) directions where the cross-section will be created. And set the gap points at the hub radius. 
    
    for j = linspace(arcLength*(i-1),(L1+arcLength*(i-1)), PropData.NoDataPointsBlade) %j is the angle of evaluation for the blade cross-section. linspace ensures the number of data points used is exactly PropData.NoDataPointsBlade
        ic=ic+1;
        x(ic)=r*cos(j);
        y(ic)=r*sin(j);
        AngleStore(ic) = j;
        BladePosition{i}(end+1) = ic;
    end
    
    for k = linspace((L1+arcLength*(i-1)),arcLength*i, PropData.NoDataPointsGap) %k is the angle of evaluation for the gap points. linspace ensures the number of data points used is exactly PropData.NoDataPointsGap
        ic=ic+1;
        x(ic)=PropData.HubRadius*cos(k);
        y(ic)=PropData.HubRadius*sin(k);
        AngleStore(ic) = k;
    end
    
end

XX = repmat(x,2,1);             % Replicate first row of x coordinates to give aerofoil upper surface values
YY = repmat(y,2,1);             % Replicate first row of y coordinates to give aerofoil upper surface values

ZZ = zeros(size(XX,1),length(YY));      % Form ZZ coordinate matrix with the same dimensions as XX and YY
ZZ(1,:)= PropData.endPos;               % Give first and second row intial z value of PropData.endPos which is 0, the blade z points will be overwritten in the following code, leaving the gap points at 0 height.
ZZ(2,:)=PropData.endPos;


% Create the aerofoil data points (yU and yL will be used as the ZZ values
% for the final blade data points), xU and xL are still needed for
% rotation. 
[xU, yU, xL, yL] = CreateAerofoil(NACA(1),NACA(2),NACA(3),PropData.NoDataPointsBlade);

%x Values are normalized from 0 - 1, and will be scaled up to the
%ArcLength. So y values must be scaled to the same length to keep
%profile proportions correct 
yU = ChordLength.*yU;
yL = ChordLength.*yL;

% Rotate aerofoil by pitch angle 
[xU, yU] = Rotate(PitchAngle, xU, yU);
[xL, yL] = Rotate(PitchAngle, xL, yL);

%Assign the ZZ values using the BladePosition store earlier saved with the
%blade indices 
for i = 1:PropData.numBlades %Loop for each blade, in the BladePosition store       
    % Assign Z values
    ZZ(1,BladePosition{i})= yL';
    ZZ(2,BladePosition{i})= yU';
end

end 

%% Function to attain x and y coordinates of the aerofoil cross-sections

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

[yc, dy_dx, theta] = CamberFunction(x, PosMaxCamber_Norm, MaxCamber_Norm);

% Thickness distribution

yt = 5*MaxThickness_Norm.*((t0.*sqrt(x)) + (t1.*x) + (t2.*x.^2) + (t3.*x.^3) + (t4.*x.^4));

% Upper surface points 
xU = x(:) - yt(:).*sin(theta);
yU = yc(:) + yt(:).*cos(theta);         

% Lower surface points 
xL = x(:) + yt(:).*sin(theta);
yL = yc(:) - yt(:).*cos(theta);


end

function [yc, dy_dx, theta] = CamberFunction(x, PosMaxCamber_Norm, MaxCamber_Norm)
    
    %Creates the camber data points and gradients. Up to the Maximum Camber
    %Position, one set of equations are used, and from the Maximum Camber
    %position another set of equations are used. 

    condition1 = (x >= 0).*(x < PosMaxCamber_Norm); %Create a matrix of 1s and 0s for all matrix where 1 is when the condition is true and 0 is when not
    yc1       = (MaxCamber_Norm./PosMaxCamber_Norm.^2).*((2.*PosMaxCamber_Norm.*x)-x.^2).*condition1; %Multiply the condition matrix by the calculated equation values to leave values only in places where it will be used
    dy_dx1   = ((2.*MaxCamber_Norm)./(PosMaxCamber_Norm.^2)).*(PosMaxCamber_Norm-x).*condition1;
    
    condition2 = (x >= PosMaxCamber_Norm) .* (x <= 1);
    yc2       = (MaxCamber_Norm./(1-PosMaxCamber_Norm).^2).*(1-(2.*PosMaxCamber_Norm)+(2.*PosMaxCamber_Norm.*x)-(x.^2)).*condition2;
    dy_dx2   = ((2.*MaxCamber_Norm)./((1-PosMaxCamber_Norm).^2)).*(PosMaxCamber_Norm-x).*condition2;
    
    yc = yc1+yc2; %Add both matrices together to get the final output matrices
    dy_dx = dy_dx1+dy_dx2;
    
    theta = atan(dy_dx); %arctan the gradients to get the theta values
end



%% Function to apply pitch angle to cross-sections

function [X_Rot, Y_Rot] = Rotate(Angle, X, Y)  

%Rotate all coordinates using equations from https://www.siggraph.org/education/materials/HyperGraph/modeling/mod_tran/2drota.htm

X_Rot = X.*(cos(degtorad(Angle))) - Y.*(sin(degtorad(Angle)));              % Equation to rotate x coordinate values
Y_Rot = Y.*(cos(degtorad(Angle))) + X.*(sin(degtorad(Angle)));              % Equation to rotate y coordinate values

end 


%% Function to vertically concatenate matrices

function OutputCoordinates = ConcatenateResults(OutputCoordinates,xx, yy, zz)
%This function concatenates xx,yy,zz values to the end of the final Output Coordinates
%matrices

OutputCoordinates.X = vertcat(OutputCoordinates.X,xx);
OutputCoordinates.Y = vertcat(OutputCoordinates.Y,yy);
OutputCoordinates.Z = vertcat(OutputCoordinates.Z,zz);
end

%% Function to create cylinder data points

function [xx,yy,zz] = CreateCylinderCoordinates(r, theta)
%Generates the cylinder coordinates for the radius r, and theta values
%calculated 

n = size(theta,2);
x = r .* cos(theta);
y = r .* sin(theta);

xx = vertcat(x,x);
yy = vertcat(y,y);

zz(1,:) = zeros(1,n);
zz(2,:) = ones(1,n);
end

%% Function to concatanete cylinder data points into output matrices

function OutputCoordinates = ConcatenateCylinderResults(OutputCoordinates,xx, yy, zz, OuterOrInner)
if  OuterOrInner %=1 if defining outer cylinder
    
    %The outer cylinder's bottom edge values to before the Aerofoil data
    %points. The outer cylinder's top edge values will be added to the end of
    %the matrix.
    
    OutputCoordinates.X = vertcat(OutputCoordinates.X,xx(2,:));
    OutputCoordinates.Y = vertcat(OutputCoordinates.Y,yy(2,:));
    OutputCoordinates.Z = vertcat(OutputCoordinates.Z,zz(2,:));
    OutputCoordinates.X = vertcat(xx(1,:),OutputCoordinates.X);
    OutputCoordinates.Y = vertcat(yy(1,:), OutputCoordinates.Y);
    OutputCoordinates.Z = vertcat(zz(1,:), OutputCoordinates.Z);
else
    
    %For the inner cylinder, only the bottom edge needs to be created and appended to the first row as this is where the mesh starts. 
    
    OutputCoordinates.X = vertcat(xx(1,:),OutputCoordinates.X);
    OutputCoordinates.Y = vertcat(yy(1,:), OutputCoordinates.Y);
    OutputCoordinates.Z = vertcat(zz(1,:), OutputCoordinates.Z);
end 
end
