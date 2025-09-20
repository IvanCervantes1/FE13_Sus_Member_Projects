% Function to find the loads through each suspension arm. Takes input of
% car parameters and g-force (NOT ACCELERATION) data entered in a matrix of 
% corresponding longitudinal (x) and lateral (y) values. The g-force data 
% must be organized in this order to produce a valid return. 

% This function can also take special cases of load table data. This data 
% must be entered in the form [WT lateral, WT longitudinal, Fz% 
% (percent of vehicle weight on wheel of interest), Fx_car (total 
% longitidunal force on the car), Fy_car (total lateral force on the car), 
% Fx (long force on wheel of interest), Fy (like Fx but lat), Fz (like Fx 
% but vertical)]. This is used for special cases like max bump or reverse 
% braking, where the load table algorithms don't apply.

% This model assumes that the lateral force of each tire is proportional
% to its share of the whole vehicle's lateral force (per the suspension
% solver sheet)

% To find ONLY the forces for each load case, run this program with the 
% fe12params structure and an array of Gs in the format:
% [forcesF, ForcesR] = suspensionSolverPlot(fe12params, Gs)

% To find the forces for each load case AND the overall max loads, 
% run this program with the fe12params structure and an array of Gs in the format:
% [forcesF, ForcesR, frontMaxes, rearMaxes] = suspensionSolverPlot(fe12params, Gs)


function [forcesF, forcesR, frontMaxes, rearMaxes] = suspensionSolverCases(carParams, Gs, specialCases)
    inboardF = carParams.inboardF;
    outboardF = carParams.outboardF;
    locationsF = carParams.tireContactPtF;
    inboardR = carParams.inboardR;
    outboardR = carParams.outboardR;
    locationsR = carParams.tireContactPtR;
    momentArm = zeros(1,3);                                                % Measurements of moments will be taken relative to the origin, this can be changed
    A = zeros(6);                                                          % Initializes a matrix to store unit vectors and moment vectors
    [loadTableF, loadTableR] = loadCases(carParams, Gs(:, 1:2)*9.81); % Creates a table of x, y, and z loads using car parameters and acceleration pairs
    exist('specialCases', 'var')
    if exist('specialCases', 'var')
        for i = 1:length(specialCases(1))
            loadTableF(size(loadTableF,1)+i, :) = specialCases;
            loadTableR(size(loadTableR,1)+i, :) = specialCases;
        end
    end
    forcesF = zeros(size(loadTableF,1),6);                                 % Initializes a matrix to store the forces with a row for each recorded acceleration pair
    forcesR = zeros(size(loadTableR,1),6);
    locationsF = locationsF-momentArm;                                     % Normalizes the force application locations to the chosen moment arm
    locationsR = locationsR-momentArm;
    for i = 1:6                                                            % For each arm
        link = (outboardF(i,:)-inboardF(i,:))';                            % Finds the distance between the inboard and outboard mounting points in terms of x, y, and z
        A(1:3,i) = link/norm(link);                                        % Turns each x, y, and z value into a unit vector and places said vector in its respective indexes within A
    end
    inboardF = inboardF - momentArm;                                       % Normalizes the inboard location to the chosen moment arm
    for i = 1:6                                                            % For each arm
        A(4:6,i) = cross(inboardF(i,:), A(1:3,i)');                        % Calculates the moment vector and places said vector in its respective indexes within A
    end
    for i = 1:size(loadTableF,1)                                            % For each lat/long acceleration pair
        fApplied = eye(3).*loadTableF(i, 6:8)';                             % Creates a diagonal matrix with the appropriate x, y, and z forces for the respective lat/long acceleration pair
        mApplied = [cross(locationsF(1,:),fApplied(1,:)); cross(locationsF(2,:),fApplied(2,:)); cross(locationsF(3,:),fApplied(3,:))]; % Creates a matrix of applied moments by crossing the applied force locations with applied forces
        x = [-fApplied(1,1); -fApplied(2,2); -fApplied(3,3); -sum(mApplied)'];      % Fills column vector x with forces applied and moments applied
        forcesF(i,:) = (A\x)';                                             % Fills the respective row of the forces matrix with the forces through each arm
    end
    frontMaxes = findMaxes(forcesF);
    frontMaxes.Properties.VariableNames = {'Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Tie Rod'};
    frontMaxes.Properties.RowNames = {'Max Compression', 'Max Tension'};
    forcesF = array2table(forcesF);
    forcesF.Properties.VariableNames = {'Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Tie Rod'};
    for i = 1:6                                                            
        link = (outboardR(i,:)-inboardR(i,:))';                            
        A(1:3,i) = link/norm(link);                                        
    end
    inboardR = inboardR - momentArm;                                  
    for i = 1:6                                                            
        A(4:6,i) = cross(inboardR(i,:), A(1:3,i)');                        
    end
    for i = 1:size(loadTableR,1)                                                 
        fApplied = eye(3).*loadTableR(i, 6:8)';                            
        mApplied = [cross(locationsR(1,:),fApplied(1,:)); cross(locationsR(2,:),fApplied(2,:)); cross(locationsR(3,:),fApplied(3,:))]; % Creates a matrix of applied moments by crossing the applied force locations with applied forces
        x = [-fApplied(1,1); -fApplied(2,2); -fApplied(3,3); -sum(mApplied)'];      
        forcesR(i,:) = (A\x)';     
    end
    rearMaxes = findMaxes(forcesR);
    rearMaxes.Properties.VariableNames = {'Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Toe Rod'};
    rearMaxes.Properties.RowNames = {'Max Compression', 'Max Tension'};
    forcesR = array2table(forcesR);
    forcesR.Properties.VariableNames = {'Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Toe Rod'};
end

% Function to find forces and weight transfer from a structure of car
% parameters and an array of lateral and longitudinal acceleration

function [loadTableF, loadTableR] = loadCases(carParams, accelData)
    loadTableF = zeros(size(accelData,1), 8);                              % Initializes a load table matrix to store the outputs of this function
    loadTableR = loadTableF;
    for i = 1:size(accelData,1)                                            % For each provided value of lat and long G's
        [loadTableF(i,1), loadTableR(i,1)] = SampoWeightTransfer(carParams, accelData(i,2));
        loadTableF(i,2) = (carParams.m*accelData(i,1)*carParams.hCG)/carParams.WB;   % Calculates longitudinal WT
        loadTableF(i,8) = carParams.m*9.81/4+loadTableF(i,1)+loadTableF(i,2)/2;    % Calculates Fz
        loadTableF(i,3) = loadTableF(i,8)/(carParams.m*9.81);              % Calculates Fz%
        loadTableF(i,4) = carParams.m*accelData(i,1);                 % Calculates Fx_car
        loadTableF(i,5) = carParams.m*accelData(i,2);                 % Calculates Fy_car
        loadTableF(i,6) = loadTableF(i,4)*loadTableF(i,3);                 % Calculates Fx
        loadTableF(i,7) = loadTableF(i,5)*loadTableF(i,3);                 % Calculates Fy
        loadTableR(i,2) = (carParams.m*accelData(i,1)*carParams.hCG)/carParams.WB;   % Calculates longitudinal WT
        loadTableR(i,8) = carParams.m*9.81/4+loadTableR(i,1)+loadTableR(i,2)/2;    % Calculates Fz
        loadTableR(i,3) = loadTableR(i,8)/(carParams.m*9.81);                       % Calculates Fz%
        loadTableR(i,4) = carParams.m*accelData(i,1);                               % Calculates Fx_car
        loadTableR(i,5) = carParams.m*accelData(i,2);                               % Calculates Fy_car
        loadTableR(i,6) = loadTableR(i,4)*loadTableR(i,3);                          % Calculates Fx
        loadTableR(i,7) = loadTableR(i,5)*loadTableR(i,3);                          % Calculates Fy
    end
end

% Function to calculate from and rear weight transfer using chassis
% flexibility

function [deltaFzFront, deltaFzRear] = SampoWeightTransfer(carParams, Ay)
    kF = carParams.kF;                                                     % For the sake of simplicity, members of carParams structure are assigned
    kR = carParams.kR;                                                     % to shorter local variables within this function
    kC = carParams.kC;
    m_uF = carParams.m_uF;
    m_uR = carParams.m_uR;
    TWf = carParams.TWf;
    TWr = carParams.TWr;
    h_uF = carParams.h_uF;
    h_uR = carParams.h_uR;
    zF = carParams.zF;
    zR = carParams.zR;
    d_sF = carParams.h_sF - zF;
    d_sR = carParams.h_sR - zR;
    m_sF = carParams.m_s*carParams.b_s/carParams.WB;
    m_sR = carParams.m_s*carParams.a_s/carParams.WB;
    deltaFzFront = ((kF*d_sF*m_sF)/(kF+(kR*kC/(kR+kC)))+((kF*kC/(kF+kC))*d_sR*m_sR)/(kF*kC/(kF+kC)+kR)+zF*m_sF+h_uF*m_uF)*Ay/TWf;
    deltaFzRear = ((kR*kC/(kR+kC)*d_sF*m_sF)/(kF+kR*kC/(kR+kC))+(kR*d_sR*m_sR)/(kF*kC/(kF+kC)+kR)+zR*m_sR+h_uR*m_uR)*Ay/TWr;
end

function maxes = findMaxes(values)
    maxVals = zeros(2,6);
    for i = 1:size(values, 1)
        for j = 1:6
            if values(i,j) > maxVals(1,j)
                maxVals(1,j) = values(i,j);
            elseif values(i,j) < maxVals(2,j)
                maxVals(2,j) = values(i,j);
            end
        end
    end
    maxes = array2table(maxVals);
end