% Function to find the loads through each suspension arm. Takes input of
% car parameters and acceleration data entered in a matrix of corresponding 
% time, longitudinal (x), and lateral (y) values. The acceleration data 
% must be organized in this order to produce a valid return

% This model assumes that the lateral force of each tire is proportional
% to its share of the whole vehicle's lateral force (per the suspension
% solver sheet)

% This function plots load vs time for both front and rear suspension and
% returns maximum tension and compression values in table format

% Run this program with the fe12params structure and an array of
% accelData in the format:
% [frontMaxes, rearMaxes] = suspensionSolverPlot(fe12params, accelData)

function [frontMaxes, rearMaxes] = suspensionSolverPlot(carParams, accelData)
    accelData(:,1) = accelData(:,1)-accelData(1,1);
    inboardF = carParams.inboardF;
    outboardF = carParams.outboardF;
    locationsF = [carParams.tireContactPtF(1,1), carParams.tireContactPtF(1,2), -8*25.4];
    locationsF = [locationsF; locationsF; locationsF; locationsF; locationsF; locationsF];
    inboardR = carParams.inboardR;
    outboardR = carParams.outboardR;
    locationsR = [carParams.tireContactPtR(1,1), carParams.tireContactPtR(1,2), -8*25.4];
    locationsR = [locationsR; locationsR; locationsR; locationsR; locationsR; locationsR];
    RF = locationsF-[carParams.tireContactPtF; carParams.tireContactPtF];
    RR = locationsR-[carParams.tireContactPtR; carParams.tireContactPtR];
    rF = outboardF - locationsF;
    rR = outboardR - locationsR;
    A = zeros(6);                                                          % Initializes a matrix to store unit vectors and moment vectors
    [loadTableF, loadTableR] = loadCases(carParams, accelData(:, 2:3)/9.81);    % Creates a table of x, y, and z loads using car parameters and acceleration pairs
    forcesF = zeros(size(loadTableF,1),6);                                 % Initializes a matrix to store the forces with a row for each recorded acceleration pair
    forcesR = zeros(size(loadTableR,1),6);
    for i = 1:6                                                            % For each arm
        link = (outboardF(i,:)-inboardF(i,:))';                            % Finds the vector components of each link
        norm(link);
        A(1:3,i) = link/norm(link);                                        % Turns each component into a unit vector and places it within A
    end
    for i = 1:6                                                            % For each arm
        A(4,i) = A(3,i)*rF(i,2)-A(2,i)*rF(i,3);                            % Calculates the moment vector and places said vector in A
        A(5,i) = A(3,i)*rF(i,1)-A(1,i)*rF(i,3);
        A(6,i) = A(2,i)*rF(i,1)-A(1,i)*rF(i,2);
    end
    for i = 1:size(loadTableF,1)                                           % For each lat/long acceleration pair
        B(1:3,1) = loadTableF(i, 3:5)';                                    % Creates a diagonal matrix with x, y, and z forces
        B(4,1) = B(3,1)*RF(1,2)-B(2,1)*RF(1,3);
        B(5,1) = B(1,1)*RF(1,3)-B(3,1)*RF(1,1);
        B(6,1) = B(2,1)*RF(1,1)-B(1,1)*RF(1,2);
        forcesF(i,:) = (A\B)';                                             % Fills the respective row of the forces matrix with the forces through each arm
    end
    
    subplot(2,1,1);
    plot(accelData(:, 1),forcesF(:,1),accelData(:,1),forcesF(:,2),accelData(:,1),forcesF(:,3),accelData(:,1),forcesF(:,4),accelData(:,1),forcesF(:,5),accelData(:,1),forcesF(:,6), 'LineWidth', 1.5); % Load vs. gg diagram angle for all 6 arms
    xlim([0,accelData(size(accelData,1),1)])
    ylim([-5000, 5000])
    xlabel('Time (s)', 'FontSize', 20)
    ylabel('Load (N)', 'FontSize', 20)
    title('Load vs. Time, Front Left', 'FontSize', 20)
    legend('Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Tie Rod')
    
    frontMaxes = findMaxes(forcesF);
    frontMaxes.Properties.VariableNames = {'Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Tie Rod'};
    frontMaxes.Properties.RowNames = {'Max Tension', 'Max Compression'};

    A = zeros(6);

    for i = 1:6                                                            
        link = (outboardR(i,:)-inboardR(i,:))';                            
        A(1:3,i) = link/norm(link);                                        
    end                                  
    for i = 1:6                                                            
        A(4,i) = A(3,i)*rR(i,2)-A(2,i)*rR(i,3);                            % Calculates the moment vector and places said vector in A
        A(5,i) = A(3,i)*rR(i,1)-A(1,i)*rR(i,3);
        A(6,i) = A(2,i)*rR(i,1)-A(1,i)*rR(i,2);                      
    end
    for i = 1:size(loadTableR,1) 
        B(1:3,1) = loadTableR(i, 3:5)';                                    % Creates a diagonal matrix with x, y, and z forces
        B(4,1) = B(3,1)*RR(1,2)-B(2,1)*RR(1,3);
        B(5,1) = B(1,1)*RR(1,3)-B(3,1)*RR(1,1);
        B(6,1) = B(2,1)*RR(1,1)-B(1,1)*RR(1,2);
        forcesR(i,:) = (A\B)';     
    end

    subplot(2,1,2);
    plot(accelData(:, 1),forcesR(:,1),accelData(:,1),forcesR(:,2),accelData(:,1),forcesR(:,3),accelData(:,1),forcesR(:,4),accelData(:,1),forcesR(:,5),accelData(:,1),forcesR(:,6), 'LineWidth', 1.5); % Load vs. gg diagram angle for all 6 arms
    xlim([0,accelData(size(accelData,1),1)])
    ylim([-5000, 5000])
    xlabel('Time (s)', 'FontSize', 20)
    ylabel('Load (N)', 'FontSize', 20)
    title('Load vs. Time, Rear Left', 'FontSize', 20)
    legend('Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Toe Rod')

    rearMaxes = findMaxes(forcesR);
    rearMaxes.Properties.VariableNames = {'Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Toe Rod'};
    rearMaxes.Properties.RowNames = {'Max Tension', 'Max Compression'};
end

% Function to find forces and weight transfer from a structure of car
% parameters and an array of lateral and longitudinal G forces

function [loadTableF, loadTableR] = loadCases(carParams, accelData)
    loadTableF = zeros(size(accelData,1), 5);                              % Initializes a load table matrix to store the outputs of this function
    loadTableR = loadTableF;
    for i = 1:size(accelData,1)                                            % For each provided value of lat and long G's
        [loadTableF(i,1), loadTableR(i,1)] = SampoWeightTransfer(carParams, accelData(i,2)*9.81);
        loadTableF(i,2) = (carParams.m*9.81*accelData(i,1)*carParams.hCG)/carParams.WB;   % Calculates longitudinal WT
        loadTableF(i,5) = -(carParams.m*9.81*carParams.PFront/2+loadTableF(i,1)-loadTableF(i,2)/2);    % Calculates Fz
        loadTableR(i,2) = (carParams.m*9.81*accelData(i,1)*carParams.hCG)/carParams.WB;   % Calculates longitudinal WT
        loadTableR(i,5) = -(carParams.m*9.81*(1-carParams.PFront)/2+loadTableR(i,1)+loadTableR(i,2)/2);    % Calculates Fz
        if accelData(i,1) > 0 && accelData(i,2) == 0
            loadTableF(i,3) = 0;
            loadTableR(i,3) = ((1.25*carParams.m*9.81*carParams.a_s)/carParams.WB)/(1-carParams.hCG/carParams.WB*1.25)/2;
        elseif accelData(i,1) < 0
            loadTableF(i,3) = 1.25*loadTableF(i,5);
            loadTableR(i,3) = 1.25*loadTableR(i,5);
        elseif accelData(i,1) == 0 && accelData(i,2) ~= 0
            loadTableF(i,3) = 0;
            loadTableR(i,3) = 0;
        elseif accelData(i,1) > 0 && accelData(i,2) ~= 0
            loadTableF(i,3) = 0;
            loadTableR(i,3) = 1.25*-loadTableR(i,5);
        end
        loadTableF(i,4) = loadTableF(i,5)*accelData(i,2);                  % Calculates Fy
        loadTableR(i,4) = loadTableR(i,5)*accelData(i,2);                  % Calculates Fy
    end
end

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