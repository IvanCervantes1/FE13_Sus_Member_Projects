% Function to find the loads through each suspension arm. Takes input of
% car parameters and g-forces entered in a matrix of corresponding 
% longitudinal (x) and lateral (y) pairs. The gg data must be organized in
% this order to produce a valid return

% This model assumes that the lateral force of each tire is propoteional
% to its share of the whole vehicle's lateral force (per the suspension
% solver sheet). It will be further expanded in the future for more complex
% cases

% This function plots load vs angle around gg diagram. The next iteration
% will plot load vs time

% Run this program with the fe12params structure and the GGV_dataComb array

function forces = suspensionSolver(carParams, Gs)
    inboard = carParams.inboard;
    outboard = carParams.outboard;
    locations = carParams.fAppLoc;
    momentArm = zeros(1,3);                                                % Measurements of moments will be taken relative to the origin, this can be changed
    A = zeros(6);                                                          % Initializes a matrix to store unit vectors and moment vectors
    forces = zeros(size(Gs,1),9);                                          % Initializes a matrix to store the forces with a row for every pair of lat/long acceleration pairs entered
    forces(:,8:9) = Gs;
    loadTable = loadCases(carParams, Gs);                                  % Creates a table of x, y, and z loads using car parameters and lat/long acceleration pairs
    locations = locations-momentArm;                                       % Normalizes the force application locations to a chosen moment arm
    for i = 1:6                                                            % For each arm
        link = (outboard(i,:)-inboard(i,:))';                              % Finds the distance between the inboard and outboard mounting points in terms of x, y, and z
        A(1:3,i) = link/norm(link);                                        % Turns each x, y, and z value into a unit vector and places said vector in its respective indexes within A
    end
    inboard = inboard - momentArm;                                         % Normalizes the inboard location to a chosen moment arm (zero in this case, but included in case a location other than the origin is chosen)
    for i = 1:6                                                            % For each arm
        A(4:6,i) = cross(inboard(i,:), A(1:3,i)');                         % Calculates the moment vector and places said vector in its respective indexes within A
    end
    for i = 1:size(Gs,1)                                                   % For each lat/long acceleration pair
        theta = cart2pol(Gs(i,2),Gs(i,1))*180/pi;                          % Finds the radial coordinate theta for each pair
        if theta < 0                                                       % Ensures all theta values are between 0 and 360 instead of -180 and 180
            theta = 360+theta;
        end
        forces(i, 1) = theta;
        fApplied = eye(3).*loadTable(i, 6:8)';                             % Creates a diagonal matrix with the appropriate x, y, and z forces for the respective lat/long acceleration pair
        mApplied = [cross(locations(1,:),fApplied(1,:)); cross(locations(2,:),fApplied(2,:)); cross(locations(3,:),fApplied(3,:))]; % Creates a matrix of applied moments by crossing the applied force locations with applied forces
        x = [-fApplied(1,1); -fApplied(2,2); -fApplied(3,3); -sum(mApplied)'];      % Fills column vector x with forces applied and moments applied
        forces(i,2:7) = (A\x)';                                            % Fills the respective row of the forces matrix with the forces through each arm
    end
    forces = sortrows(forces);                                             % Sorts forces by angle around gg diagram
    plot(forces(:, 1),forces(:,2),forces(:,1),forces(:,3),forces(:,1),forces(:,4),forces(:,1),forces(:,5),forces(:,1),forces(:,6),forces(:,1),forces(:,7)); % Load vs. gg diagram angle for all 6 arms
    xlim([0,360]);
    xlabel('gg Diagram Angle (Degrees, by polar coordinate system)')
    ylabel('Load (N)')
    title('Load vs. gg Diagram Angle')
    legend('Up-Fore','Up-Aft','Low-Fore','Low-Aft','Pushrod','Tie/Toe')
    clear A;
    clear x;
    clear fApplied;
    clear mApplied;
    clear link;
end

% Function to find forces and weight transfer from a structure of car
% parameters and an array of lateral and longitudinal G forces

function loadTable = loadCases(carParams, Gs)
    loadTable = zeros(size(Gs,1), 8);                                      % Initializes a load table matrix to store the outputs of this function
    for i = 1:size(Gs,1)                                                          % For each provided value of lat and long G's
        loadTable(i,1) = (carParams.m*Gs(i,2)*9.81*carParams.hCG)/carParams.TWf;  % Calculates lateral WT
        loadTable(i,2) = (carParams.m*Gs(i,1)*9.81*carParams.hCG)/carParams.WB;   % Calculates longitudinal WT
        loadTable(i,8) = carParams.m*9.81/4 + (loadTable(i,1)+loadTable(i,2))/2;  % Calculates Fz
        loadTable(i,3) = loadTable(i,8)/(carParams.m*9.81);                       % Calculates Fz%
        loadTable(i,4) = carParams.m*Gs(i,1)*9.81;                                % Calculates Fx_car
        loadTable(i,5) = carParams.m*Gs(i,2)*9.81;                                % Calculates Fy_car
        loadTable(i,6) = loadTable(i,4)*loadTable(i,3);                           % Calculates Fx
        loadTable(i,7) = loadTable(i,5)*loadTable(i,3);                           % Calculates Fy
    end
end