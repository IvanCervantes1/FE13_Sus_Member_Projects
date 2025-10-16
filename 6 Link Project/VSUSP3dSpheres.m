function VSUSP3dSpheres(carParams, maxDroop, maxJounce)
    inboard = carParams.inboardF([1:4,6],:);
    outboard = carParams.outboardF([1:4,6],:);
    midTire = carParams.midTireF;
    tireOuterUpper = carParams.tireOuterUpperF;
    tireOuterFront = carParams.tireOuterFrontF;
    
    obLtoU = outboard(1,:)-outboard(3,:);
    obTtoU = outboard(1,:)-outboard(5,:);
    obTtoL = outboard(3,:)-outboard(5,:);
    obTtoM = midTire-outboard(5,:);

    upright = zeros(1,3);
    upright(1) = norm(obLtoU);
    upright(2) = norm(obTtoU);
    upright(3) = norm(obTtoL);

    camber0 = -1.693;

    toTireMid = zeros(1,3);
    toTireMid(1) = norm(midTire-outboard(1,:));
    toTireMid(2) = norm(midTire-outboard(3,:));
    toTireMid(3) = norm(midTire-outboard(5,:));

    toTireUpper = zeros(1,3);
    toTireUpper(1) = norm(tireOuterUpper-outboard(1,:));
    toTireUpper(2) = norm(tireOuterUpper-outboard(5,:));
    toTireUpper(3) = norm(tireOuterUpper-midTire);

    toTireFront = zeros(1,3);
    toTireFront(1) = norm(tireOuterFront-outboard(1,:));
    toTireFront(2) = norm(tireOuterFront-outboard(3,:));
    toTireFront(3) = norm(tireOuterFront-midTire);

    linksVec = outboard-inboard;
    linksMag = zeros(1,5);
    
    for i = 1:5
        linksMag(i) = norm(linksVec(i,:));
    end

    KPI0 = acosd(obLtoU(2)/norm(obLtoU(2:3)));
    caster0 = acosd(obLtoU(1)/norm(obLtoU([1,3])))-90;

    droopRange = 0:-0.1:maxDroop;
    jounceRange = 0.1:0.1:maxJounce;

    deltaKPI = zeros(1, length(droopRange)+length(jounceRange));
    caster = zeros(1, length(droopRange)+length(jounceRange));
    toe = zeros(1, length(droopRange)+length(jounceRange));
    mechTrail = zeros(1, length(droopRange)+length(jounceRange));
    scrub = zeros(1, length(droopRange)+length(jounceRange));
    camber1 = zeros(1, length(droopRange)+length(jounceRange));
    
    coordsNew = [outboard(1,:); outboard(3,:); outboard(5,:); midTire; tireOuterUpper; tireOuterFront];
    iLinkL = norm(outboard(3,[2,3])-inboard(3,[2,3]));

    for i = 1:length(droopRange)
        coordsNew(2,3) = outboard(3,3)-droopRange(i)*25.4;                % Adds bump to lower control arm z coordinate
        coordsNew(2,2) = inboard(3,2)-sqrt(iLinkL^2-(coordsNew(2,3)-inboard(3,3))^2);   % Finds new lower control arm y coordinate
        [x, y, z] = solveSpheres(inboard(1,:), inboard(2,:), coordsNew(2,:), linksMag(1), linksMag(2), upright(1), coordsNew(1,:));  % Finds new upper control arm coordinates
        coordsNew(1,:) = [x, y, z];
        KPI = acosd((coordsNew(1,2)-coordsNew(2,2))/norm(obLtoU(2:3)));   % Finds new KPI value
        deltaKPI(length(droopRange)+1-i) = KPI-KPI0;                                            % Finds change in KPI
        caster(length(droopRange)+1-i) = acosd((coordsNew(1,1)-coordsNew(2,1))/norm(obLtoU([1,3])))-90;      % Finds caster value
        [x, y, z] = solveSpheres(coordsNew(1,:), coordsNew(2,:), inboard(5,:), upright(2), upright(3), linksMag(5), coordsNew(3,:));
        coordsNew(3,:) = [x, y, z];
        [x, y, z] = solveSpheres(coordsNew(1,:), coordsNew(2,:), coordsNew(3,:), toTireMid(1), toTireMid(2), toTireMid(3), coordsNew(4,:));
        coordsNew(4,:) = [x, y, z];
        v = coordsNew(2,:)-coordsNew(1,:);
        t = (0-coordsNew(1,3))/v(3);
        x = coordsNew(1,1)+t*v(1);
        y = coordsNew(1,2)+t*v(2);
        mechTrail(length(droopRange)+1-i) = (x-coordsNew(4,1))/25.4;
        [x, y, z] = solveSpheres(coordsNew(1,:), coordsNew(3,:), coordsNew(4,:), toTireUpper(1), toTireUpper(2), toTireUpper(3), coordsNew(5,:));
        coordsNew(5,:) = [x, y, z];
        camber1(length(droopRange)+1-i) = acosd((coordsNew(4,3)-coordsNew(5,3))/norm((coordsNew(4,:)-coordsNew(5,:))));
        [x, y, z] = solveSpheres(coordsNew(1,:), coordsNew(2,:), coordsNew(4,:), toTireFront(1), toTireFront(2), toTireFront(3), coordsNew(6,:));
        coordsNew(6,:) = [x, y, z];
        toe(length(droopRange)+1-i) = acosd((coordsNew(6,1)-coordsNew(4,1))/norm((coordsNew(6,1:2)-coordsNew(4,1:2))));
    end
    
    coordsNew = [outboard(1,:); outboard(3,:); outboard(5,:); midTire; tireOuterUpper; tireOuterFront];

    for i = 1:length(jounceRange)
        coordsNew(2,3) = outboard(3,3)-jounceRange(i)*25.4;                % Adds bump to lower control arm z coordinate
        coordsNew(2,2) = inboard(3,2)-sqrt(iLinkL^2-(coordsNew(2,3)-inboard(3,3))^2);   % Finds new lower control arm y coordinate
        [x, y, z] = solveSpheres(inboard(1,:), inboard(2,:), coordsNew(2,:), linksMag(1), linksMag(2), upright(1), coordsNew(1,:));  % Finds new upper control arm coordinates
        coordsNew(1,:) = [x, y, z];
        KPI = acosd((coordsNew(1,2)-coordsNew(2,2))/norm(obLtoU(2:3)));   % Finds new KPI value
        deltaKPI(i+length(droopRange)) = KPI-KPI0;                                          % Finds change in KPI
        caster(i+length(droopRange)) = acosd((coordsNew(1,1)-coordsNew(2,1))/norm(obLtoU([1,3])))-90;      % Finds caster value
        [x, y, z] = solveSpheres(coordsNew(1,:), coordsNew(2,:), inboard(5,:), upright(2), upright(3), linksMag(5), coordsNew(3,:));
        coordsNew(3,:) = [x, y, z];
        [x, y, z] = solveSpheres(coordsNew(1,:), coordsNew(2,:), coordsNew(3,:), toTireMid(1), toTireMid(2), toTireMid(3), coordsNew(4,:));
        coordsNew(4,:) = [x, y, z];
        v = coordsNew(2,:)-coordsNew(1,:);
        t = (0-coordsNew(1,3))/v(3);
        x = coordsNew(1,1)+t*v(1);
        y = coordsNew(1,2)+t*v(2);
        mechTrail(i+length(droopRange)) = (x-coordsNew(4,1))/25.4;
        [x, y, z] = solveSpheres(coordsNew(1,:), coordsNew(3,:), coordsNew(4,:), toTireUpper(1), toTireUpper(2), toTireUpper(3), coordsNew(5,:));
        coordsNew(5,:) = [x, y, z];
        camber1(i+length(droopRange)) = acosd((coordsNew(4,3)-coordsNew(5,3))/norm((coordsNew(4,:)-coordsNew(5,:))));
        [x, y, z] = solveSpheres(coordsNew(1,:), coordsNew(2,:), coordsNew(4,:), toTireFront(1), toTireFront(2), toTireFront(3), coordsNew(6,:));
        coordsNew(6,:) = [x, y, z];
        toe(i+length(droopRange)) = acosd((coordsNew(6,1)-coordsNew(4,1))/norm((coordsNew(6,1:2)-coordsNew(4,1:2))));
    end

    camber = deltaKPI+camber0
    camber1
    caster
    mechTrail
    toe
    subplot(2,1,1);
    plot([flip(droopRange), jounceRange], camber, 'LineWidth', 2)
    xlabel('Front Left Bump (in)', 'FontSize', 20)
    ylabel('Front Left Camber (degrees)', 'FontSize', 20)
    title('Camber vs. Bump, Front Left', 'FontSize', 20)
    subplot(2,1,2);
    plot([flip(droopRange), jounceRange], caster, 'LineWidth', 2)
    xlabel('Front Left Bump (in)', 'FontSize', 20)
    ylabel('Front Left Caster (degrees)', 'FontSize', 20)
    title('Caster vs. Bump, Front Left', 'FontSize', 20)
end

function [x, y, z] = solveSpheres(p1, p2, p3, v1, v2, v3, p0)
    syms x y z
    a = (2*p2(1)-2*p1(1));
    b = (2*p2(2)-2*p1(2));
    c = (2*p2(3)-2*p1(3));
    d = v1^2-v2^2-p1(1)^2+p2(1)^2-p1(2)^2+p2(2)^2-p1(3)^2+p2(3)^2;
    A = (2*p3(1)-2*p1(1));
    B = (2*p3(2)-2*p1(2));
    C = (2*p3(3)-2*p1(3));
    D = v1^2-v3^2-p1(1)^2+p3(1)^2-p1(2)^2+p3(2)^2-p1(3)^2+p3(3)^2;
    A = [a b c d; A B C D];
    A = rref(A);
    eq1 = x*A(1,1)+y*A(1,2)+z*A(1,3) == A(1,4);
    eq2 = x*A(2,1)+y*A(2,2)+z*A(2,3) == A(2,4);
    x = solve(eq1, x);
    y = solve(eq2, y);
    eq3 = (x-p1(1))^2+(y-p1(2))^2+(z-p1(3))^2 == v1^2;
    zcoord = solve(eq3, z);
    z1 = double(zcoord(1));
    z2 = double(zcoord(2));
    x1 = double((A(1,4)-y*A(1,2)-z1*A(1,3))/A(1,1));
    x2 = double((A(1,4)-y*A(1,2)-z2*A(1,3))/A(1,1));
    y1 = double((A(2,4)-x*A(2,1)-z1*A(2,3))/A(2,2));
    y2 = double((A(2,4)-x*A(2,1)-z2*A(2,3))/A(2,2));
    if abs(norm([x1, y1, z1]-p0))<abs(norm([x2, y2, z2]-p0))
        x = x1;
        y = y1;
        z = z1;
    else
        x = x2;
        y = y2;
        z = z2;
    end
end