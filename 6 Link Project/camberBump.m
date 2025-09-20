function camberBump(carParams, bumpRange)
    inboard = carParams.inboardR;
    outboard = carParams.outboardR;
    inboardU = [(inboard(1,2)+inboard(2,2))/2, (inboard(1,3)+inboard(2,3))/2];
    inboardL = [(inboard(3,2)+inboard(4,2))/2, (inboard(3,3)+inboard(4,3))/2];
    outboardU = [(outboard(1,2)+outboard(2,2))/2, (outboard(1,3)+outboard(2,3))/2];
    outboardL = [(outboard(3,2)+outboard(4,2))/2, (outboard(3,3)+outboard(4,3))/2];
    bvec0 = outboardL-inboardL;
    cvec = inboardU-inboardL;
    dvec0 = outboardU-inboardU;
    fvec0 = outboardU-outboardL;
    theta0 = -atand(fvec0(2)/fvec0(1));
    b = norm(bvec0);
    c = norm(cvec);
    d = norm(dvec0);
    f = norm(fvec0);
    camberValues = zeros(size(bumpRange, 1), 1);
    for i = 1:size(bumpRange, 1)
        bvec = [0, bvec0(2)-bumpRange(i)*25.4];
        bvec(1) = -sqrt(b^2-bvec(2)^2);
        A = acosd(dot(bvec, cvec)/(b*c));
        a = sqrt(b^2+c^2-2*b*c*cosd(A));
        alpha = acosd((a^2-d^2-f^2)/(-2*d*f));
        D = asind(d*sind(alpha)/a);
        avec = inboardU-(bvec+inboardL);
        theta = D-atand(avec(2)/avec(1));
        deltaTheta = theta-theta0;
        newCamber = -1.693+deltaTheta;
        camberValues(i) = newCamber;
    end
    plot(bumpRange, camberValues)
    xlabel('Rear Left Bump (in)', 'FontSize', 20)
    ylabel('Rear Left Camber Angle (deg)', 'FontSize', 20)
    title('Rear Left Camber Angle vs. Rear Left Bump', 'FontSize', 20)
end
    