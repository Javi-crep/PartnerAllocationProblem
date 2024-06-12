function [x, y, z, r] = earth_circ_sphere(x1, y1, z1, x2, y2, z2, x3, y3, z3)
    R = earthRadius;
    F = @(x) [(x(1) - x1)^2 + (x(2) - y1)^2 + (x(3) - z1)^2 - x(4)^2; ...
              (x(1) - x2)^2 + (x(2) - y2)^2 + (x(3) - z2)^2 - x(4)^2; ...
              (x(1) - x3)^2 + (x(2) - y3)^2 + (x(3) - z3)^2 - x(4)^2; ...
               x(1)^2 + x(2)^2 + x(3)^2 - R^2];

    options = optimoptions('fsolve','Display','none');
    [sol, ~, exitFlag] = fsolve(F, [x1 y1 z1 1000], options);
    if exitFlag < 1
        if (x1 == x2) && (y1 == y2) && (z1 == z2)
            x = x1;
            y = y1;
            z = z1;
            r = 0;
        else
            error("Fsolve in earth_circ_sphere didn't converge or had some problem")
        end
    else
        x = sol(1);
        y = sol(2);
        z = sol(3);
        r = sol(4);
    end    
end

    
