function [rho, h0, h1] = vertical_profile(W0, W1, AC_perf)  

    S = AC_perf(4); % In m^2    
    Cd0 = AC_perf(7);
    b = AC_perf(10);

    e = 0.75; % Oswald factor
    k = 1/(pi*b^2/S*e);
    Cl = sqrt(Cd0/k)*1; % Cl Optimal

    M = AC_perf(13);
    height = AC_perf(1)*100*unitsratio("meter", "feet"); % Initial height guess in meters
    height0 = height;
    height1 = height;

    iterLoop0 = 1;
    iterLoop1 = 1;
    while (iterLoop0 + iterLoop1) ~= 0
        [~, a0] = atmosisa(height0);
        [~, a1] = atmosisa(height1);
        rho0 = 2*W0/(Cl*(M*a0)^2*S); % Air density at the beginning of the segment
        rho1 = 2*W1/(Cl*(M*a1)^2*S); % Air density at the end of the segment 
        
        rho_tropopause = 0.3639;
        if rho0 > rho_tropopause
            h0 = (1-(rho0/1.225)^(1/4.2259))/22.558e-06;
        else
            h0 = -log(rho0/0.3639)/157.69e-06 + 11000;
        end
        if rho1 > rho_tropopause
            h1 = (1-(rho1/1.225)^(1/4.2259))/22.558e-06;
        else
            h1 = -log(rho1/0.3639)/157.69e-06 + 11000;
        end

        if (h0 - height0) > 1e-06
            height0 = h0;
        else
            iterLoop0 = 0;
        end
        if (h1 - height1) > 1e-06
            height1 = h1;
        else
            iterLoop1 = 0;
        end
    end

    % Check that the altitude is within operating margin for the AC
    h_max = AC_perf(11)*unitsratio("meter", "feet");
    if h0 > h_max
        h0 = h_max;
        [~,~,~,rho0] = atmosisa(h0);
    end
    if h1 > h_max
        h1 = h_max;
        [~,~,~,rho1] = atmosisa(h1);
    end

    rho = (rho0 + rho1)/2;  % Mean air density for the segment

end