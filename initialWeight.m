function W0 = initialWeight(AC_perf, f, R, rho, W1)
    % Assumed overall savings given that follower AC get a 20% Saving
    switch f
        case 1
            lambda = 1;
        case 2
            lambda = 0.9;
        case 3
            lambda = 0.85;
    end

    % Inputs from plane type BADA
     rho_tropopause = 0.3639;
        if rho > rho_tropopause
            h = (1-(rho/1.225)^(1/4.2259))/22.558e-06;
        else
            h = -log(rho/0.3639)/157.69e-06 + 11000;
        end
    [~, a] = atmosisa(h);
    Vtas = AC_perf(13)*a; % in m/s
    
    S = AC_perf(4); % In m^2
    
    Cf1 = AC_perf(5)/(60*1000);
    Cf2 = convvel(AC_perf(6), "kts", "m/s");
    
    Cd0 = AC_perf(7);
    Cd2 = AC_perf(8);

    b = AC_perf(10);

    e = 0.75;

    % Environment constants
    g = 9.81;
    
    % AC performance calculations
    Ct = Cf1*(1 + Vtas/Cf2); % In Kg/(min*kN)
    k = 1/(pi*b^2/S*e);
    Cl = sqrt(Cd0/k)*1; % Cl Optimal
    Cd = Cd0 + Cd2*Cl^2;
    gamma = sqrt(rho*S/2)*Ct/(sqrt(Cl)/Cd);

    W0 = (sqrt(W1) + lambda*gamma/2*R)^2; % Final fuel calculation
end