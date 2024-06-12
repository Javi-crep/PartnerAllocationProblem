function W0 = initialWeight_fixedV(AC_perf, f, R, rho, W1, Vtas)
    % Assumed overall savings given that follower AC get a 20% Saving
    switch f
        case 1
            lambda = 1;
        case 2
            lambda = 0.9;
        case 3
            lambda = 0.85;
    end
    
    S = AC_perf(4); % In m^2
    
    Cf1 = AC_perf(5)/(60*1000);
    Cf2 = convvel(AC_perf(6), "kts", "m/s");
    
    Cd0 = AC_perf(7);
    Cd2 = AC_perf(8);
    
    % AC performance calculations
    Ct = Cf1*(1 + Vtas/Cf2); % In Kg/(min*kN)
    Cl = 2*W1/(rho*Vtas^2*S); % Constant Cl assumed equal to the final Cl of the segment as a conservative approach
    Cd = Cd0 + Cd2*Cl^2;
    gamma = sqrt(rho*S/2)*Ct/(sqrt(Cl)/Cd);

    W0 = (sqrt(W1) + lambda*gamma/2*R)^2; % Final fuel calculation
end