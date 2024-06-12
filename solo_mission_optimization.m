function [m_fuel, h, t] = solo_mission_optimization(R_solo, AC_perf)
    % Pre-allocation of variables
    h = zeros(1, 2);
    
    % Calculation of reserve fuel
    R = convvel(230, "kts", "m/s")*45*60; % Range covered by holding 45 mins at 230 KTAS
    [~,~,~,rho] = atmosisa(457.2); % Density at 1500 ft (Holding altitude)
    W1 = (AC_perf(3) + AC_perf(12))*9.81; % Empty weight of the AC + Payload
    W0 = initialWeight_fixedV(AC_perf, 1, R, rho, W1, convvel(230, "kts", "m/s"));
    W_reserve = (W0 - W1)*1.1; % +10% to account for the extra fuel on top of the extended mission fuel

    % Optimize solo mission
    W1 = W1 + W_reserve; % Final weight of the mission
    [m_fuel, h(1), h(2), v] = segment_optimization(AC_perf, 1, R_solo, W1);
    t = sqrt(R_solo^2 + (h(2)-h(1))^2)/v;

    h = h*unitsratio("feet", "meter"); % Height in feet
end