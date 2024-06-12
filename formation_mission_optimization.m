function [m_fuel, h, t, v] = formation_mission_optimization(R_merge, R_formation, R_diverge, AC_perf, Vtas_merge)
    % Pre-allocation of variables
    m_fuel = zeros(1, 3);
    h = zeros(1, 4);
    t = zeros(1, 3);
    v = zeros(1, 3);
    
    % Calculation of reserve fuel
    R = convvel(230, "kts", "m/s")*45*60; % Range covered by holding 45 mins at 230 KTAS <  230 KEAS
    [~,~,~,rho] = atmosisa(457.2); % Density at 1500 ft (Holding altitude)
    W1 = (AC_perf(3) + AC_perf(12))*9.81; % Empty weight of the AC + Payload
    W0 = initialWeight_fixedV(AC_perf, 1, R, rho, W1, convvel(230, "kts", "m/s"));

    % Optimize diverge segment
    W1 = W0*1.1; % Final weight of the mission +10% to account for the extra fuel on top of the extended mission fuel
    [m_fuel(3), h(3), h(4), v(3)] = segment_optimization(AC_perf, 1, R_diverge, W1);
    t(3) = sqrt(R_diverge^2 + (h(4)-h(3))^2)/v(3);

    % Optimize formation segment
    W1 = W1 + m_fuel(3)*9.81; % Weight at divergence point
    [m_fuel(2), h(2), h(3), v(2)] = segment_optimization(AC_perf, 2, R_formation, W1);
    t(2) = sqrt(R_formation^2 + (h(3)-h(2))^2)/v(2);

    % Optimize formation segment
    W1 = W1 + m_fuel(2)*9.81; % Weight at divergence point
    [m_fuel(1), h(1)] = segment_optimization_fixedV(AC_perf, 1, R_merge, W1, Vtas_merge, h(2));
    t(1) = sqrt(R_merge^2 + (h(2)-h(1))^2)/Vtas_merge;
    v(1) = Vtas_merge;

    h = h*unitsratio("feet", "meter"); % Height in feet
end