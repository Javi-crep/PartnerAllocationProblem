function [m_fuel, v] = segment_optimization_fixedH(AC_perf, f, R, W1, h0, h1)
    h = (h0 + h1)/2;
    height = h*unitsratio("meter", "feet"); % Mean height in meters
    [~,a,~,rho] = atmosisa(height);

    W0 = initialWeight_fixedH(AC_perf, f, R, rho, W1, a);

    m_fuel = (W0 - W1)/9.81;

    v = AC_perf(13)*a; % Average segment speed
end

