function m_fuel = segment_optimization_fixedVH(AC_perf, f, R, W1, Vtas, h0, h1)
    h = (h0 + h1)/2;
    height = h*unitsratio("meter", "feet"); % Mean height in meters
    [~,~,~,rho] = atmosisa(height);

    W0 = initialWeight_fixedV(AC_perf, f, R, rho, W1, Vtas);

    m_fuel = (W0 - W1)/9.81;

end

