function [m_fuel, h0, h1, v] = segment_optimization(AC_perf, f, R, W1)
    height = AC_perf(1)*100*unitsratio("meter", "feet"); % Initial height guess in meters
    [~,~,~,rho] = atmosisa(height);

    iterLoop = 1;
    while iterLoop
        W0 = initialWeight(AC_perf, f, R, rho, W1);
        [rho_update, h0, h1] = vertical_profile(W0, W1, AC_perf);
        if (rho - rho_update) > 1e-06
            rho = rho_update;
        else
            iterLoop = 0;
        end
    end

    m_fuel = (W0 - W1)/9.81;

    [~, a0] = atmosisa(h0);
    [~, a1] = atmosisa(h1);
    v = AC_perf(13)*(a0 + a1)/2; % Average segment speed
end

