function [m_fuel, h0] = segment_optimization_fixedV(AC_perf, f, R, W1, Vtas, h1)
    [~,~,~,rho] = atmosisa(h1);

    W0 = initialWeight_fixedV(AC_perf, f, R, rho, W1, Vtas);
    h0 = vertical_profile_fixedV(W0, W1, AC_perf, Vtas, rho);

    m_fuel = (W0 - W1)/9.81;

end


