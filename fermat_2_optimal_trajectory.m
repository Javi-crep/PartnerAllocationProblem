% Fermat Problem for 2-formation flight
function optimalData = fermat_2_optimal_trajectory(dep_coord, arr_coord, AC1_perf, AC2_perf, lambda_f)

    % Initial speed guess in m/s
    Vtas_1 = convvel(AC1_perf(2), "kts", "m/s");
    Vtas_2 = convvel(AC2_perf(2), "kts", "m/s");
    
    % BADA Coefficients
    Cf1_1 = AC1_perf(5)/(60*1000);
    Cf2_1 = convvel(AC1_perf(6), "kts", "m/s");
    Cf1_2 = AC2_perf(5)/(60*1000);
    Cf2_2 = convvel(AC2_perf(6), "kts", "m/s");
    
    iterLoop0 = 1;
    iterLoop1 = 1;
    while (iterLoop0 + iterLoop1) ~= 0
        % Thrust specific fuel consumption calculation
        Ct_1 = Cf1_1*(1 + Vtas_1/Cf2_1);
        Ct_2 = Cf1_2*(1 + Vtas_2/Cf2_2);
        
        % Arc weights
        w1 = Ct_1; % A/C 1 arc weight (SFC)
        w2 = Ct_2; % A/C 2 arc weight (SFC)
        w12 = (w1 + w2)*lambda_f; % Formation segment arc weight
        
        % Geometrical calculation of the optimal ground track
        [f1_merge, f1_diverge, f_formation, f2_merge, f2_diverge, f1_solo, f2_solo] = groundTrack_optimization(dep_coord, arr_coord, w1, w2, w12);
        
        % Formation mission optimization
        [~, index] = max([f1_merge f2_merge]);
        if index == 1
            Vtas_merge_1 = Vtas_1;
            Vtas_merge_2 = Vtas_1*f2_merge/f1_merge;
        else
            Vtas_merge_1 = Vtas_2*f1_merge/f2_merge;
            Vtas_merge_2 = Vtas_2;
        end

        % Stall speeds in CAS we treat also as TAS
        Vstall_1 = convvel(AC1_perf(14), "kts", "m/s");
        Vstall_2 = convvel(AC2_perf(14), "kts", "m/s");

        if Vtas_1 <= Vstall_1 || Vtas_2 <= Vstall_2
            f1_m_fuel = NaN;
            f1_h = NaN;
            f1_t = NaN;
            f1_v = NaN;
            f2_m_fuel = NaN;
            f2_h = NaN;
            f2_t = NaN;
            f2_v = NaN;
            break;
        end
        
        [f1_m_fuel, f1_h, f1_t, f1_v] = formation_mission_optimization(f1_merge, f_formation, f1_diverge, AC1_perf, Vtas_merge_1);
        [f2_m_fuel, f2_h, f2_t, f2_v] = formation_mission_optimization(f2_merge, f_formation, f2_diverge, AC2_perf, Vtas_merge_2);

        % % See if the flights follow the same altitude
        if f1_h(2:3) ~= f2_h(2:3)
            h_formation = (f1_h*w1 + f2_h*w2)/(w1+w2);
            f1_h = h_formation;
            f2_h = h_formation;
            % Fix performance parameters equal for formation
            AC_perf_min = min(AC1_perf, AC2_perf);
            AC_perf_max = max(AC1_perf, AC2_perf);
        
            AC1_perf(3) = AC_perf_max(3); % Take maximum of minimum mass
            AC2_perf(3) = AC_perf_max(3);
        
            AC1_perf(12) = AC_perf_min(12); % Take minimum of payload mass
            AC2_perf(12) = AC_perf_min(12);
        
            AC1_perf(13) = AC_perf_min(13); % Take minimum Mach number
            AC2_perf(13) = AC_perf_min(13);
        
            % Optimized fixed altitude mission
            [f1_m_fuel, f1_t, f1_v] = formation_mission_optimization_fixedH(f1_merge, f_formation, f1_diverge, AC1_perf, h_formation, Vtas_merge_1);
            [f2_m_fuel, f2_t, f2_v] = formation_mission_optimization_fixedH(f2_merge, f_formation, f2_diverge, AC2_perf, h_formation, Vtas_merge_2);
        end
    
        if (Vtas_1 - f1_v(1)) > 1e-06
                Vtas_1 = f1_v(1);
        else
                iterLoop0 = 0;
        end
        if (Vtas_2 - f2_v(1)) > 1e-06
                Vtas_2 = f2_v(1);
        else
                iterLoop1 = 0;
        end
    end
    
    % Solo mission optimization
    [f1_solo_m_fuel, f1_solo_h, f1_solo_t] = solo_mission_optimization(f1_solo, AC1_perf);
    [f2_solo_m_fuel, f2_solo_h, f2_solo_t] = solo_mission_optimization(f2_solo, AC2_perf);

    % Save results as structure
    f1_range = [f1_merge f_formation f1_diverge];
    f2_range = [f2_merge f_formation f2_diverge];
    f1_solo_range = f1_solo;
    f2_solo_range = f2_solo;
    optimalData = struct( ...
        "f1_fuel", f1_m_fuel, ...
        "f2_fuel", f2_m_fuel, ...
        "f1_solo_fuel", f1_solo_m_fuel, ...
        "f2_solo_fuel", f2_solo_m_fuel, ...
        "f1_altitude", f1_h, ...
        "f2_altitude", f2_h, ...
        "f1_solo_altitude", f1_solo_h, ...
        "f2_solo_altitude", f2_solo_h, ...
        "f1_range", f1_range, ...
        "f2_range", f2_range, ...
        "f1_solo_range", f1_solo_range, ...
        "f2_solo_range", f2_solo_range, ...
        "f1_time", f1_t/60, ...
        "f2_time", f2_t/60, ...
        "f1_solo_time", f1_solo_t/60, ...
        "f2_solo_time", f2_solo_t/60, ...
        "f1_speed", f1_v, ...
        "f2_speed", f2_v, ...
        "f1_solo_speed", f1_solo_range/f1_solo_t, ...
        "f2_solo_speed", f2_solo_range/f2_solo_t);
end