function h0 = vertical_profile_fixedV(W0, W1, AC_perf, Vtas, rho)  
    S = AC_perf(4); % In m^2
    
    Cl = 2*W1/(rho*Vtas^2*S); % Average Cl of the segment 

    rho0 = 2*W0/(Cl*(Vtas)^2*S); % Air density at the end of the segment considering Cl constant
        
    rho_tropopause = 0.3639;
    if rho0 > rho_tropopause
        h0 = (1-(rho0/1.225)^(1/4.2259))/22.558e-06;
    else
        h0 = -log(rho0/0.3639)/157.69e-06 + 11000;
    end

    
end