function posForms = possible_formations(flight, compSet, n, nRows, t_comp_index)
    % Data allocated in matrix "compSet" and "flight" for convenience 
    % rows --> One for each flight in the set
    % columns --> [1] = t1 latitude
    %             [2] = t2 longitude
    %             [3] = t1 latitude
    %             [4] = t2 longitude
    %             [5] = bound
    %             [6] = t1 azimuth
    %             [7] = t2 azimuth
    %             [8] = t1 FL
    %             [9] = t2 FL
    
    posForms = zeros(nRows, 1);

    % East and West bound filter
    index_bound = (compSet(:,5) == flight(:,5));

    % Departure and Destination Azimuth filter
    aTol_dep = 45; % Angular tolerance in degrees at departure
    aTol_arr = 2.5; % Angular tolerance in degrees at arrival
        % Azimuth filter at t1        
    index_azimuthDEP = ((compSet(:,6) >= (flight(:,6) - aTol_dep)) .* (compSet(:,6) <= (flight(:,6) + aTol_dep)));
        % Azimuth filter at t2
    index_azimuthDES = ((compSet(:,7) >= (flight(:,7) - aTol_arr)) .* (compSet(:,7) <= (flight(:,7) + aTol_arr)));
        % Azimuth index
    index_azimuth = index_azimuthDEP .* index_azimuthDES;

    % Distance filter
    dTol_dep = 5; % Angular tolerance in degrees at departure
    dTol_arr = 2.5; % Angular tolerance in degrees at arrival
        % Separation filter at t1
    dep_distance = sqrt((flight(:,1) - compSet(:,1)).^2 + (flight(:,2) - compSet(:,2)).^2);
        % Separation filter at t2
    arr_distance = sqrt((flight(:,3) - compSet(:,3)).^2 + (flight(:,4) - compSet(:,4)).^2);
    %    % Fancy distance filter
    % theta = flight(:,6);
    % rot_lat = (compSet(:,1) - flight(:,1)).*cosd(theta) + (compSet(:,2) - flight(:,2)).*cosd(theta);
    % rot_lon = -(compSet(:,1) - flight(:,1)).*sind(theta) + (compSet(:,2) - flight(:,2)).*cosd(theta);
    % 
    % b = 5;
    % m = 100;
    % 
    % if rot_lon >= 0 
    %     idx1 = (rot_lat <= b+m*rot_lon);
    %     idx2 = (rot_lat >= -b-m*rot_lon);
    % else
    %     idx1 = (rot_lat <= b-m*rot_lon);
    %     idx2 = (rot_lat >= -b+m*rot_lon);
    % end
    % index_separation_fancy = idx1.*idx2;
        % Separation index
    index_separation = (dep_distance <= dTol_dep) .* (arr_distance <= dTol_arr);

   
    kT = 60; % At least 60 mins of formation flight
    t_overlap = t_comp_index(:,2) - t_comp_index(:,1);
    index_overlap = (t_overlap >= kT); 
    % Filter now accounts for time of formation as absolute value, do we
    % prefer a percentage of the total flight?

    % Flight level filter
    % hTol = 2000; % Flight level tolerance in feet
    %     % Flight level mean difference
    % fl_mean_diff = (flight(:,8) + flight(:,9))/2 - (compSet(:,8) + compSet(:,9))/2;
    %     % Flight level index
    % index_flightLevel = (abs(fl_mean_diff) <= hTol);
    % 
    %     % Flight level filter at t1
    % index_FL_DEP = ((compSet(:,8) >= (flight(:,8) - hTol)) .* (compSet(:,8) <= (flight(:,8) + hTol)));
    %     % Flight level filter at t2
    % index_FL_DES = ((compSet(:,9) >= (flight(:,9) - hTol)) .* (compSet(:,9) <= (flight(:,9) + hTol)));
    %     % Flight level index
    % index_flightLevel = index_FL_DEP .* index_FL_DES;  

    % Assignment
    posForms(n+1:nRows) = index_bound .* index_azimuth .* index_separation .* index_overlap;

end