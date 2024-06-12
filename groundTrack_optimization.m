function [f1_merge, f1_diverge, f_formation, f2_merge, f2_diverge, f1_solo, f2_solo] = groundTrack_optimization(dep_coord, arr_coord, w1, w2, w12)

    f1_dep_lat = dep_coord(1);
    f1_dep_lon = dep_coord(2);
    f2_dep_lat = dep_coord(3);
    f2_dep_lon = dep_coord(4);
    
    f1_arr_lat = arr_coord(1);
    f1_arr_lon = arr_coord(2);
    f2_arr_lat = arr_coord(3);
    f2_arr_lon = arr_coord(4);
    
    % X and Y location calculation by similarity
    wgs84 = wgs84Ellipsoid("m"); % Earth Model in meters
    AB = distance(f1_dep_lat,f1_dep_lon, f2_dep_lat,f2_dep_lon, wgs84); % Calculate great circle distance between flights departures
    CD = distance(f1_arr_lat,f1_arr_lon, f2_arr_lat,f2_arr_lon, wgs84); % Calculate great circle distance between flights arrivals
    
    AX = w1/w12*AB/1000; % In km
    BX = w2/w12*AB/1000; % In km
    
    CY = w1/w12*CD/1000; % In km
    DY = w2/w12*CD/1000; % In km
    
    if dep_coord(1:2) == dep_coord(3:4)
        X_lat = dep_coord(1);
        X_lon = dep_coord(2);
    else
        [X_lat, X_lon] = scxsc(f1_dep_lat, f1_dep_lon, km2deg(AX), f2_dep_lat, f2_dep_lon, km2deg(BX), "degrees"); % Note there is 2 possible solutions
    end
    if arr_coord(1:2) == arr_coord(3:4)
        Y_lat = arr_coord(1);
        Y_lon = arr_coord(2);
    else
        [Y_lat, Y_lon] = scxsc(f1_arr_lat, f1_arr_lon, km2deg(CY), f2_arr_lat, f2_arr_lon, km2deg(DY), "degrees"); % Note there is 2 possible solutions
    end
    
    % Plot points for visual representation
    figure()
    geoplot(f1_dep_lat, f1_dep_lon, "rx", "MarkerSize", 10, "LineWidth", 2)
    hold on
    geoplot(f2_dep_lat, f2_dep_lon, "bx", "MarkerSize", 10, "LineWidth", 2)
    
    geoplot(f1_arr_lat, f1_arr_lon, "ro", "MarkerSize", 10, "LineWidth", 2)
    geoplot(f2_arr_lat, f2_arr_lon, "bo", "MarkerSize", 10, "LineWidth", 2)
    
    geoplot(X_lat, X_lon, "*", "MarkerSize", 10, "LineWidth", 2)
    geoplot(Y_lat, Y_lon, "*", "MarkerSize", 10, "LineWidth", 2)
    
    % Loop to fin the combination of back nodes X and Y from the 4 possibilities
    formation_segment = 0;
    for i = 1:length(X_lat)
        for j = 1:length(Y_lat)
            dist = distance(X_lat(i), X_lon(i), Y_lat(j), Y_lon(j), wgs84);
            if dist > formation_segment
                formation_segment = dist;
                i_index = i;
                j_index = j;
            end
        end
    end
    
    % Keep only the correct combination of X and Y for this case
    X_lat = X_lat(i_index);
    X_lon = X_lon(i_index);
    Y_lat = Y_lat(j_index);
    Y_lon = Y_lon(j_index);
    
    % Plot Similar triangle used to find the back vertices X and Y
    geoplot([f1_dep_lat f2_dep_lat X_lat f1_dep_lat], [f1_dep_lon f2_dep_lon X_lon f1_dep_lon])
    geoplot([f1_arr_lat f2_arr_lat Y_lat f1_arr_lat], [f1_arr_lon f2_arr_lon Y_lon f1_arr_lon])
    
    % Plot formation segment great circle defined by the back vertices X and Y
    [formation_gc_lat, formation_gc_lon] = track2(X_lat, X_lon, Y_lat, Y_lon, wgs84);
    geoplot(formation_gc_lat, formation_gc_lon)
    
    % Transform to Earth Centered Earth Fixed frame to find the circunscribed
    % sphere to the similar triangle in cartesian coordinates
    [Ax, Ay, Az] = geodetic2ecef(wgs84, f1_dep_lat, f1_dep_lon, 0);
    [Bx, By, Bz] = geodetic2ecef(wgs84, f2_dep_lat, f2_dep_lon, 0);
    [Xx, Xy, Xz] = geodetic2ecef(wgs84, X_lat, X_lon, 0);
    % Function call to find the circunscribed sphere
    [x, y, z, r] = earth_circ_sphere(Ax, Ay, Az, Bx, By, Bz, Xx, Xy, Xz);
    % Transform back to geodetic coordinates
    [lat0, lon0] = ecef2geodetic(wgs84, x, y, z);
    % Calculate and plot the small circle resultant from the intersection of
    % the sphere and the Earth's surface
    [lat,lon] = scircle1(lat0, lon0, r, [0 360], wgs84);
    geoplot(lat, lon)
    % Calculate the azimuth of the great circle between both back vertices
    az = azimuth(X_lat, X_lon, Y_lat, Y_lon, wgs84);
    % Find the intersection point P where the formation will merge/diverge
    [P_dep_lat,P_dep_lon] = gcxsc(X_lat, X_lon, az, lat0, lon0, km2deg(r/1000)); % Note that there is 2 possible intersections
    
    % Repeat the process for the opposite back vertex
    [Cx, Cy, Cz] = geodetic2ecef(wgs84, f1_arr_lat, f1_arr_lon, 0);
    [Dx, Dy, Dz] = geodetic2ecef(wgs84, f2_arr_lat, f2_arr_lon, 0);
    [Yx, Yy, Yz] = geodetic2ecef(wgs84, Y_lat, Y_lon, 0);
    [x, y, z, r] = earth_circ_sphere(Cx, Cy, Cz, Dx, Dy, Dz, Yx, Yy, Yz);
    [lat0, lon0] = ecef2geodetic(wgs84, x, y, z);
    [lat,lon] = scircle1(lat0, lon0, r, [0 360], wgs84);
    geoplot(lat, lon)
    az = azimuth(Y_lat, Y_lon, X_lat, X_lon, wgs84);
    [P_arr_lat,P_arr_lon] = gcxsc(Y_lat, Y_lon, az, lat0, lon0, km2deg(r/1000));
    
    % Given the 4 possible intersections, find out the correct one for this
    % case
    formation_segment = Inf;
    for i = 1:length(P_dep_lat)
        for j = 1:length(P_arr_lat)
            dist = distance(P_dep_lat(i), P_dep_lon(i), P_arr_lat(j), P_arr_lon(j), wgs84);
            if dist < formation_segment
                formation_segment = dist;
                i_index = i;
                j_index = j;
            end
        end
    end
    
    % Keep only the correct merging and diverging nodes
    P_dep_lat = P_dep_lat(i_index);
    P_dep_lon = P_dep_lon(i_index);
    P_arr_lat = P_arr_lat(j_index);
    P_arr_lon = P_arr_lon(j_index);
    
    % Extra comprobation to move the merging point to the closest departure
    % point if the distance is minimized
    [minDistance, idx] = min([distance(f1_dep_lat, f1_dep_lon ,P_arr_lat, P_arr_lon, wgs84), distance(f2_dep_lat, f2_dep_lon ,P_arr_lat, P_arr_lon, wgs84)]);
    if minDistance < formation_segment
        formation_segment = minDistance;
        if idx == 1
            P_dep_lat = f1_dep_lat;
            P_dep_lon = f1_dep_lon;
        else
            P_dep_lat = f2_dep_lat;
            P_dep_lon = f2_dep_lon;
        end
    end

    % Extra comprobation to move the diverging point to the closest arrival
    % point if the distance is minimized
    [minDistance, idx] = min([distance(f1_arr_lat, f1_arr_lon ,P_dep_lat, P_dep_lon, wgs84), distance(f2_arr_lat, f2_arr_lon ,P_dep_lat, P_dep_lon, wgs84)]);
    if minDistance < formation_segment
        formation_segment = minDistance;
        if idx == 1
            P_arr_lat = f1_arr_lat;
            P_arr_lon = f1_arr_lon;
        else
            P_arr_lat = f2_arr_lat;
            P_arr_lon = f2_arr_lon;
        end
    end
    
    % Plot the merging and diverging nodes
    geoplot([P_dep_lat, P_arr_lat], [P_dep_lon, P_arr_lon], "g^", "MarkerSize", 10, "LineWidth", 2)
    
    % Calculate the track points for the optimal formation trajectories
    optTrack_f1 = [track2(f1_dep_lat, f1_dep_lon, P_dep_lat, P_dep_lon, wgs84); ...
                   track2(P_dep_lat, P_dep_lon, P_arr_lat, P_arr_lon, wgs84); ...
                   track2(P_arr_lat, P_arr_lon, f1_arr_lat, f1_arr_lon, wgs84)];
    
    optTrack_f2 = [track2(f2_dep_lat, f2_dep_lon, P_dep_lat, P_dep_lon, wgs84); ...
                   track2(P_dep_lat, P_dep_lon, P_arr_lat, P_arr_lon, wgs84); ...
                   track2(P_arr_lat, P_arr_lon, f2_arr_lat, f2_arr_lon, wgs84)];
    
    % Calculate the track points for the optimal solo trajectories
    optTrack_f1_solo = track2(f1_dep_lat, f1_dep_lon, f1_arr_lat, f1_arr_lon, wgs84);
    
    optTrack_f2_solo = track2(f2_dep_lat, f2_dep_lon, f2_arr_lat, f2_arr_lon, wgs84);
    
    % Plot the optimized trajectories for the flights
    geoplot(optTrack_f1(:,1), optTrack_f1(:,2), "--r", "LineWidth", 2)
    geoplot(optTrack_f2(:,1), optTrack_f2(:,2), "--b", "LineWidth", 2)
    
    % Plot the optimized trajectories for the solo flights
    geoplot(optTrack_f1_solo(:,1), optTrack_f1_solo(:,2), ":r", "LineWidth", 2)
    geoplot(optTrack_f2_solo(:,1), optTrack_f2_solo(:,2), ":b", "LineWidth", 2)
    
    hold off
    
    % Calculate great circle distances covered by flight segment
    f1_merge = distance(f1_dep_lat, f1_dep_lon, P_dep_lat, P_dep_lon, wgs84);
    f1_diverge = distance(P_arr_lat, P_arr_lon, f1_arr_lat, f1_arr_lon, wgs84);
    f_formation = distance(P_dep_lat, P_dep_lon, P_arr_lat, P_arr_lon, wgs84);
    f2_merge = distance(f2_dep_lat, f2_dep_lon, P_dep_lat, P_dep_lon, wgs84);
    f2_diverge = distance(P_arr_lat, P_arr_lon, f2_arr_lat, f2_arr_lon, wgs84);

    % Calculate great circle distance for solo flights without formation
    f1_solo = distance(f1_dep_lat, f1_dep_lon, f1_arr_lat, f1_arr_lon, wgs84);
    f2_solo = distance(f2_dep_lat, f2_dep_lon, f2_arr_lat, f2_arr_lon, wgs84);

    close all
end