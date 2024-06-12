clear; clc; close all; warning off

tic %%%


rawData = readtable("EUROCONTROL\FlightDataJUN2019.xlsx", "Sheet", "DataSet0"); % Baseline One World Scenario
% rawData = readtable("EUROCONTROL\FlightDataJUN2019.xlsx", "Sheet", "DataSet1"); % USA Carriers Scenario
% rawData = readtable("EUROCONTROL\FlightDataJUN2019.xlsx", "Sheet", "DataSet2"); % Low Cost Scenario
% rawData = readtable("EUROCONTROL\FlightDataJUN2019.xlsx", "Sheet", "DataSet3"); % All together Scenario

fprintf("\n--- Data exported from excel ---\n")
toc %%%

%% MAP PLOT OF ALL FLIGHTS
tic %%%

figure(1)
gx = geoaxes; % Use geographical coordinate axis

% window = [700 800]; % Flight departure window
% index = logical((data.TDEP >= window(1)).*(data.TDEP < window(2))); % Filter departure range window
% data = data(index, :); % Filter data values of interest

[nRows, ~] = size(rawData); 
depLat = rawData.ADEPLatitude;
depLon = rawData.ADEPLongitude;
desLat = rawData.ADESLatitude;
desLon = rawData.ADESLongitude;
color = rawData.Color;

nPointsMax = 48*60; % For maximum operation window of 48h with precision by minute (this could be smaller)
gC_path_lat = NaN(nRows, nPointsMax);
gC_path_lon = NaN(nRows, nPointsMax);

tDep = rawData.FILEDOFFBLOCKTIME;
% tDep = rawData.ACTUALOFFBLOCKTIME;
tDep_min = ((day(tDep)-1)*24 + hour(tDep))*60 + minute(tDep);
tArr = rawData.FILEDARRIVALTIME;
% tArr = rawData.ACTUALARRIVALTIME;
tArr_min = ((day(tArr)-1)*24 + hour(tArr))*60 + minute(tArr);
fTime = (tArr_min - tDep_min); % Flight time in minutes

nPoints = fTime + 1; % One track point every minute of flight

wgs84 = wgs84Ellipsoid("m"); % Earth Model in meters

for i = 1:nRows
    [gC_path_lat(i,tDep_min(i):tArr_min(i)),gC_path_lon(i,tDep_min(i):tArr_min(i))] = track2(depLat(i), depLon(i), desLat(i), desLon(i), wgs84, "degrees", nPoints(i)); % Calculate 100 track points along great circle path
end

plot(gx, gC_path_lat', gC_path_lon') 
title("OneWorld Alliance flights in June 1^{st} 2019")
colororder(color)

fprintf("\n--- Map plot and calculation of the great circle trajectories ---\n")
toc %%%
%% HEURISTIC FILTER LOOP
tic %%%

posibleFormations = zeros(nRows);

% Correction for initial climb
initialClimbTimes = [15 25 30]; % Average initial climb margins for SH, MH, LH flights in minutes

sh_index = (fTime/60 < 2);
mh_index = (fTime/60 <= 4) .* (fTime/60 >= 2);
lh_index= (fTime/60 > 4);

initialClimbMargin = sh_index*initialClimbTimes(1) + mh_index*initialClimbTimes(2) + lh_index*initialClimbTimes(3);
tDep_min = tDep_min + initialClimbMargin;

% Correction for final descent and landing
finalDescentTimes = [15 25 35]; % Average final descent margins for SH, MH, LH flights in minutes

finalDescentMargin = sh_index*finalDescentTimes(1) + mh_index*finalDescentTimes(2) + lh_index*finalDescentTimes(3);
tArr_min = tArr_min - finalDescentMargin;

for i = 1:nRows-1
    t_comp_index = [max(tDep_min(i), tDep_min(i+1:end)), min(tArr_min(i), tArr_min(i+1:end))];
    % Still need to add the margin time for takeoff and landing
    
    % Position coordinate extraction for flight to compare to
    dep_flight_lat = gC_path_lat(i, t_comp_index(:,1));
    dep_flight_lon = gC_path_lon(i, t_comp_index(:,1));
    arr_flight_lat = gC_path_lat(i, t_comp_index(:,2));
    arr_flight_lon = gC_path_lon(i, t_comp_index(:,2));

    bound_flight = sign(arr_flight_lon - dep_flight_lon); % 1 for East and -1 for West
    dep_flight_azimuth = azimuth(dep_flight_lat, dep_flight_lon, arr_flight_lat, arr_flight_lon, wgs84, "degrees");
    arr_flight_azimuth = azimuth(arr_flight_lat, arr_flight_lon, dep_flight_lat, dep_flight_lon, wgs84, "degrees");
    
    % Data allocated in matrix for convenience [departure latitude,
    % departure longitude, arrival latitude, arrival longitude, bound,
    % departure azimuth, arrival azimuth]
    flight = [dep_flight_lat', dep_flight_lon' , arr_flight_lat', arr_flight_lon', ...
        bound_flight', dep_flight_azimuth', arr_flight_azimuth'];

    % Position coordinate extraction for flights in the comparison set
    dep_compSet_lat = gC_path_lat(sub2ind(size(gC_path_lat), i+1:nRows, t_comp_index(:,1)'));
    dep_compSet_lon = gC_path_lon(sub2ind(size(gC_path_lon), i+1:nRows, t_comp_index(:,1)'));
    arr_compSet_lat = gC_path_lat(sub2ind(size(gC_path_lat), i+1:nRows, t_comp_index(:,2)'));
    arr_compSet_lon = gC_path_lon(sub2ind(size(gC_path_lon), i+1:nRows, t_comp_index(:,2)'));

    bound_compSet = sign(arr_compSet_lon - dep_compSet_lon); % 1 for East and -1 for West
    dep_compSet_azimuth = azimuth(dep_compSet_lat, dep_compSet_lon, arr_compSet_lat, arr_compSet_lon, wgs84, "degrees");
    arr_compSet_azimuth = azimuth(arr_compSet_lat, arr_compSet_lon, dep_compSet_lat, dep_compSet_lon, wgs84, "degrees");

    % Data allocated in matrix for convenience [departure latitude,
    % departure longitude, arrival latitude, arrival longitude, bound,
    % departure azimuth, arrival azimuth]
    compSet = [dep_compSet_lat', dep_compSet_lon', arr_compSet_lat', arr_compSet_lon', ...
        bound_compSet', dep_compSet_azimuth', arr_compSet_azimuth'];

    posibleFormations(:,i) = possible_formations(flight, compSet, i, nRows, t_comp_index);
end

figure(2)
G = graph(posibleFormations, 'lower');
p = plot(G, "NodeLabel", 1:nRows);
title("Graph of inter-flight compatibility")

fprintf("\n--- Heuristic filters applied ---\n")
toc %%%

%% MAP PLOT OF SELECTED FLIGHT POSIBLE FORMATIONS
disp("--------------------------------------------")
target = 1;
while target ~= 0

    target = input("\nEnter the number ID of the flight to plot the posible formations [0 to skip loop]: ");
    if target == 0
        break
    end

    index = [target; neighbors(G,target)];

    tic %%%

    figure(3)
    clf % Clear figure
    gx = geoaxes; % Use geographical coordinate axis
    plot(gx, gC_path_lat(index,:)', gC_path_lon(index,:)') 

    title("OneWorld Alliance Selected Flights in June 1^{st} 2019")
    legend(index+ " || " + rawData.ADEP(index) + "-" + rawData.ADES(index) +  " || " + rawData.ACOperator(index) + " || " + rawData.ACType(index) + " || DEP:" + string(rawData.FILEDOFFBLOCKTIME(index)) + " || ARR:" + string(rawData.FILEDARRIVALTIME(index)))

    highlight(p, index, "NodeColor", "g", "MarkerSize", 4)
    highlight(p, target, "NodeColor", "r", "MarkerSize", 4)

    toc %%%

end

%% CASE COMBINATION IDENTIFICATION
disp("--------------------------------------------")
tic %%%

% Subgraph count and size
[bins, binsizes] = conncomp(G);
nSubGraphs = max(bins); % Number of subgraphs

subgraphsOfTwo = nonzeros((1:max(bins)).*(binsizes == 2)); % Vector with the ID of the subgraphs that have 2 nodes
subgraphsOfThree = nonzeros((1:max(bins)).*(binsizes == 3)); % Vector with the ID of the subgraphs that have 3 nodes
subgraphsOfMore = nonzeros((1:max(bins)).*(binsizes > 3)); % Vector with the ID of the subgraphs that have more than 3 nodes

flightID_2forms = (bins == subgraphsOfTwo).*(1:nRows); % Matrix with flight IDs 
flightID_3forms = (bins == subgraphsOfThree).*(1:nRows);
flightID_4Plusforms = (bins == subgraphsOfMore).*(1:nRows);

fprintf("\nThere are the following number of subgraphs:" + ...
    " \nSubgraphs of [1] = %i \nSubgraphs of [2] = %i \nSubgraphs of [3] = %i \nSubgraphs of [4+] = %i \nTotal subgraphs = %i\n", ...
    nSubGraphs-(length(subgraphsOfMore)+length(subgraphsOfThree)+length(subgraphsOfTwo)), length(subgraphsOfTwo), length(subgraphsOfThree), length(subgraphsOfMore), nSubGraphs);

% Case for 2-A/C-formations
two_formations = zeros(length(subgraphsOfTwo), 2); % Pre-allocation
    % Flght ID extraction of subgraphs of 2 nodes
for i = 1:length(subgraphsOfTwo)
        two_formations(i,:) = nonzeros(flightID_2forms(i,:));
end

restOfCases = [flightID_3forms;flightID_4Plusforms];
for j = 1:(size(restOfCases, 1))
    subgraphIDs = nonzeros(restOfCases(j,:));
    s = subgraph(G, subgraphIDs);
    node_degrees = degree(s);
    [max_degree, cutNode] = max(node_degrees);
    while max_degree > 0
        partners = neighbors(s,cutNode);
        combinations = [ones(size(partners))*subgraphIDs(cutNode), subgraphIDs(partners)];
        two_formations = [two_formations;combinations];
        s = rmnode(s, cutNode);
        subgraphIDs(cutNode) = [];
        node_degrees = degree(s);
        [max_degree, cutNode] = max(node_degrees);
    end
end

% Case for 3-A/C-formations
three_formations = cell2mat(allcycles(G, "MinCycleLength", 3,'MaxCycleLength', 3));
fprintf("\n--- Possible formation combinations created ---\n")
toc %%%

% Command Window Output data on filters
total_combinations = nRows + length(two_formations) + length(three_formations);

nofilter_2combinations = nchoosek(nRows, 2);
nofilter_3combinations = nchoosek(nRows, 3);
nofilter_total_combinations = nRows + nofilter_2combinations + nofilter_3combinations;

filter_ratio_1 = (1-nRows/nRows)*100;
filter_ratio_2 = (1-length(two_formations)/nofilter_2combinations)*100;
filter_ratio_3 = (1-length(three_formations)/nofilter_3combinations)*100;
filter_ratio_total = (1-total_combinations/nofilter_total_combinations)*100;

fprintf("\nTotal number of formation combinations after filters:\n Individual flights = %i [%i Originally] Filtered ratio = %.2f\n" + ...
    " 2-Formations = %i [%i Originally] Filtered ratio = %.2f\n" + ...
    " 3-Formations = %i [%i Originally] Filtered ratio = %.2f\n" + ...
    " Total Combinations = %i [%i Originally] Filtered ratio = %.2f\n", ...
    nRows, nRows, filter_ratio_1, ...
    length(two_formations), nofilter_2combinations, filter_ratio_2 , ...
    length(three_formations), nofilter_3combinations, filter_ratio_3, ...
    total_combinations, nofilter_total_combinations, filter_ratio_total);

% Plot with the final combinations to be computed

flightsIn2 = unique(two_formations);
flightsIn3 = unique(three_formations);


figure(4)
gx = geoaxes; % Use geographical coordinate axis

f1 = plot(gx, gC_path_lat', gC_path_lon', "k");
hold on
f2 = plot(gx, gC_path_lat(flightsIn2,:)', gC_path_lon(flightsIn2,:)', "Color", "#77AC30");
hold off

title("Possible flight formations after filters")
legend([f1(1), f2(1)], {"Individual flights", "Flights in 2-formation"})

%% GEOMETRIC TRAJECTORY OPTIMIZATION 2-FORMATIONS
disp("--------------------------------------------")
tic %%%

data2Formations = zeros(length(two_formations), 4);
plane2Formations = string(zeros(length(two_formations), 2));
time2Formations = zeros(length(two_formations), 4);

perf_data = readtable("BADA_data.xlsx");

parfor k = 1:length(two_formations)
    warning off
    % Extract index of flights
    f1 = two_formations(k, 1); % Flight one identifier
    f2 = two_formations(k, 2); % Flight two identifier

    % Flight 1 initial and final cruise locations
    f1_dep_lat = gC_path_lat(f1,tDep_min(f1));
    f1_dep_lon = gC_path_lon(f1,tDep_min(f1));
    f1_arr_lat = gC_path_lat(f1,tArr_min(f1));    
    f1_arr_lon = gC_path_lon(f1,tArr_min(f1));

    % Flight 2 initial and final cruise locations
    f2_dep_lat = gC_path_lat(f2,tDep_min(f2));
    f2_dep_lon = gC_path_lon(f2,tDep_min(f2));
    f2_arr_lat = gC_path_lat(f2,tArr_min(f2));    
    f2_arr_lon = gC_path_lon(f2,tArr_min(f2));

    % Organize coordinates in vector for optimization
    dep_coord = [f1_dep_lat, f1_dep_lon, f2_dep_lat, f2_dep_lon];
    arr_coord = [f1_arr_lat, f1_arr_lon, f2_arr_lat, f2_arr_lon];

    % Import performace data
    f1_planeType = rawData.PlaneType{f1};
    f2_planeType = rawData.PlaneType{f2};

    AC1_perf = perf_data.(f1_planeType);
    AC2_perf = perf_data.(f2_planeType);
    
    optData = fermat_2_optimal_trajectory(dep_coord, arr_coord, AC1_perf, AC2_perf, 0.9);

    % Total fuel burn and time for comparison solo vs formation
    f_fuel_formation = sum(optData.f1_fuel + optData.f2_fuel);
    f_fuel_solo = sum(optData.f1_solo_fuel + optData.f2_solo_fuel);
    f_time_formation = sum(optData.f1_time + optData.f2_time);
    f_time_solo = sum(optData.f1_solo_time + optData.f2_solo_time);

    data2Formations(k,:) = [f_fuel_formation, f_fuel_solo, f_time_formation, f_time_solo];
    plane2Formations(k,:) = [string(f1_planeType) string(f2_planeType)];
    time2Formations(k,:) = [sum(optData.f1_time), sum(optData.f2_time), sum(optData.f1_solo_time), sum(optData.f2_solo_time)];

    % % Visualization of results
    % figure()
    % subplot(2,2,1)
    % plot(optData.f1_altitude, "--*r")
    % hold on
    % plot(optData.f2_altitude, "--ob")
    % plot([1 4], optData.f1_solo_altitude, "-*r")
    % plot([1 4], optData.f2_solo_altitude, "-ob")
    % xticks([1 2 3 4])
    % xticklabels(["Cruise Start", "Formation Start", "Formation End", "Cruise End"])
    % title("Vertical optimal profiles")
    % legend(["F1 in formation", "F2 in formation", "F1 solo", "F2 solo"])
    % xlabel("Milestone")
    % ylabel("Altitude [ft]")
    % hold off
    % 
    % subplot(2,2,2)
    % plot(optData.f1_speed, "--*r")
    % hold on
    % plot(optData.f2_speed, "--ob")
    % plot([1 3], [1 1]*optData.f1_solo_speed, "-*r")
    % plot([1 3], [1 1]*optData.f2_solo_speed, "-ob")
    % xticks([1 2 3])
    % xticklabels(["Merging", "Formation", "Diverging"])
    % title("Speed optimal profiles")
    % legend(["F1 in formation", "F2 in formation", "F1 solo", "F2 solo"])
    % xlabel("Milestone")
    % ylabel("Speed [m/s]")
    % hold off
    % 
    % subplot(2,2,3)
    % b = bar(["Flight 1 in Formation", "Flight 2 in Formation"],[optData.f1_fuel; optData.f2_fuel], "stacked");
    % hold on
    % b_solo = bar(["Flight 1 Solo", "Flight 2 Solo"],[optData.f1_solo_fuel; optData.f2_solo_fuel], "stacked");
    % b_total = bar(["Formation", "Solo"],[f_fuel_formation; f_fuel_solo], "stacked");
    % title("Fuel Consumption")
    % ylabel("Fuel consumed [kg]")
    % hold off
    % 
    % p = [b b_solo b_total];
    % for j = 1:length(p)
    %     b = p(j);
    %     for i = 1:length(b)
    %     xtips1 = b(i).XEndPoints;
    %     ytips1 = b(i).YEndPoints;
    %     labels1 = string(b(i).YData);
    %     text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    %         'VerticalAlignment','cap')
    %     end
    % end
    % 
    % subplot(2,2,4)
    % b = bar(["Flight 1 in Formation", "Flight 2 in Formation"],[optData.f1_time; optData.f2_time], "stacked");
    % hold on
    % b_solo = bar(["Flight 1 Solo", "Flight 2 Solo"],[optData.f1_solo_time; optData.f2_solo_time], "stacked");
    % b_total = bar(["Formation", "Solo"],[f_time_formation; f_time_solo], "stacked");
    % title("Flight Duration")
    % ylabel("Time [min]")
    % hold off
    % 
    % p = [b b_solo b_total];
    % for j = 1:length(p)
    %     b = p(j);
    %     for i = 1:length(b)
    %     xtips1 = b(i).XEndPoints;
    %     ytips1 = b(i).YEndPoints;
    %     labels1 = string(b(i).YData);
    %     text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    %         'VerticalAlignment','cap')
    %     end
    % end
    % close all
end

% clearvars -except data2Formations % Delete all varibales except

% Determine operating cost valid formations
fuel_cost = 0.8; % In euros per kg of kerosene
crew_cost_L = 24; % In euros per hour
crew_cost_AM = 52; % In euros per hour
crew_cost_AH = 36; % In euros per hour
crew_cost_AJ = 40; % In euros per hour
crew_cost_BM = 16; % In euros per hour
crew_cost_BH = 72; % In euros per hour
crew_cost_BJ = 44; % In euros per hour

crew_cost = zeros(size(time2Formations));
for i = 1:length(plane2Formations)
    for j = 1:2
        switch plane2Formations(i,j)
            case "L"
            crew_cost(i, j) = crew_cost_L;
            crew_cost(i, j+2) = crew_cost_L;
            case "AM"
            crew_cost(i, j) = crew_cost_AM;
            crew_cost(i, j+2) = crew_cost_AM;
            case "AH"
            crew_cost(i, j) = crew_cost_AH;
            crew_cost(i, j+2) = crew_cost_AH;
            case "AJ"
            crew_cost(i, j) = crew_cost_AJ;
            crew_cost(i, j+2) = crew_cost_AJ;
            case "BM"
            crew_cost(i, j) = crew_cost_BM;
            crew_cost(i, j+2) = crew_cost_BM;
            case "BH"
            crew_cost(i, j) = crew_cost_BH;
            crew_cost(i, j+2) = crew_cost_BH;
            case "BJ"
            crew_cost(i, j) = crew_cost_BJ;
            crew_cost(i, j+2) = crew_cost_BJ;
        end
    end
end

total_cost_formation = fuel_cost*data2Formations(:,1) + sum(crew_cost(:,1:2).*time2Formations(:,1:2)/60, 2);
total_cost_formation_solo = fuel_cost*data2Formations(:,2) + sum(crew_cost(:,3:4).*time2Formations(:,3:4)/60, 2);
valid2Formations_TOC = total_cost_formation < total_cost_formation_solo;

% n_of_valid2Formations = sum(valid2Formations_TOC);
% efficiency_2Formations = n_of_valid2Formations/numel(valid2Formations_TOC)*100; % Efficiency of the filters for the optimization
% cost_saved_2Formations = (total_cost_solo(valid2Formations_TOC) - total_cost_formation(valid2Formations_TOC))./ total_cost_solo(valid2Formations_TOC); % TOC saved as percentage of solo TOC
% average_TOC_saved_2Formations = sum(cost_saved_2Formations)/numel(cost_saved_2Formations)*100;
% fuel_saved_2Formations = (data2Formations(valid2Formations_TOC, 2) - data2Formations(valid2Formations_TOC, 1))./data2Formations(valid2Formations_TOC, 2); % Fuel saved as percentage of solo fuel
% average_fuel_saved_2Formations = sum(fuel_saved_2Formations)/numel(fuel_saved_2Formations)*100;
% time_saved_2Formations = (data2Formations(valid2Formations_TOC, 4) - data2Formations(valid2Formations_TOC, 3))./data2Formations(valid2Formations_TOC, 4); % time saved as percentage of solo time
% average_time_saved_2Formations = sum(time_saved_2Formations)/numel(time_saved_2Formations)*100;
% 
% fprintf("\n--- Total Operating Costs Criteria ---\n")
% fprintf("There is a total of %d valid formations out of %d [Filter Efficiency = %.2f%%] \n", n_of_valid2Formations, numel(valid2Formations_TOC), efficiency_2Formations)
% fprintf("The TOC saving by the 2-formations is a %.2f%% of the TOC in solo mission\n", average_TOC_saved_2Formations)
% fprintf("The fuel saving by the 2-formations is a %.2f%% of the fuel in solo mission\n", average_fuel_saved_2Formations)
% fprintf("The time saving by the 2-formations is a %.2f%% of the time in solo mission\n", average_time_saved_2Formations)
% 
% % Table organization of the results
% results2Formation_TOC = table(identifier2Formations(valid2Formations_TOC,:), plane2Formations(valid2Formations_TOC,:), data2Formations(valid2Formations_TOC,:), cost_saved_2Formations*100, fuel_saved_2Formations*100, time_saved_2Formations*100);
% results2Formation_TOC = splitvars(results2Formation_TOC);
% results2Formation_TOC = renamevars(results2Formation_TOC, ["Var1_1", "Var1_2", "Var2_1", "Var2_2", "Var3_1", "Var3_2", "Var3_3", "Var3_4", "Var4", "Var5", "Var6"], ...
%     ["f1_ID", "f2_ID", "f1_type", "f2_type", "fuel_formation", "fuel_solo", "time_formation", "time_solo", "TOC_saved", "fuel_saved", "time_saved"]);


% Determine fuel/emission valid formations
valid2Formations_CO2 = data2Formations(:,1) < data2Formations(:,2);

% n_of_valid2Formations = sum(valid2Formations_CO2);
% efficiency_2Formations = n_of_valid2Formations/numel(valid2Formations_CO2)*100; % Efficiency of the filters for the optimization
% fuel_saved_2Formations = (data2Formations(valid2Formations_CO2, 2) - data2Formations(valid2Formations_CO2, 1))./data2Formations(valid2Formations_CO2, 2); % Fuel saved as percentage of solo fuel
% average_fuel_saved_2Formations = sum(fuel_saved_2Formations)/numel(fuel_saved_2Formations)*100;
% time_saved_2Formations = (data2Formations(valid2Formations_CO2, 4) - data2Formations(valid2Formations_CO2, 3))./data2Formations(valid2Formations_CO2, 4); % time saved as percentage of solo time
% average_time_saved_2Formations = sum(time_saved_2Formations)/numel(time_saved_2Formations)*100;
% 
% fprintf("\n--- Fuel consumption and CO2 emission Criteria ---\n")
% fprintf("There is a total of %d valid formations out of %d [Filter Efficiency = %.2f%%] \n", n_of_valid2Formations, numel(valid2Formations_CO2), efficiency_2Formations)
% fprintf("The fuel saving by the 2-formations is a %.2f%% of the fuel in solo mission\n", average_fuel_saved_2Formations)
% fprintf("The time saving by the 2-formations is a %.2f%% of the time in solo mission\n", average_time_saved_2Formations)
% 
% % Table organization of the results
% results2Formation_CO2 = table(identifier2Formations(valid2Formations_CO2,:), plane2Formations(valid2Formations_CO2,:), data2Formations(valid2Formations_CO2,:), fuel_saved_2Formations*100, time_saved_2Formations*100);
% results2Formation_CO2 = splitvars(results2Formation_CO2);
% results2Formation_CO2 = renamevars(results2Formation_CO2, ["Var1_1", "Var1_2", "Var2_1", "Var2_2", "Var3_1", "Var3_2", "Var3_3", "Var3_4", "Var4", "Var5"], ...
%     ["f1_ID", "f2_ID", "f1_type", "f2_type", "fuel_formation", "fuel_solo", "time_formation", "time_solo", "fuel_saved", "time_saved"]);

fprintf("\n--- 2-Formations Optimized ---\n")
toc %%%

%% GEOMETRIC TRAJECTORY OPTIMIZATION SOLO FLIGHTS
tic %%%

data_solo = zeros(nRows, 2);
plane_solo = string(zeros(nRows, 1));
crew_cost_solo = zeros(nRows,1);

for k = 1:nRows

    % Flight initial and final cruise locations
    f_dep_lat = gC_path_lat(k,tDep_min(k));
    f_dep_lon = gC_path_lon(k,tDep_min(k));
    f_arr_lat = gC_path_lat(k,tArr_min(k));    
    f_arr_lon = gC_path_lon(k,tArr_min(k));

    % Import performace data
    f_planeType = rawData.PlaneType{k};
    AC_perf = perf_data.(f_planeType);

    % Calculate great circle distance for solo flight without formation
    R_solo = distance(f_dep_lat, f_dep_lon, f_arr_lat, f_arr_lon, wgs84);

    % mission optimization
    [m_fuel, h, t] = solo_mission_optimization(R_solo, AC_perf);

    % Total fuel burn and time 
    data_solo(k,:) = [m_fuel, t/60]; % Fuel in Kg, time in min
    plane_solo(k) = string(f_planeType);

    switch plane_solo(k)
        case "L"
        crew_cost_solo(k) = crew_cost_L;
        case "AM"
        crew_cost_solo(k) = crew_cost_AM;
        case "AH"
        crew_cost_solo(k) = crew_cost_AH;
        case "AJ"
        crew_cost_solo(k) = crew_cost_AJ;
        case "BM"
        crew_cost_solo(k) = crew_cost_BM;
        case "BH"
        crew_cost_solo(k) = crew_cost_BH;
        case "BJ"
        crew_cost_solo(k) = crew_cost_BJ;
    end
end

total_cost_solo = fuel_cost*data_solo(:,1) + crew_cost_solo.*data_solo(:,2)/60;

fprintf("\n--- Solo Flights Optimized ---\n")

toc %%%
%% MILP SCHEDULE SOLVER TOC CRITERIA
tic %%%

cost_formation = total_cost_formation(valid2Formations_TOC);
cost_formation_solo = total_cost_formation_solo(valid2Formations_TOC);
cost_solo = total_cost_solo;
valid_two_formations = two_formations(valid2Formations_TOC, :);

i = nRows; % Number of flights
j = nRows + length(cost_formation); % Number of possible formations (including solo)

f = [cost_solo; cost_formation]; % TOC per formation (first solo then formations)
intcon = 1:j; % All variables are integers
b = []; % No inequalities
A = []; % No inequalities

Aeq = zeros(i,j);
Aeq(:, 1:nRows) = eye(nRows);
for k = 1:length(valid_two_formations)
    Aeq(valid_two_formations(k, 1), nRows+k) = 1;
    Aeq(valid_two_formations(k, 2), nRows+k) = 1;
end

beq = ones(i, 1);
lb = zeros(j, 1); % Lover bound = 0 for all variables
ub = ones(j, 1); % Upper bound = 1 for all variables
x0 = []; % Default initial guess

options = optimoptions('intlinprog');
[formations_ID_TOC, formation_TOC] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, x0, options);

idx = logical(formations_ID_TOC(nRows+1:end)); % ID of 2-formations selected
TOC_saved_2Formations = (sum(cost_formation_solo(idx)) - sum(cost_formation(idx)))/sum(cost_formation_solo(idx))*100;
TOC_saved_total = (sum(cost_solo) - formation_TOC)/sum(cost_solo)*100;

fprintf("\n--- Total Operating Costs Criteria ---\n")
fprintf("There is a total of %d valid formations out of %d [Filter Efficiency = %.2f%%] \n", sum(valid2Formations_TOC), numel(valid2Formations_TOC), sum(valid2Formations_TOC)/numel(valid2Formations_TOC)*100)
fprintf("After optimizing the schedule, a total of %d formations have been use out of the %d valid formations\n", sum(idx), sum(valid2Formations_TOC))
fprintf("The TOC saving by the 2-formations is a %.2f%% of the TOC in solo missions\n", TOC_saved_2Formations)
fprintf("The Total Network TOC saving by the 2-formations is a %.2f%% of the TOC in solo missions\n", TOC_saved_total)

% Aircraft Type Plots for result analysis
figure()
aux = plane2Formations(valid2Formations_TOC,:);
planes_select = aux(idx,:);
T = array2table(planes_select);
T = renamevars(T, ["planes_select1", "planes_select2"], ["Aircraft_1", "Aircraft_2"]);
heatmap(T, "Aircraft_1", "Aircraft_2")
title("Heatmap of Aircraft Categories in 2-Formations with TOC Criteria")
xlabel("Aircraft 1")
ylabel("Aircraft 2")

%geoplots
aux = two_formations(valid2Formations_TOC,:);
flights_select = aux(idx,:);

figure()
depLat_select = [depLat(flights_select(:,1)); depLat(flights_select(:,2))];
depLon_select = [depLon(flights_select(:,1)); depLon(flights_select(:,2))];
geoplot(depLat_select, depLon_select, "o")
hold on

arrLat_select = [desLat(flights_select(:,1)); desLat(flights_select(:,2))];
arrLon_select = [desLon(flights_select(:,1)); desLon(flights_select(:,2))];
geoplot(arrLat_select, arrLon_select, "x")

[lat, lon] = track2(depLat_select, depLon_select, arrLat_select, arrLon_select);

gx = geoaxes; % Use geographical coordinate axis
plot(gx, lat, lon) 
title("Flights selected for 2-formations with TOC Criteria")
colororder(color([flights_select(:,1);flights_select(:,2)]))

%  Band time Plots for result analysis
depTime_select = (tDep_min(flights_select(:,1)) + tDep_min(flights_select(:,2)))/(60*2);
arrTime_select = (tArr_min(flights_select(:,1)) + tArr_min(flights_select(:,2)))/(60*2);

flightTime = arrTime_select - depTime_select;
flightType = (flightTime < 2)*1 + ((flightTime <= 4).*(flightTime >= 2))*2 + (flightTime > 4)*3;

% Correct to mean meridian local time
depTime_select = depTime_select - (timezone(depLon(flights_select(:,1))) + timezone(depLon(flights_select(:,2))))/2;
arrTime_select = arrTime_select - (timezone(desLon(flights_select(:,1))) + timezone(desLon(flights_select(:,2))))/2;

bound = sign(desLon(flights_select(:,1)) - depLon(flights_select(:,1))); % 1 East -1 West
bound(bound == -1) = ones(nnz(bound == -1), 1)*2;

slot_count = zeros(3,6,2);
for i = 1:length(depTime_select)
    j = bound(i);
    k = flightType(i);

    if depTime_select(i) < 8
        slot_count(k,1,j) = slot_count(k,1,j) + 1;
    elseif depTime_select(i) < 12
        slot_count(k,2,j) = slot_count(k,2,j) + 1;
    elseif depTime_select(i) < 16
        slot_count(k,3,j) = slot_count(k,3,j) + 1;
    elseif depTime_select(i) < 20
        slot_count(k,4,j) = slot_count(k,4,j) + 1;
    elseif depTime_select(i) < 24
        slot_count(k,5,j) = slot_count(k,5,j) + 1;
    else
        slot_count(k,6,j) = slot_count(k,6,j) + 1;
    end
end

figure()
for g = 1:2
    subplot(2,2,g)
    b = bar(slot_count(:,:,g)');
    ylabel("Number of formations")
    xticklabels(["Early Morning", "Morning", "Midday", "Afternoon", "Evening", "Night"])
    for k = 1:3
        xtips1 = b(k).XEndPoints;
        ytips1 = b(k).YEndPoints;
        labels1 = string(b(k).YData);
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
    ylim([0 i/3])
    legend("Short Haul", "Medium Haul", "Long Haul")
end

slot_count = zeros(3,6,2);
for i = 1:length(arrTime_select)
    j = bound(i);
    k = flightType(i);

    if arrTime_select(i) < 8
        slot_count(k,1,j) = slot_count(k,1,j) + 1;
    elseif arrTime_select(i) < 12
        slot_count(k,2,j) = slot_count(k,2,j) + 1;
    elseif arrTime_select(i) < 16
        slot_count(k,3,j) = slot_count(k,3,j) + 1;
    elseif arrTime_select(i) < 20
        slot_count(k,4,j) = slot_count(k,4,j) + 1;
    elseif arrTime_select(i) < 24
        slot_count(k,5,j) = slot_count(k,5,j) + 1;
    else
        slot_count(k,6,j) = slot_count(k,6,j) + 1;
    end
end

for g = 1:2
    subplot(2,2,g+2)
    b = bar(slot_count(:,:,g)');
    ylabel("Number of formations")
    xticklabels(["Early Morning", "Morning", "Midday", "Afternoon", "Evening", "Night"])
    for k = 1:3
        xtips1 = b(k).XEndPoints;
        ytips1 = b(k).YEndPoints;
        labels1 = string(b(k).YData);
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
    ylim([0 i/3])
    legend("Short Haul", "Medium Haul", "Long Haul")
end

subplot(2,2,1)
title("Departure LMT time band of East Bound 2-formations")
subplot(2,2,2)
title("Departure LMT time band of West Bound 2-formations")
subplot(2,2,3)
title("Arrival LMT time band of East Bound 2-formations")
subplot(2,2,4)
title("Arrival LMT time band of West Bound 2-formations")

sgtitle("Aircraft Bands of 2-formations with TOC Criteria")

toc %%%
%% MILP SCHEDULE SOLVER CO2 CRITERIA
tic %%%

cost_formation = data2Formations(valid2Formations_CO2,1);
cost_formation_solo = data2Formations(valid2Formations_CO2,2); % Cost of solo missions of the formation
cost_solo = data_solo(:,1);
valid_two_formations = two_formations(valid2Formations_CO2, :);

i = nRows; % Number of flights
j = nRows + length(cost_formation); % Number of possible formations (including solo)

f = [cost_solo; cost_formation]; % TOC per formation (first solo then formations)
intcon = 1:j; % All variables are integers
b = []; % No inequalities
A = []; % No inequalities

Aeq = zeros(i,j);
Aeq(:, 1:nRows) = eye(nRows);
for k = 1:length(valid_two_formations)
    Aeq(valid_two_formations(k, 1), nRows+k) = 1;
    Aeq(valid_two_formations(k, 2), nRows+k) = 1;
end

beq = ones(i, 1);
lb = zeros(j, 1); % Lover bound = 0 for all variables
ub = ones(j, 1); % Upper bound = 1 for all variables
x0 = []; % Default initial guess

options = optimoptions('intlinprog');
[formations_ID_CO2, formation_CO2] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, x0, options);

idx = logical(formations_ID_CO2(nRows+1:end)); % ID of 2-formations selected
CO2_saved_2Formations = (sum(cost_formation_solo(idx)) - sum(cost_formation(idx)))/sum(cost_formation_solo(idx))*100;
CO2_saved_Total = (sum(cost_solo) - formation_CO2)/sum(cost_solo)*100;

fprintf("\n--- Total CO2 Emmision Criteria ---\n")
fprintf("There is a total of %d valid formations out of %d [Filter Efficiency = %.2f%%] \n", sum(valid2Formations_CO2), numel(valid2Formations_CO2), sum(valid2Formations_CO2)/numel(valid2Formations_CO2)*100)
fprintf("After optimizing the schedule, a total of %d formations have been use out of the %d valid formations\n", sum(idx), sum(valid2Formations_CO2))
fprintf("The CO2 emmision saving by the 2-formations is a %.2f%% of the emissions in solo missions\n", CO2_saved_2Formations)
fprintf("The Total Network CO2 emmision saving by the 2-formations is a %.2f%% of the emissions in solo missions\n", CO2_saved_Total)

%  Plane type Plots for result analysis
figure()
aux = plane2Formations(valid2Formations_CO2,:);
planes_select = aux(idx,:);
T = array2table(planes_select);
T = renamevars(T, ["planes_select1", "planes_select2"], ["Aircraft_1", "Aircraft_2"]);
heatmap(T, "Aircraft_1", "Aircraft_2")
title("Heatmap of Aircraft Categories in 2-Formations with CO2 Criteria")
xlabel("Aircraft 1")
ylabel("Aircraft 2")

% geoplots
aux = two_formations(valid2Formations_CO2,:);
flights_select = aux(idx,:);
figure()
depLat_select = [depLat(flights_select(:,1)); depLat(flights_select(:,2))];
depLon_select = [depLon(flights_select(:,1)); depLon(flights_select(:,2))];
% geoplot(depLat_select, depLon_select, "o")
% hold on

arrLat_select = [desLat(flights_select(:,1)); desLat(flights_select(:,2))];
arrLon_select = [desLon(flights_select(:,1)); desLon(flights_select(:,2))];
% geoplot(arrLat_select, arrLon_select, "x")

[lat, lon] = track2(depLat_select, depLon_select, arrLat_select, arrLon_select);

gx = geoaxes; % Use geographical coordinate axis
plot(gx, lat, lon) 
title("Flights selected for 2-formations with CO2 Criteria")
colororder(color([flights_select(:,1);flights_select(:,2)]))

%  Band time Plots for result analysis
depTime_select = (tDep_min(flights_select(:,1)) + tDep_min(flights_select(:,2)))/(60*2);
arrTime_select = (tArr_min(flights_select(:,1)) + tArr_min(flights_select(:,2)))/(60*2);

flightTime = arrTime_select - depTime_select;
flightType = (flightTime < 2)*1 + ((flightTime <= 4).*(flightTime >= 2))*2 + (flightTime > 4)*3;

% Correct to mean meridian local time
depTime_select = depTime_select - (timezone(depLon(flights_select(:,1))) + timezone(depLon(flights_select(:,2))))/2;
arrTime_select = arrTime_select - (timezone(desLon(flights_select(:,1))) + timezone(desLon(flights_select(:,2))))/2;

bound = sign(desLon(flights_select(:,1)) - depLon(flights_select(:,1))); % 1 East -1 West
bound(bound == -1) = ones(nnz(bound == -1), 1)*2;

slot_count = zeros(3,6,2);
for i = 1:length(depTime_select)
    j = bound(i);
    k = flightType(i);

    if depTime_select(i) < 8
        slot_count(k,1,j) = slot_count(k,1,j) + 1;
    elseif depTime_select(i) < 12
        slot_count(k,2,j) = slot_count(k,2,j) + 1;
    elseif depTime_select(i) < 16
        slot_count(k,3,j) = slot_count(k,3,j) + 1;
    elseif depTime_select(i) < 20
        slot_count(k,4,j) = slot_count(k,4,j) + 1;
    elseif depTime_select(i) < 24
        slot_count(k,5,j) = slot_count(k,5,j) + 1;
    else
        slot_count(k,6,j) = slot_count(k,6,j) + 1;
    end
end

figure()
for g = 1:2
    subplot(2,2,g)
    b = bar(slot_count(:,:,g)');
    ylabel("Number of formations")
    xticklabels(["Early Morning", "Morning", "Midday", "Afternoon", "Evening", "Night"])
    for k = 1:3
        xtips1 = b(k).XEndPoints;
        ytips1 = b(k).YEndPoints;
        labels1 = string(b(k).YData);
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
    ylim([0 i/3])
    legend("Short Haul", "Medium Haul", "Long Haul")
end

slot_count = zeros(3,6,2);
for i = 1:length(arrTime_select)
    j = bound(i);
    k = flightType(i);

    if arrTime_select(i) < 8
        slot_count(k,1,j) = slot_count(k,1,j) + 1;
    elseif arrTime_select(i) < 12
        slot_count(k,2,j) = slot_count(k,2,j) + 1;
    elseif arrTime_select(i) < 16
        slot_count(k,3,j) = slot_count(k,3,j) + 1;
    elseif arrTime_select(i) < 20
        slot_count(k,4,j) = slot_count(k,4,j) + 1;
    elseif arrTime_select(i) < 24
        slot_count(k,5,j) = slot_count(k,5,j) + 1;
    else
        slot_count(k,6,j) = slot_count(k,6,j) + 1;
    end
end

for g = 1:2
    subplot(2,2,g+2)
    b = bar(slot_count(:,:,g)');
    ylabel("Number of formations")
    xticklabels(["Early Morning", "Morning", "Midday", "Afternoon", "Evening", "Night"])
    for k = 1:3
        xtips1 = b(k).XEndPoints;
        ytips1 = b(k).YEndPoints;
        labels1 = string(b(k).YData);
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
    ylim([0 i/3])
    legend("Short Haul", "Medium Haul", "Long Haul")
end

subplot(2,2,1)
title("Departure LMT time band of East Bound 2-formations")
subplot(2,2,2)
title("Departure LMT time band of West Bound 2-formations")
subplot(2,2,3)
title("Arrival LMT time band of East Bound 2-formations")
subplot(2,2,4)
title("Arrival LMT time band of West Bound 2-formations")


sgtitle("Aircraft Bands of 2-formations with CO2 Criteria")


toc %%%