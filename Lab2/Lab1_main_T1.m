clear all; close all; clc;

%% Load Data
load('run_7');
load('run_9');
load('run_11');

[T, a, P, rho] = atmoscoesa(1624);

% Parse Data from Run with 7 degree temperature difference
time7 = run_7(:,1);
pressures7 = run_7(:,2)*6894.76 + P;
topOfTop7 = run_7(:,3)+273.15;
bottomOfTop7 = run_7(:,4)+273.15;
topOfBottom7 = run_7(:,5)+273.15;
bottomOfBottom7 = run_7(:,6)+273.15;
current7 = run_7(:,7);
optSwitch7 = run_7(:,8);

% Parse Data from Run with 9 degree temperature difference
time9 = run_9(:,1);
pressures9 = run_9(:,2)*6894.76 + P;
topOfTop9 = run_9(:,3)+273.15;
bottomOfTop9 = run_9(:,4)+273.15;
topOfBottom9 = run_9(:,5)+273.15;
bottomOfBottom9 = run_9(:,6)+273.15;
current9 = run_9(:,7);
optSwitch9 = run_9(:,8);

% Parse Data from Run with 11 degree temperature difference
time11 = run_11(:,1);
pressures11 = run_11(:,2)*6894.76 + P;
topOfTop11 = run_11(:,3)+273.15;
bottomOfTop11 = run_11(:,4)+273.15;
topOfBottom11 = run_11(:,5)+273.15;
bottomOfBottom11 = run_11(:,6)+273.15;
current11 = run_11(:,7);
optSwitch11 = run_11(:,8);

%% Find Extreme Temperatures

tempHigh7 = max((topOfBottom7 + bottomOfBottom7)./2);
tempLow7  = min((topOfTop7 + bottomOfTop7)./2);

tempHigh9 = max((topOfBottom9 + bottomOfBottom9)./2);
tempLow9  = min((topOfTop9 + bottomOfTop9)./2);

tempHigh11 = max((topOfBottom11 + bottomOfBottom11)./2);
tempLow11  = min((topOfTop11 + bottomOfTop11)./2);

%% Efficiency based on reversible cycle

thermalEfficiency7  = (1 - (tempLow7 / tempHigh7))*100;
thermalEfficiency9  = (1 - (tempLow9 / tempHigh9))*100;
thermalEfficiency11 = (1 - (tempLow11 / tempHigh11))*100;

%% From Motion Study

piston_disp = xlsread('MoStudy7C.xlsx');

r_power = 15.5/2000; %power piston raduis (diameter in mm /2 /1000 = radius in m)
h = 10/1000; % power piston max height (mm/1000 = m) unused

piston_disp = piston_disp(:,2:3); % getting absolute power piston displacements
heights = -(min(piston_disp(:,2)) - piston_disp(:,2))/1000;

vol_power = @(height) pi*r_power^2*height; %vector of power piston volumes
vol_main = 1.7267*10^-4; %main chamber volume [m^3]
vol = vol_power(heights); 

vol_total = vol + vol_main;

total_cycle_piston = 177;
vol_tot = vol_total(1:total_cycle_piston);

%% Run 7

%indices for a full cycle from run_7 
start_index7 = 1511;
end_index7 = 3465;

% pressures on the index of that cycle
pressures7 = pressures7(start_index7:end_index7);

%found interparc online to interpolate cyclical functions
vol_tot7 = interparc(length(pressures7),piston_disp(1:total_cycle_piston,1),vol_tot,'spline');

figure
plot(vol_tot7(1:end-5,2),pressures7(1:end-5)) %pressure psi -> Pa
xlabel('Volume (m^3)')
ylabel('Pressure (Pa)')
title('P-V Diagram for 7 Degree Temperature Difference')

% finding min and max volumes
[minVol,idxMin]= min(vol_tot7(:,2));
[maxVol,idxMax] = max(vol_tot7(:,2));

% getting top half of curve for work out calculation
topHalf7(1:idxMin-1,1) = vol_tot7(1:idxMin-1,2);
topHalf7(idxMax+1:length(vol_tot7),1) = vol_tot7(idxMax+1:end,2);
topHalf7(1:idxMin-1,2) = pressures7(1:idxMin-1);
topHalf7(idxMax+1:length(vol_tot7),2) = pressures7(idxMax+1:end);
topHalf7 = flip(topHalf7);
topHalf7 = [topHalf7(1466:end,:);topHalf7(1:485,:)];

% Work under top half of curve
workOut7 = polyarea([topHalf7(1,1);topHalf7(:,1);topHalf7(end,1)], [0;topHalf7(:,2);0]);

% Work net and work in
workNet7 = polyarea(vol_tot7(1:end-5,2),pressures7(1:end-5));
wIn7 = workOut7-workNet7;

% Heat transfer
qIn7 = trapz(run_7(start_index7:end_index7,1),48.*run_7(start_index7:end_index7,7)); % Joules

% Thermal Efficiency
efficiency7 = (workNet7/qIn7)*100;
fprintf('Net work for 7 degree temp difference is %.6f \nEfficiency for 7 degree temp difference is %.4f%%\n',workNet7,efficiency7)


%% Run 9

% Indexing the cycle
start_index9 = 911;
end_index9 = 2143;

% pressures on the index of that cycle
pressures9 = pressures9(start_index9:end_index9);

%found interparc online to interpolate cyclical functions
vol_tot9 = interparc(length(pressures9),piston_disp(1:total_cycle_piston,1),vol_tot,'spline');

figure
plot(vol_tot9(1:end-5,2),pressures9(1:end-5)) %pressure psi -> Pa
xlabel('Volume (m^3)')
ylabel('Pressure (Pa)')
title('P-V Diagram for 9 Degree Temperature Difference')

% finding min and max volumes
[minVol,idxMin]= min(vol_tot9(:,2));
[maxVol,idxMax] = max(vol_tot9(:,2));

% getting top half of curve for work out calculation
topHalf9(1:idxMin-1,1) = vol_tot9(1:idxMin-1,2);
topHalf9(idxMax+1:length(vol_tot9),1) = vol_tot9(idxMax+1:end,2);
topHalf9(1:idxMin-1,2) = pressures9(1:idxMin-1);
topHalf9(idxMax+1:length(vol_tot9),2) = pressures9(idxMax+1:end);
topHalf9 = flip(topHalf9);
topHalf9 = [topHalf9(1466:end,:);topHalf9(1:485,:)];

% Work under top half of curve
workOut9 = polyarea([topHalf9(1,1);topHalf9(:,1);topHalf9(end,1)], [0;topHalf9(:,2);0]);

% Work net and work in
workNet9 = polyarea(vol_tot9(1:end-5,2),pressures9(1:end-5));
wIn9 = workOut9-workNet9;

% Heat transfer
qIn9 = trapz(run_9(start_index9:end_index9,1),48.*run_9(start_index9:end_index9,7)); % Joules

% Thermal Efficiency
efficiency9 = (workNet9/qIn9)*100;
fprintf('Net work for 9 degree temp difference is %.6f \nEfficiency for 9 degree temp difference is %.4f%%\n',workNet9,efficiency9)

%% Run 11
% Finding the index of one cycle
start_index11 = 507;
end_index11 = 1452;

% pressures on the index of that cycle
pressures11 = pressures11(start_index11:end_index11);

%found interparc online to interpolate cyclical functions
vol_tot11 = interparc(length(pressures11),piston_disp(1:total_cycle_piston,1),vol_tot,'spline');

figure
plot(vol_tot11(1:end-5,2),pressures11(1:end-5)) %pressure psi -> Pa
xlabel('Volume (m^3)')
ylabel('Pressure (Pa)')
title('P-V Diagram for 11 Degree Temperature Difference')

% finding min and max volumes
[minVol,idxMin]= min(vol_tot11(:,2));
[maxVol,idxMax] = max(vol_tot11(:,2));

% getting top half of curve for work out calculation
topHalf11(1:idxMin-1,1) = vol_tot11(1:idxMin-1,2);
topHalf11(idxMax+1:length(vol_tot11),1) = vol_tot11(idxMax+1:end,2);
topHalf11(1:idxMin-1,2) = pressures11(1:idxMin-1);
topHalf11(idxMax+1:length(vol_tot11),2) = pressures11(idxMax+1:end);
topHalf11 = flip(topHalf11);
topHalf11 = [topHalf11(1466:end,:);topHalf11(1:485,:)];

% Work under top half of curve
workOut11 = polyarea([topHalf11(1,1);topHalf11(:,1);topHalf11(end,1)], [0;topHalf11(:,2);0]);

% Work net and work in
workNet11 = polyarea(vol_tot11(1:end-5,2),pressures11(1:end-5));
wIn11 = workOut11-workNet11;

% Heat transfer
qIn11 = trapz(run_11(start_index11:end_index11,1),48.*run_11(start_index11:end_index11,7)); % Joules

% Thermal Efficiency
efficiency11 = (workNet11/qIn11)*100;
fprintf('Net work for 11 degree temp difference is %.6f \nEfficiency for 11 degree temp difference is %.4f%%\n',workNet11,efficiency11)
