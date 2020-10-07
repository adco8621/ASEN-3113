clear;clc;close all;

%% data load

load('run_7')
load('run_9')
load('run_11')

[T, a, P, rho] = atmoscoesa(1624);

%% rpm calc

% 7 degree diff
figure
plot(run_7(:,1),run_7(:,2))               % plotted the pressure vs time and pulled
time7 = [0 1.165 2.359 3.532 4.701 5.88]; % the times out manually and hardcoded
diffs7 = diff(time7); % finding difference in times
rpm7 = mean((diffs7./60).^-1); % converting [sec/rev] to rpm

% 9 degree diff
figure
plot(run_9(:,1),run_9(:,2));
time9 = [0 0.7484 1.488 2.233 2.985 3.741 4.475 5.223 5.977];
diffs9 = diff(time9);
rpm9 = mean((diffs9./60).^-1);

% 11 degree diff
figure
plot(run_11(:,1),run_11(:,2));
time11 = [0.4698 1.042 1.61 2.185 2.755 3.326];
diffs11 = diff(time11);
rpm11 = mean((diffs11./60).^-1);

% found the average temperatures of the entire plates and then subtracted
% the top avg from bottom
% avg = -((top of top + bot of top)/2) + ((bot of top + top of top)/2)
avg7 = -(mean(run_7(:,3))+mean(run_7(:,4)))/2 + (mean(run_7(:,5))+mean(run_7(:,6)))/2;
avg9 = -(mean(run_9(:,3))+mean(run_9(:,4)))/2 + (mean(run_9(:,5))+mean(run_9(:,6)))/2;
avg11 = -(mean(run_11(:,3))+mean(run_11(:,4)))/2 + (mean(run_11(:,5))+mean(run_11(:,6)))/2;

% temps for efficeincey ca;cs
temps(1) = mean(run_7(:,3))+mean(run_7(:,4))/2;
temps(2) = mean(run_7(:,5))+mean(run_7(:,6))/2;
temps(3) = mean(run_9(:,3))+mean(run_9(:,4))/2;
temps(4) = mean(run_9(:,5))+mean(run_9(:,6))/2;
temps(5) = mean(run_11(:,3))+mean(run_11(:,4))/2;
temps(6) = mean(run_11(:,5))+mean(run_11(:,6))/2;

piston_disp = xlsread('MoStudy7C.xlsx');

r_power = 7.75/2000; %power piston raduis (diameter in mm /2 /1000 = radius in m)
h = 10/1000; % power piston max height (mm/1000 = m) unused

piston_disp = piston_disp(:,2:3); % getting absolute power piston displacements
heights = -(min(piston_disp(:,2)) - piston_disp(:,2))/1000;

vol_power = @(height) pi*r_power^2*height; %vector of power piston volumes
vol_main = 1.7267*10^-4; %main chamber volume [m^3]
vol = vol_power(heights);

vol_tot = vol + vol_main;


%% Run 7

%the '-225' below was just me messing around trying to fix the graphs and
%it worked. The original numbers should be the proper place to start and
%end the cycle but it isnt for some reason
start_index = 1511;
end_index = 3465; %indices for a full cycle from data
total_cycle = end_index-start_index; %size of a cycle

% used to find indices for motion study cycle
%next = find((piston_disp(2:end,2) > (piston_disp(1,2) - 0.05)) & (piston_disp(2:end,2) < (piston_disp(1,2) + 0.05)));

total_cycle_piston = 178;
pressures = run_7(start_index:end_index,2);
vol_tot = vol_tot(1:total_cycle_piston);

% interp doesnt work on non-functions (i.e. doesnt work on cyclical plots)
%vol_tot = interp1(piston_disp(1:total_cycle_piston,1),vol_tot,linspace(1,piston_disp(total_cycle_piston,1),1955));

%found interparc online to interpolate cyclical functions
vol_tot7 = interparc(1955,piston_disp(1:total_cycle_piston,1),vol_tot,'spline');


Vf7 = vol_tot7(1:end-10,2); %=======removing the overlap===========
Vf7 = [Vf7; Vf7(1)];
Pf7 = pressures(1:end-10)*6894.76 + P;
Pf7 = [Pf7; Pf7(1)];
figure
plot(Vf7,Pf7) %pressure psi -> Pa

Imin = find(Vf7 == min(Vf7));
Imax = find(Vf7 == max(Vf7));

Int1 = flip(1:Imin);
Int2 = Imin:Imax;
Int3 = flip(Imax:length(Vf7));
%Wnet7 = abs(trapz(Vf7,Pf7))
Win7 = trapz(Vf7(Int2),Pf7(Int2));
Wout7 = trapz(Vf7(Int1),Pf7(Int1))+trapz(Vf7(Int3),Pf7(Int3));
Wnet7 = Wout7-Win7

%Wc = cumtrat
  
q_in7 = trapz(run_7(:,1),48.*run_7(:,7)); % Joules

efficiency7 = (Wnet7/q_in7)*100
%% Run 9

%the '-225' below was just me messing around trying to fix the graphs and
%it worked. The original numbers should be the proper place to start and
%end the cycle but it isnt for some reason
start_index = 909;
end_index = 2142; %indices for a full cycle from data
total_cycle = end_index-start_index; %size of a cycle

% used to find indices for motion study cycle
%next = find((piston_disp(2:end,2) > (piston_disp(1,2) - 0.05)) & (piston_disp(2:end,2) < (piston_disp(1,2) + 0.05)));

total_cycle_piston = 178;
pressures = run_9(start_index:end_index,2);
vol_tot = vol_tot(1:total_cycle_piston);

% interp doesnt work on non-functions (i.e. doesnt work on cyclical plots)
%vol_tot = interp1(piston_disp(1:total_cycle_piston,1),vol_tot,linspace(1,piston_disp(total_cycle_piston,1),1955));

%found interparc online to interpolate cyclical functions
vol_tot9 = interparc(total_cycle+1,piston_disp(1:total_cycle_piston,1),vol_tot,'spline');

Vf9 = vol_tot9(1:end-10,2); %=======removing the overlap===========
Vf9 = [Vf9; Vf9(1)];
Pf9 = pressures(1:end-10)*6894.76 + P;
Pf9 = [Pf9; Pf9(1)];
figure
plot(Vf9,Pf9) %pressure psi -> Pa

Imin = find(Vf9 == min(Vf9));
Imax = find(Vf9 == max(Vf9));
Int1 = flip(1:Imin);
Int2 = Imin:Imax;
Int3 = flip(Imax:length(Vf9));
%Wnet7 = abs(trapz(Vf7,Pf7))
Wout9 = trapz(Vf9(Int2),Pf9(Int2));
Win9 = trapz(Vf9(Int1),Pf9(Int1))+trapz(Vf9(Int3),Pf9(Int3));
Wnet9 = Wout9-Win9

q_in9 = trapz(run_9(:,1),48.*run_9(:,7)); % Joules

efficiency9 = (Wnet9/q_in9)*100
%% Run 11

%the '-225' below was just me messing around trying to fix the graphs and
%it worked. The original numbers should be the proper place to start and
%end the cycle but it isnt for some reason
start_index = 506;
end_index = 1451; %indices for a full cycle from data
total_cycle = end_index-start_index; %size of a cycle

% used to find indices for motion study cycle
%next = find((piston_disp(2:end,2) > (piston_disp(1,2) - 0.05)) & (piston_disp(2:end,2) < (piston_disp(1,2) + 0.05)));

total_cycle_piston = 178;
pressures = run_11(start_index:end_index,2);
vol_tot = vol_tot(1:total_cycle_piston);

% interp doesnt work on non-functions (i.e. doesnt work on cyclical plots)
%vol_tot = interp1(piston_disp(1:total_cycle_piston,1),vol_tot,linspace(1,piston_disp(total_cycle_piston,1),1955));

%found interparc online to interpolate cyclical functions
vol_tot11 = interparc(total_cycle+1,piston_disp(1:total_cycle_piston,1),vol_tot,'spline');


Vf11 = vol_tot11(1:end-10,2); %=======removing the overlap===========
Vf11 = [Vf11; Vf11(1)];
Pf11 = pressures(1:end-10)*6894.76 + P;
Pf11 = [Pf11; Pf11(1)];
figure
plot(Vf11,Pf11) %pressure psi -> Pa

Imin = find(Vf11 == min(Vf11));
Imax = find(Vf11 == max(Vf11));
Int1 = flip(1:Imin);
Int2 = Imin:Imax;
Int3 = flip(Imax:length(Vf11));
%Wnet7 = abs(trapz(Vf7,Pf7))
Wout11 = trapz(Vf11(Int2),Pf11(Int2));
Win11 = trapz(Vf11(Int1),Pf11(Int1))+trapz(Vf11(Int3),Pf11(Int3));
Wnet11 = Wout11-Win11

q_in11 = trapz(run_11(:,1),48.*run_11(:,7)); % Joules

efficiency11 = (Wnet11/q_in11)*100