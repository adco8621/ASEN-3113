clear;clc;close all;

solar_e = 1361;
rw = 147.1*10^6;
re = 149.598*10^6;
rs = 152.1*10^6;
solar_w = solar_e*(re/rw)^2;
solar_s = solar_e*(re/rs)^2;
inc = deg2rad(23.5);
area_s = 0.0838;
alpha = 0.1;
emiss = 0.8;
rea = 6371;
ro = 42164;
op = 20;
back_s = 63;
back_w = 88;
back_eq = 75.5;
back_ec = 11;
sigma = 5.670374419*10^-8;

eclipse = 2*asin(rea/(2*ro));

eclipseSection = [pi-eclipse/2 pi+eclipse/2];

totalOrbit = linspace(0,2*pi,1000);
time = totalOrbit*24/(2*pi);

% solar * sin(angle) * area * alpha + op + heater + backload*alpha*area = emiss *
% area_s*sigma*T^4

flux_s = zeros(1,1000);
flux_s(1:500) = solar_s*cos(inc)*sin(totalOrbit(1:500))*area_s*alpha;

flux_e = zeros(1,1000);
flux_e(1:488) = solar_e*sin(totalOrbit(1:488))*area_s*alpha;

flux_w = zeros(1,1000);
flux_w(1:500) = solar_w*cos(inc)*sin(totalOrbit(1:500))*area_s*alpha;

figure
hold on
grid on
plot(time,flux_s,time,flux_e,time,flux_w)
legend('summer','equinox','winter')
xlabel('Time [hr]')
ylabel('Solar Radiation Flux [W]')
xlim([0 24])

% solar * sin(angle) * area * alpha + op + heater + backload*alpha*area = emiss *
% area_s*sigma*T^4

T_unheated_s = zeros(1,1000);
T_unheated_s(1:500) = nthroot((solar_s*cos(inc)*sin(totalOrbit(1:500))*area_s*alpha...
    +op+back_s*area_s*alpha)/(emiss*area_s*sigma),4);
T_unheated_s(501:1000) = nthroot((op+back_s*area_s*alpha)/(emiss*area_s*sigma),4);

T_unheated_w = zeros(1,1000);
T_unheated_w(1:500) = nthroot((solar_w*cos(inc)*sin(totalOrbit(1:500))*area_s*alpha...
    +op+back_w*area_s*alpha)/(emiss*area_s*sigma),4);
T_unheated_w(501:1000) = nthroot((op+back_w*area_s*alpha)/(emiss*area_s*sigma),4);

T_unheated_e = zeros(1,1000);
T_unheated_e(1:488) = nthroot((solar_e*sin(totalOrbit(1:488))*area_s*alpha...
    +op+back_eq*area_s*alpha)/(emiss*area_s*sigma),4);
T_unheated_e(489:512) = nthroot((op+back_ec*area_s*alpha)/(emiss*area_s*sigma),4);
T_unheated_e(513:1000) = nthroot((op+back_eq*area_s*alpha)/(emiss*area_s*sigma),4);

T_unheated_s = T_unheated_s - 273.15;
T_unheated_e = T_unheated_e - 273.15;
T_unheated_w = T_unheated_w - 273.15;


% solar * sin(angle) * area * alpha + op + heater + backload*alpha*area = emiss *
% area_s*sigma*T^4

need_heat_w_on = find(T_unheated_w < 20);
need_heat_w_sun_on = need_heat_w_on(find(need_heat_w_on <= 500));
need_heat_w_no_sun_on = need_heat_w_on(find(need_heat_w_on > 500));
heat_w_on = zeros(1,1000);
heat_w_on(need_heat_w_sun_on) = emiss*area_s*sigma*(20+273.15)^4 - op - solar_w*cos(inc)*...
    sin(totalOrbit(need_heat_w_sun_on))*area_s*alpha - back_w*alpha*area_s;
heat_w_on(need_heat_w_no_sun_on) = emiss*area_s*sigma*(20+273.15)^4 - op - back_w*alpha*area_s;

need_heat_s_on = find(T_unheated_s < 20);
need_heat_s_sun_on = need_heat_s_on(find(need_heat_s_on <= 500));
need_heat_s_no_sun_on = need_heat_s_on(find(need_heat_s_on > 500));
heat_s_on = zeros(1,1000);
heat_s_on(need_heat_s_sun_on) = emiss*area_s*sigma*(20+273.15)^4 - op - solar_s*cos(inc)*...
    sin(totalOrbit(need_heat_s_sun_on))*area_s*alpha - back_s*alpha*area_s;
heat_s_on(need_heat_s_no_sun_on) = emiss*area_s*sigma*(20+273.15)^4 - op - back_s*alpha*area_s;

need_heat_e_on = find(T_unheated_e < 20);
need_heat_e_sun_on = need_heat_e_on(find(need_heat_e_on <= 488));
need_heat_e_ec_on = need_heat_e_on(find(need_heat_e_on > 488 & need_heat_e_on < 513));
need_heat_e_no_sun_on = need_heat_e_on(find(need_heat_e_on >= 513));
heat_e_on = zeros(1,1000);
heat_e_on(need_heat_e_sun_on) = emiss*area_s*sigma*(20+273.15)^4 - op - solar_e*cos(inc)*...
    sin(totalOrbit(need_heat_e_sun_on))*area_s*alpha - back_eq*alpha*area_s;
heat_e_on(need_heat_e_ec_on) = emiss*area_s*sigma*(20+273.15)^4 - op - back_ec*alpha*area_s;
heat_e_on(need_heat_e_no_sun_on) = emiss*area_s*sigma*(20+273.15)^4 - op - back_eq*alpha*area_s;


need_heat_w_off = find(T_unheated_w < -40);
need_heat_w_sun_off = need_heat_w_off(find(need_heat_w_off <= 500));
need_heat_w_no_sun_off = need_heat_w_off(find(need_heat_w_off > 500));
heat_w_off = zeros(1,1000);
heat_w_off(need_heat_w_sun_off) = emiss*area_s*sigma*(20+273.15)^4 - op - solar_w*cos(inc)*...
    sin(totalOrbit(need_heat_w_sun_off))*area_s*alpha - back_w*alpha*area_s;
heat_w_off(need_heat_w_no_sun_off) = emiss*area_s*sigma*(20+273.15)^4 - op - back_w*alpha*area_s;

need_heat_s_off = find(T_unheated_s < -40);
need_heat_s_sun_off = need_heat_s_off(find(need_heat_s_off <= 500));
need_heat_s_no_sun_off = need_heat_s_off(find(need_heat_s_off > 500));
heat_s_off = zeros(1,1000);
heat_s_off(need_heat_s_sun_off) = emiss*area_s*sigma*(20+273.15)^4 - op - solar_s*cos(inc)*...
    sin(totalOrbit(need_heat_s_sun_off))*area_s*alpha - back_s*alpha*area_s;
heat_s_off(need_heat_s_no_sun_off) = emiss*area_s*sigma*(20+273.15)^4 - op - back_s*alpha*area_s;

need_heat_e_off = find(T_unheated_e < -40);
need_heat_e_sun_off = need_heat_e_off(find(need_heat_e_off <= 488));
need_heat_e_ec_off = need_heat_e_off(find(need_heat_e_off > 488 & need_heat_e_off < 513));
need_heat_e_no_sun_off = need_heat_e_off(find(need_heat_e_off >= 513));
heat_e_off = zeros(1,1000);
heat_e_off(need_heat_e_sun_off) = emiss*area_s*sigma*(20+273.15)^4 - op - solar_e*cos(inc)*...
    sin(totalOrbit(need_heat_e_sun_off))*area_s*alpha - back_eq*alpha*area_s;
heat_e_off(need_heat_e_ec_off) = emiss*area_s*sigma*(20+273.15)^4 - op - back_ec*alpha*area_s;
heat_e_off(need_heat_e_no_sun_off) = emiss*area_s*sigma*(20+273.15)^4 - op - back_eq*alpha*area_s;



figure
hold on
grid on
plot(time,T_unheated_s,time,heat_s_on,time,heat_s_off)
legend('Unheated Temperature [C]','Operational Heater Power Draw [W]','Survival Heater Power Draw [W]')
xlabel('Time [hr]')
title('Summer')
xlim([0 24])

figure
hold on
grid on
plot(time,T_unheated_w,time,heat_w_on,time,heat_w_off)
legend('Unheated Temperature [C]','Operational Heater Power Draw [W]','Survival Heater Power Draw [W]')
xlabel('Time [hr]')
title('Winter')
xlim([0 24])

figure
hold on
grid on
plot(time,T_unheated_e,time,heat_e_on,time,heat_e_off)
legend('Unheated Temperature [C]','Operational Heater Power Draw [W]','Survival Heater Power Draw [W]')
xlabel('Time [hr]')
title('Equinox')
xlim([0 24])