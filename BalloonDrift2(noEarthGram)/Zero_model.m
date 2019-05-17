%% Purdue Orbital
%
% HELIUM GAS LIFTING MODEL V2.0
%
% Written by Matthew Powell, Drew Sherman, Ethan Wahl, Michael Bailey, Dante Sanaei, and Josh Fitch
% Started 10/25/2016
% Updated 11/11/2017

%% Description of model
%
% Uses the basic buoyancy force equation to find force on a payload lifted
% by helium gas
%
% Time-steps through force, velocity, and position equations to find values
%
% Plots acceleration, velocity, and position vs time
% Plots volume vs altitude
%
% ---ASSUMPTIONS---
%   ideal gas
%   onion shaped balloon
%   standard atmosphere as used in AAE 251 at Purdue (see
%   	Standard_Atmosphere.m)
%   Temperature of helium is approximately the same as the ambient air at
%   initial conditiaon.
%   No loss of helium during ascent (constant mass)
%   Mass of entire system constant during ascent
%   	Mass of helium and balloon needs to be included in weight fighting
%   	the buoyancy force - recursive, may be possible to estimate
%   Helium is initially @ STP

%3,000 m^3 with a membrane thickness of 80 micrometers
%100 kg payload
% kg mass of balloon
% Test launch:
% .5216 kg balloon	130 cubic feet volume (3.86 )	5 lbs payload  .62 grams helium
% mass?
%% Helium Lifting Model V2.0
% 88023.05 g
%Zero_model(1, 3, 10, 5, 1)

function [altitude_array] = Zero_model(payload_mass, volume_balloon, mass_balloon, A_vent, mass_helium, plot_suppression)
% inputs :
%   payload_mass:   	wanted payload mass in kg
%   volume:         	Volume of the Balloon (m^3)
%   mass_balloon:   	Mass of balloon in kg
%   A_vent:         	Area of the zero pressure vent in m^2
%   mass_helium:    	Mass of helium in kg
%   plot_suppression:   0: no plots, anything else: plots


clc
close all

% set arrays for use later

delta_t = 0.25;                     	% Change in time for each iteration of the loop [s]
Cd = 0.3;                           	% default drag coefficient for spherical balloons
Tg = [];                            	% Gas temp inside balloon (assumed to be Tatm initially)
P_a = [];                           	% air pressure at altitude
    
R_He = 2077;                        	% Gas constant of helium J/Kg
Cp_he = 5.300*1000;                 	% Specific heat of helium J/Kg*C
rho_g = 0.1786;
g = 9.81;                           	% [m/s^2]
G = 6.67408 * 10^-11;               	% Gravitation Constant
mol_mass = .004002602;              	% Molar mass of heliumkg/mol
R = 2078.5;                         	% Specific Gas constant of Helium [J/kg-K]
cylinder_coef = 1.233;              	% kg.(He)/Cylinder(244 ft.^3)
C_nozzle = 0.5;

M_earth = 5.972 * 10^24;            	% mass of earth [kg]
R_earth = 6371000;                  	% radius of earth [m]

mol_mass_atmo = .029;               	% kg/mol
mol_ratio = mol_mass_atmo / mol_mass;   % Molar mass ratio [N/A]




%% Find intial conditions (namely He volume)
% diameter = 1.383 * volume^(1/3);	%Use diameter equation from Helium Ascent 3D paper
%
% Ab = pi * (diameter(end) / 2) ^ 2;              	%Cross-sectional Area of the balloon

Msys = payload_mass + mass_balloon;   %Total Mass of the whole system

Mtot = Msys + mass_helium;

num_moles = mass_helium / mol_mass;   %mol ratio (constant, super pressure balloon)

num_cylinders = mass_helium / cylinder_coef;

[~, P_air, rho_air, T_air] = Standard_Atmosphere(0, 0);
volume = mass_helium*R*T_air/P_air;
F_bouy = rho_air * g * volume - mass_helium * g;
F_bouy_array  = [F_bouy];                    	% buoyancy force [N]
accel_array = [0];                         	% acceleration seen by system [m/s^2]
velocity_array = [0];               	% velocity of the system [m/s]
altitude_array = [0];                   	% vertical position of system [m]
balloon_vol =[0];                    	% volume of helium [m^3]
F_tot_array  = [F_bouy-g*Mtot];
F_drag_array = [0];
t_array = [0];
delta_m_array = [0];
t = 0;

%% Loop Through Kinetics
% begin time_stepping
%while (velocity(end)>=0 && system_alt(end) < final_alt && system_alt(end) >= 0 )
i = 1;
velocity = 0.0000000000000001;
altitude = 0;
max_alt = 0;

 while (t_array(end) < 50000 )%&& altitude < 10000) %&& velocity > 0)   
	[~, P_air, rho_air, T_air] = Standard_Atmosphere(altitude, 0);   %Standard Atmosphere function from AAE 251
	P_a = [P_a, P_air];
    
	%Update g
	if i ~= 1
    	g = G * (M_earth) / (altitude + R_earth)^2;
	end
    
	volume = mass_helium*R*T_air/P_air;
	if volume > volume_balloon
    	volume = volume_balloon;
	end
	diameter = 2.23 * power(0.75 * volume / pi,1/3);
	proj_area = pi / 4 * diameter^2;
	delta_m=0;
    
	if max_alt == 0
    	F_bouy = rho_air * g * volume - mass_helium * g;
    	F_drag = 0.5 * Cd * rho_air * velocity^2 * proj_area * (velocity/abs(velocity));
    	F_grav = g * Msys;
    	F_tot  = F_bouy - F_drag - F_grav;
   	 
    	accel = F_tot / Mtot;
   	 
    	velocity = accel * delta_t + velocity_array(end);
   	 
    	altitude = velocity * delta_t + 0.5 * accel * delta_t^2 + altitude_array(end);
   	 
	else
    	F_bouy = 0;
    	F_drag = 0;
    	F_grav = 0;
    	F_tot  = 0;
    	accel = 0;
    	velocity = 0;
    	altitude = 25000;
	end
	if altitude >= 25000
    	max_alt = 1;
	end
    
	if volume >= volume_balloon && accel >= 0
    	rho_he = mass_helium / volume_balloon;
    	b = g * (rho_air - mass_helium/volume_balloon);
    	delta_p = b * 0.517 * power(volume_balloon,1/3);
    	delta_m = -A_vent * C_nozzle * sqrt(2 * delta_p * rho_he);
    
    	mass_helium = mass_helium + delta_m;
    	num_moles = mass_helium / mol_mass;   
    	fprintf('delta m = %d   volume = %d 	velocity = %d\n', delta_m,volume,velocity)
	end
    
	t = t + delta_t;
    
	F_tot_array = [F_tot_array, F_tot];
	F_bouy_array = [F_bouy_array, F_bouy];
	F_drag_array = [F_drag_array, F_drag];
	accel_array = [accel_array, accel];
	velocity_array = [velocity_array, velocity];
	altitude_array = [altitude_array, altitude];
	delta_m_array = [delta_m_array,delta_m];
	balloon_vol = [balloon_vol, volume];
	t_array = [t_array, t];
    
% 	disp(F_tot);
% 	disp(F_bouy);
% 	disp(F_drag);
% 	disp(F_grav);
% 	disp(accel);
% 	disp(velocity);
% 	disp(altitude);
% 	disp(volume);
% 	fprintf('NEW LINE\n');
    
	i = i + 1;
 end
    
 %% Plotting
if (plot_suppression ~= 0)
	hold on

	subplot(3,2,1)
	plot(t_array, altitude_array);
	xlabel('Time (s)');
	ylabel('Altitude (m)');
	title('Altitude vs Time for given Payload Mass')
   
	subplot(3,2,2)
	plot(t_array, accel_array);
	xlabel('Time (s)');
	ylabel('Acceleration (m/s^2)');
	title('Acceleration vs Time for given Payload Mass')

	subplot(3,2,3)
	plot(t_array, velocity_array);grid
	xlabel('Time (s)');
	ylabel('Velocity (m/s)');
	title('Velocity vs Time for given Payload Mass')
    
	subplot(3,2,4)
	plot(t_array, F_tot_array);
	xlabel('time(s)');
	ylabel('Total Force(F)');
	title('Total Force vs Time')
    
	subplot(3,2,5)
	plot(t_array, F_bouy_array)
	xlabel('time(s)');
	ylabel('Fb(F)');
	title('Bouyant Force vs Time')
    
	subplot(3,2,6)
	plot(t_array, F_drag_array)
	xlabel('Time (s)')
	ylabel('Drag (N)')
	title('Drag vs Time')
    
	fprintf('The total number of helium cansiters needed are: %f\n', num_cylinders);
	fprintf('The total number of moles of helium needed are: %f\n', num_moles);
	fprintf('The final altitude is : %f\n', altitude_array(end));

end
