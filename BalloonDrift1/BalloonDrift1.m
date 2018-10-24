%% Purdue Orbital 
% 
% BALLOON DRIFT MODEL V1.0
%
% Written by Dante Sanaei
% Written: September 5, 2018
% Update 1: 10/15/18

%% Description of model
% 
% Uses an inputted array of altitude vs time and sounding data from the
% University of Wyoming (http://weather.uwyo.edu/upperair/sounding.html) 
%
% Time-steps through force, velocity, and position to find new
% coordninates at each altitude represented on the sounding data.
%
% For every modeled altitude within each sounding height, 1 second is
% counted and the total amount of seconds is the delta T in the time step
% 
% 
% ---ASSUMPTIONS---
%   standard atmosphere as used in AAE 251 at Purdue (see
%       Standard_Atmosphere.m)
%   Wind speed and direction is the same as the sounding data at the ILN
%   sounding station
%
%'aug31.rtf'

%% Function
function [klmData] = BalloonDrift1(alt,soundingFile, area, mass,startLat,startLon, timeInterval)
% inputs :
% alt:                array of altitudes (m)
% soundingFile:       .txt file of wyoming sounding data. (remeber to
%                     delete any columns or rows that have empty spots)
% area:               Surface area of balloon (m^2) 2
% mass:               Total mass of balloon(g) 3500
% startLat:           Latitude of launch site (purdue: 40.41279385)
% startLon:           Longitude of launch site (Purdue: -86.93937395889625)
% timeInterval:       the interval of time between each altitude
%                     measurement(Matt Powell: 1, Astra: 3)

%% Initializations
clc
data = load(soundingFile);                  %Input sounding data
height = data(:,2);                         %Initialize sounding heights (m)
windSpeed = data(:,8) * .5144444;           %Initialize sounding windspeeds (m/s)
windDirection = data(:,7);                  %Initialize sounding wind directions (degrees true)
[alt, t] = Zero_model(1, 5, 1.5, .5, 0);    %Zero_model(payload_mass, volume, mass_balloon, mass_helium, plot_suppression);
%alt = xlsread('astratest1.xls');
location = [0 0];                           %Locations (distances in m)
locs = []; 
latlon = [startLat, startLon]; %Starting coordniates (Purdue Airport)
latlon1 = latlon;                           
initV = [0 0];
heights = [0];

%% Time Step
for i = 2:length(height)                                                % loops through each height in the sounding data
    time = sum(alt < height(i) & alt > height(i-1));                    % the amount of altitude points in each interval is equal to time (s)
    time = time * timeInterval;
    Xvelocity = -windSpeed(i)*cosd(windDirection(i));                   % finding x component of wind velocity (m/s) 
    yvelocity = -windSpeed(i)*sind(windDirection(i));                   % finding y component of wind velocity (m/s) 
    windV = [Xvelocity, yvelocity];                                     % wind velocity vector (m/s)
    [GeopotentialAlt, DesiredPressure, rho, DesiredTemperature] = Standard_Atmosphere(height(i), 0);    %finds pressure using 1976 Standard Atmosphere
    vel = ((((.5)*(rho ).*([Xvelocity, yvelocity].^2).*area)/mass).* time + initV); %finds balloon velocity vector using wind load. 
    initV = windV;                                                      
    distance = vel .* time;                                             %finds distance (m)
    %distance = [Xvelocity, yvelocity] .* time;                         %if the previous method does not work comment it out and use this line only
    location = location + distance;
    locs = [locs;location];
    lat = latlon(1) + (distance(1) * 0.000621371) / 68.703;             %convert x distance (m) to latitude
    long = latlon(2) + (distance(2) * 0.000621371) /  55.2428;          %convert y distance (m) to longitude
    latlon = [lat long ];
    latlon1 = [latlon1; latlon];
    heights = [heights; height(i)];
end


%% Exporting Data to .CSV
klmData = [heights latlon1];
%dlmwrite('bestone.csv', xml, 'delimiter', ',', 'precision', 12,'roffset',1);
cHeader = {'Altitude' 'Latitude' 'Longitude'}; %dummy header
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
%write header to file
fid = fopen('KLMFILE.csv','w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
%write data to end of file
dlmwrite('KLMFILE.csv', klmData, 'delimiter', ',', 'precision', 12,'-append');


