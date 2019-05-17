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
function [klmData] = BalloonDrift3(alt,fileName,area,mass,phi,theta, timeInterval)
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
%alt = linspace(1,25000,5400);
%phi = 40.8500; 
%theta = -119.1238;
%month = 10;
%area = 10;
%mass = 30000;
%timeInterval = 1;
latlon = [phi, theta]; %Starting coordniates (Purdue Airport)
latlon1 = latlon;    
heights = 0;
initV = [0 0];
totalD = 0;
%% Time Step
for i = 1:length(alt)
    time = timeInterval;
    height =  alt(i) / 1000;    
    [windSpeed, windDir] = NOAAData(fileName, height, phi, 180+(180-abs(theta)));
    Xvelocity = -windSpeed*cosd(windDir);                   % finding x component of wind velocity (m/s) 
    yvelocity = -windSpeed*sind(windDir);
    windV = [Xvelocity, yvelocity];                                    % wind velocity vector (m/s)
    [GeopotentialAlt, DesiredPressure, rho, DesiredTemperature] = Standard_Atmosphere(height, 0);    %finds pressure using 1976 Standard Atmosphere
    vel = ((((.5)*(rho ).*([Xvelocity, yvelocity].^2).*area)/mass).* time + initV); %finds balloon velocity vector using wind load. 
    initV = windV;                                                      
    distance = vel .* time;                                             %finds distance (m)
    %distance = [Xvelocity, yvelocity] .* time  ; 
   %if the previous method does not work comment it out and use this line only
    theta = theta + (distance(1) * 0.000621371) / 55.2428;             %convert x distance (m) to longitutde
    phi = phi + (distance(2) * 0.000621371) / 68.703 ;          %convert y distance (m) to latitude
    latlon = [phi theta];
    latlon1 = [latlon1; latlon];
    height = height * 1000;
    heights = [heights; height];
    totalD = totalD + distance;
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


