%% Purdue Orbital 
% 
% BALLOON DRIFT MODEL V2.0
%
% Written by Dante Sanaei and Michael Baily
% Written: October 22, 2018

%% Description of model
% This model uses average wind data from EARTHGRAM in order to model a high
% altitude balloon
%
% 
% 
% ---ASSUMPTIONS---
%   standard atmosphere as used in AAE 251 at Purdue (see
%       Standard_Atmosphere.m)
%   
%


%% Function
function [klmData] = BalloonDrift2(alt, area, mass,phi,theta, timeInterval, month)
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
clear
alt = Zero_model(68.0389,2831.6847,27.2155,0.25,15.5,1);

alt = alt(1:20000) + 1191;
alt = alt(1:100:length(alt))
phi = 40.8500; 
theta = -119.1238;
month = 3;
area = 10;
mass = 68.0389 * 1000;
timeInterval = 100;
day = 1;
year = 2018;
hour = 12;
minute = 0;
second = 0;
fileFolder = 'My Test';
latlon = [phi, theta]; %Starting coordniates (Purdue Airport)
latlon1 = latlon;    
heights = 0;
initV = [0 0];
totalD = 0;
dx =0;
dy = 0;
%% Time Step
for i = 1:length(alt)
    i
    time = timeInterval;
    height =  alt(i);    
    [east_west,north_south] = windDirection(height,phi,theta,month,day,year,hour,minute,second,fileFolder);
    second = second + 1;
    if second == 60
        second = 0;
        minute = minute +1;
        if minute == 60
            minute = 0;
            hour = hour +1;
            if hour == 24
                hour = 0;
                day = day+1;
            end
        end
    end
    Xvelocity = east_west;                   % finding x component of wind velocity (m/s) 
    yvelocity = north_south;                   % finding y component of wind velocity (m/s) 
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
    height = height;
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
latlon1;
totalD / 1000;
heights;


