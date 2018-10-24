%% Purdue Orbital 
% 
% NOAA Data
%
% 
% Created: October 24, 2018

%% Description of model
% 
% Exports wind speed and direction given inputs of altitude, latitude, and
% longitude from a pre-loaded NOAA .csv file.
% 

%% Function
function [windSpeed, windDir] = NOAAData(alt,lat, lon)
% inputs :
% alt: altitude (m)
% lat: latitude
% lon: longitude

% outputs: 
% windSpeed: wind speed (m/s)
% windDir: wind direction (degrees from N)

%% Initializations

