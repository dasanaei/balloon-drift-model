%% Purdue Orbital 
% 
% NOAA Data
%
% 
% Created: October 24, 2018
%          Matthew Popplewell

%% Description of model
% 
% Exports wind speed and direction given inputs of altitude, latitude, and
% longitude from a pre-loaded NOAA .csv file.
% 

%% Function
function [windSpeed, windDir] = NOAAData(fileName, alt, lat, lon)
% inputs :
% alt: altitude (m)
% lat: latitude
% lon: longitude

% outputs: 
% windSpeed: wind speed (m/s)
% windDir: wind direction (degrees from N)
%% Initializations
altBars = [10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000]'; % altitudes in millibars
altMeters = -log(altBars * 100/101325) * 7000;
%% Find Nearest Altitude to NOAA Points
if alt < altMeters(length(altMeters))
    altRange = [0;altMeters(length(altMeters))] % if less than the first value (~92.14m) set to between sea level and first
    altRangeLoc = [length(altMeters) + 1, length(altMeters)] % set beyond, will be used when parsing data set
else
    resultSortAlt = sort(abs(alt - altMeters)); % find difference between each alt and sort by lowest
    altRangeLoc = [find(abs(alt - altMeters) == resultSortAlt(1)), find(abs(alt - altMeters) == resultSortAlt(2))] % location of two closest alts
    altRange = sort(altMeters(altRangeLoc))
end

%% Parse Data for Altitude
startDir1 = (400 + (altRangeLoc(2) - 1) * 398) - 1; % starting row for csvread (upper bound, direction)
startSpd1 = (599 + (altRangeLoc(2) - 1) * 398) - 1; % starting row for csvread (upper bound, speed)
endDir1 = (597 + (altRangeLoc(2) - 1) * 398) - 1; % ending row for csvread (upper bound, direction)
endSpd1 = (796 + (altRangeLoc(2) - 1) * 398) - 1; % ending row for csvread (upper bound, speed)
startDir2 = (400 + (altRangeLoc(1) - 1) * 398) - 1; % starting row for csvread (lower bound, direction)
startSpd2 = (599 + (altRangeLoc(1) - 1) * 398) - 1; % starting row for csvread (lower bound, speed)
endDir2 = (597 + (altRangeLoc(1) - 1) * 398) - 1; % ending row for csvread (lower bound, direction)
endSpd2 = (796 + (altRangeLoc(1) - 1) * 398) - 1; % ending row for csvread (lower bound, speed)

upperBoundDirM = csvread(fileName, startDir1, 0, [startDir1, 0, endDir1, 2]); % direction data for upper bound alt range
lowerBoundDirM = csvread(fileName, startDir2, 0, [startDir2, 0, endDir2, 2]); % direction data for lower bound alt range
upperBoundSpdM = csvread(fileName, startSpd1, 0, [startSpd1, 0, endSpd1, 2]); % speed data for upper bound alt range
lowerBoundSpdM = csvread(fileName, startSpd2, 0, [startSpd2, 0, endSpd2, 2]); % speed data for lower bound alt range

%% Parse Data for Lat/Lon
possibleLats = 40:.25:42.5; % vector of possible latitude values
possibleLons = 239.25:.25:243.5; % vector of possible longitude values

resultSortLat = sort(abs(lat - possibleLats));
latRangeLoc = [find(abs(lat - possibleLats) == resultSortLat(1)), find(abs(lat - possibleLats) == resultSortLat(2))]; % location of two closest lats
latRange = sort(possibleLats(latRangeLoc)) % two closest NOAA lat values to current lat

resultSortLon = sort(abs(lon - possibleLons));
lonRangeLoc = [find(abs(lon - possibleLons) == resultSortLon(1)), find(abs(lon - possibleLons) == resultSortLon(2))]; % location of two closest Lons
lonRange = sort(possibleLons(lonRangeLoc)) % two closest NOAA lon values to current lon

dataLatLonLoc = zeros(1,4); % vector to store row locations of lat, lon values in NOAA data; ["lowerBoundLat- lowerBound Lon", "lowerBoundLat- upperBound Lon", "upperBoundLat- lowerBound Lon", "upperBoundLat- lowerBound Lon"]

for c = 1:length(upperBoundDirM) % find upper and lower bound lon row locations for lower bound lat
   if upperBoundDirM(c, 2) == latRange(1)
       for v = 0:length(possibleLons) - 1
           if upperBoundDirM(c+v, 1) == lonRange(1);
               if dataLatLonLoc(1) == 0
                   dataLatLonLoc(1) = c+v;
               end
           end
           if upperBoundDirM(c+v, 1) == lonRange(2);
              if dataLatLonLoc(2) == 0
                   dataLatLonLoc(2) = c+v;
              end
           end
       end
   end
   if upperBoundDirM(c, 2) == latRange(2) % find upper and lower bound lon row locations for upper bound lat
       for v = 0:length(possibleLons) - 1
           if upperBoundDirM(c+v, 1) == lonRange(1);
              if dataLatLonLoc(3) == 0
                  dataLatLonLoc(3) = c+v;
              end
           end
           if upperBoundDirM(c+v, 1) == lonRange(2);
              if dataLatLonLoc(4) == 0
                  dataLatLonLoc(4) = c+v;
               end
           end
       end
   end
end

%% Interpolate Data
% Upper Bound Direction
lowerBoundLatDirU = Interpolate(lonRange(1), lon, lonRange(2), upperBoundDirM(dataLatLonLoc(1), 3), upperBoundDirM(dataLatLonLoc(2), 3));
upperBoundLatDirU = Interpolate(lonRange(1), lon, lonRange(2), upperBoundDirM(dataLatLonLoc(3), 3), upperBoundDirM(dataLatLonLoc(4), 3));

upperBoundDir = Interpolate(latRange(1), lat, latRange(2), lowerBoundLatDirU, upperBoundLatDirU)

% Lower Bound Direction
lowerBoundLatDirL = Interpolate(lonRange(1), lon, lonRange(2), lowerBoundDirM(dataLatLonLoc(1), 3), lowerBoundDirM(dataLatLonLoc(2), 3));
upperBoundLatDirL = Interpolate(lonRange(1), lon, lonRange(2), lowerBoundDirM(dataLatLonLoc(3), 3), lowerBoundDirM(dataLatLonLoc(4), 3));

lowerBoundDir = Interpolate(latRange(1), lat, latRange(2), lowerBoundLatDirL, upperBoundLatDirL)

% Upper Bound Speed
lowerBoundLatSpdU = Interpolate(lonRange(1), lon, lonRange(2), upperBoundSpdM(dataLatLonLoc(1), 3), upperBoundSpdM(dataLatLonLoc(2), 3))
upperBoundLatSpdU = Interpolate(lonRange(1), lon, lonRange(2), upperBoundSpdM(dataLatLonLoc(3), 3), upperBoundSpdM(dataLatLonLoc(4), 3))

upperBoundSpd = Interpolate(latRange(1), lat, latRange(2), lowerBoundLatSpdU, upperBoundLatSpdU)

% Lower Bound Speed
lowerBoundLatSpdL = Interpolate(lonRange(1), lon, lonRange(2), lowerBoundSpdM(dataLatLonLoc(1), 3), lowerBoundSpdM(dataLatLonLoc(2), 3))
upperBoundLatSpdL = Interpolate(lonRange(1), lon, lonRange(2), lowerBoundSpdM(dataLatLonLoc(3), 3), lowerBoundSpdM(dataLatLonLoc(4), 3))

lowerBoundSpd = Interpolate(latRange(1), lat, latRange(2), lowerBoundLatSpdL, upperBoundLatSpdL)

% Final Values
windSpeed = Interpolate(altRange(1), alt, altRange(2), lowerBoundSpd, upperBoundSpd);
windDir = Interpolate(altRange(1), alt, altRange(2), lowerBoundDir, upperBoundDir);