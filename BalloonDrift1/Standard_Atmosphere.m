function[GeopotentialAlt, DesiredPressure, DesiredDensity, DesiredTemperature] = Standard_Atmosphere(GeometricAlt, Units)
% This function takes as input the geometric altitude desired
% and a 0 if you want to use SI units, anything else if you want to use
% English units. Input in feet if in English, meters if in SI.

%% Establishing base values for layers
BaseAltitudes = [0,11000, 25000, 47000, 53000, 79000, 90000, 105000]; % meters
Slopes = [-6.5 * 10^(-3), 0, 3 * 10 ^ (-3), 0, -4.5 * 10 ^ (-3), 0, 4 * 10^(-3),0]; 
BaseTemps = [288.16, 216.66, 216.66, 282.66, 282.66, 165.66, 165.66, 225.66]; % Kelvin
BasePressures = [1.01325e5, 0, 0, 0, 0, 0, 0, 0]; % N/m^2
BaseDensities = [1.2250, 0, 0, 0, 0, 0, 0, 0]; % kg/m^3
R = 287.058;

%% Other miscellaneous constants
g0 = 9.8;
REarth = 6.3781 * 10 ^ 6;
R = 287;

%% Adjust for units
if Units == 0
    GeometricAlt = GeometricAlt;
else
    GeometricAlt = GeometricAlt / 3.28084;
end
    
%% Calculate the Geopotential altitude given the geometric altitude as input
GeopotentialAlt = REarth / (REarth + GeometricAlt) * GeometricAlt;

%% Discover which layer the given altitude is in
for i = 1:8
    if GeopotentialAlt >= BaseAltitudes(i)
        Layer = i; 
    end
end
%if the layer is the base
if (GeometricAlt == 0)
    DesiredPressure =  1.01325e5;
    DesiredDensity =  1.2250;
    DesiredTemperature = 288.16;
    return
end
%% Calculate successive base layer values until reaching the 
%  layer that the altitude is in at which point you
%  calculate that altitude's properties
for i = 1:Layer
    % Skip the first iteration -- we already have values for this 
    if i == 1
    
    % In other iterations, calculate the base values for layer
    else
        % Isothermic
        if Slopes(i) == 0
            BasePressures(i) = BasePressures(i - 1) * (BaseTemps(i) / BaseTemps(i - 1)) ^ (- g0 / (Slopes(i - 1) * R));
            BaseDensities(i) = BaseDensities(i - 1) * (BaseTemps(i) / BaseTemps(i - 1)) ^ (- ((g0 / (Slopes(i - 1) * R)) + 1));
        % Gradient
        else
            BasePressures(i) = BasePressures(i - 1) * exp(- g0 / (R * BaseTemps(i)) * (BaseAltitudes(i) - BaseAltitudes(i - 1)));
            BaseDensities(i) = BaseDensities(i - 1) * BasePressures(i) / BasePressures(i - 1);
        end
    end
    % Once we reach the desired layer and have calculated its base values,
    % we then use them to find the desired altitude's properties
    if i == Layer
        % Isothermic
        if Slopes(i) == 0
            T = BaseTemps(i);
            DesiredPressure = BasePressures(i) * exp(- g0 /(R * T) * (GeopotentialAlt - BaseAltitudes(i)));
            DesiredDensity = BaseDensities(i) * DesiredPressure / BasePressures(i);
        % Gradient
        else
            T = BaseTemps(i) + Slopes(i) * (GeopotentialAlt - BaseAltitudes(i));
            DesiredPressure = BasePressures(i) * ((T / BaseTemps(i)) ^ (- g0 / (Slopes(i) * R)));
            DesiredDensity = BaseDensities(i) * (T / BaseTemps(i)) ^ (- (g0 / (Slopes(i) * R) + 1));
        end
    end
end

%% The following changes the units based on the input
if Units == 0 
    DesiredPressure = DesiredPressure;
    DesiredDensity = DesiredDensity;
    GeopotentialAlt = GeopotentialAlt;
    % fprintf('Geopotential Altitude = %1.3f m\nPressure = %1.5f N/m^2\nDensity = %1.5f kg/m^3\n', GeopotentialAlt, DesiredPressure, DesiredDensity);
else 
    DesiredPressure = DesiredPressure * 0.0208854342;
    DesiredDensity = DesiredDensity * 0.00194032033;
    GeopotentialAlt = GeopotentialAlt * 3.28084;
    %fprintf('Geopotential Altitude = %1.3f ft\nPressure = %1.5f lb/m^2\nDensity = %1.5f slugs/ft^3\n', GeopotentialAlt, DesiredPressure, DesiredDensity);
end

DesiredTemperature = DesiredPressure / DesiredDensity / R;

return