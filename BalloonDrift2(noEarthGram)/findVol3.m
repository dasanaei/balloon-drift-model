function [  Vol, diameter ] = findVol3(M_helium, T_g, V_max, alt)
%findVol
%   Finds the balloon volume for the given altitude and payload
%   mass range.  1.796


R = 2077.2;               %J/kgK specific gas const of helium


[~, P_a,~,~] = Standard_Atmosphere(alt, 0);

% quadratic equation to solve for the volume at each altitude
delta_p = M_helium * R * T_g / V_max - P_a;
delta_p = g*(rho_a - rho_g) * .834 * (.75 * M_gas / (pi * rho_g))^(1/3);
Pg = P_a - delta_p;
rho_g = Pg / (R * T_g);
diameter = 2.296 * ((.75 * M_helium) / (pi * rho_g))^(1/3); %Equation from Balloon Ascent Manual
Vol = (diameter / 1.424)^3;                                 %Equation from BalloonAscent: 3-D symmulator 







end
