function [delta_m]=find_mass(alt1,alt2,m)
%find_mass
%calculate the change in mass of ballon

[~, P_a, rho1, T_a] = Standard_Atmosphere(alt1, 0);
[~, P_a, rho2, T_a] = Standard_Atmosphere(alt2, 0);
delta_m = (1-rho2/rho1)*m


end