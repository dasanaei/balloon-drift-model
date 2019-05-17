function [y2] = Interpolate(x1,x2,x3,y1,y3)
%Interpolate Interpolates for y2 using points (x1,y1) and (x3,y3)
y2 = (x2-x1)*(y3-y1)/(x3-x1) + y1;
end

