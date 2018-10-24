function [east_west,north_south] = windDirection(height,phi,theta,month,day,year,hour,minute,second,fileFolder)
changeFile(height,phi,theta,month,day,year,hour,minute,second,fileFolder);
EarthGRAMExecuter();
[east_west,north_south] = readOutput(fileFolder);
if east_west == 'NaN'
    east_west = 0;
end
if north_south == 'NaN'
    north_south = 0;
end
end