%You shouldn't need to touch this

function [east_west, north_south] = readOutput(fileFolder)
fileLocation = strcat(fileFolder,'\output.txt');
fid = fopen(fileLocation);

lineNum = 0;
lineData = '0';
while lineNum ~= 21
    lineNum = lineNum + 1;
    lineData=fgetl(fid);
end

x = 83;

while x >= 2
    x = x - 1;
    if (x < 58) || (x > 73)
        lineData(x) = [];
    elseif (x < 57) || (x > 65)
        if lineData(x) == [' ']
            lineData(x) = [];
        end
    end
end

for x = 1:6
    east_west(x) = lineData(x);
end

for x = 7:12
    z = x - 6;
    north_south(z) = lineData(x);
end

east_west = str2double(east_west);
north_south = str2double(north_south);

fclose(fid);
end