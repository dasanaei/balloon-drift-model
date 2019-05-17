%Author: Mike Bailey

%Reads the h1 value in NameRef.txt, change the h1 value, 
%then run the GRAM2016 application 

%To run this code, change the fileFolder variable to the location of 
%the test to the one you want to use. Also do this in the response.txt file
%found in this directory. After that, change the starting variables to
%numbers that apply to your specific simulation. The output will be found 
%in output.txt in which ever folder you changed fileFolder to.

clear;clc; %clears the workspace and command window
fileFolder = 'My Test'; %folder that includes necessary files DO NOT INCLUDE THE FINAL BACKSLASH
x = 1; %index value for arrays
height = 0; %starting altitude in kilometers
time = 0; %starting time in seconds
phi = 40.422; %starting latitude in degrees
theta = -86.931; %starting longitude in degrees
delta_h = 0.5; %change in height in kilometers
delta_phi = 0.25*delta_h; %change in latitude in degrees
delta_theta = -0.125*delta_h; %change in longitude in degrees
changeFile(height,phi,theta,10,10,2018,10,0,30,fileFolder); %runs the changeFile function to initialize EarthGRAM

while x <= 10
    EarthGRAMExecuter(); %runs EarthGRAM
    [ew(x),ns(x)] = readOutput(fileFolder); %finds the west to east and south to north wind speeds
    fprintf('h = %.3f km\n', height); 
    fprintf('ew = %.2f m/s\n',ew(x));
    fprintf('ns = %.2f m/s\n',ns(x));
    fprintf('lat = %.3f\n',phi);
    fprintf('long = %.3f\n\n',theta);
    height = height + delta_h; %iterates the heights
    phi = phi + delta_phi; %iterates the latitude
    theta = theta + delta_theta; %iterates the longitude
    changeFile(height,phi,theta,10,10,2018,10,0,30,fileFolder); %changes the reference file for EarthGRAM
    x = x+1; %iterates the index value
end