## How to use BalloonDrift1.m
1) Obtain payload mass, volume, balloon mass, helium mass, total mass, and impacted surface area of the balloon.
2) Determine starting location and find the coordinates (with as much significant figures as possible)
3) Go to http://weather.uwyo.edu/upperair/sounding.html and change the parameters to the desired date and time of launch (only 2 days in the future is available)
4) Click on the closest sounding station on the map, and copy the data provided on a .txt file.
5) Remove all columns or rows with missing data (so matlab can read it as a matrix)
6) Using zero_pressure.m and the parameters found in step 1, obtain an array of altitudes. Know what the time interval is between each measured altitude. 
7) Plug in all the required inputs in that you have gathered above. 
8) Go to http://www.gpsvisualizer.com/map_input?form=googleearth 
9) Change the Output File Type to .kml
10) Change the altitude mode to “extruded”
11) Upload the KMLFILE.csv that is created when running the code
12) Download and open the new .kml file (must have google earth)





