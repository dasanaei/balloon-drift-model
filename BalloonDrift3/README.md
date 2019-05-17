	HOW TO USE BalloonDrift1.m
Obtain payload mass, volume, balloon mass, helium mass, total mass, and impacted surface area of the balloon.
Determine starting location and find the coordinates (with as much significant figures as possible)
Go to http://weather.uwyo.edu/upperair/sounding.html and change the parameters to the desired date and time of launch (only 2 days in the future is available)
click on the closest sounding station on the map, and copy the data provided on a .txt file.
remove all columns or rows with missing data (so matlab can read it as a matrix)
using zero_pressure.m and the parameters found in step 1, obtain an array of altitudes. Know what the time interval is between each measured altitude. 
Plug in all the required inputs in that you have gathered above. 
Go to http://www.gpsvisualizer.com/map_input?form=googleearth 
change the Output File Type to .kml
Change the altitude mode to “extruded”
upload the KMLFILE.csv that is created when running the code
download and open the new .kml file (must have google earth)





