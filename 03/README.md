## additional benefits

Managing voids with script tools brings along a lot of benefits:
* reduction of time wasted on manual comparing
* reduction of human errors
* data consistency
* reliable base for additional scripts - examples: 
    * void tagging 
    * void dimensioning


### scripted tagging

Below we see a script run to tag the specific void family. 
Other than a regular tag-all, this script checks the height 
of the voids to be tagged and places the tag always with a 
pre-defined distance either above or below 
(respectively left or right) of a wall.
Different tag types can be chosen for specific tags
This leaves only graphical corrections needed to be done 
manually.

![updates](rvt_void_tagging.gif)


### scripted dimensioning

Here we see scripted dimensioning. Similar to the script 
above, it takes into account the height of the void. 
The voids get dimensioned along the wall direction with a 
pre-defined distance grouped by the void height.
This leaves only the ends of the dimension strings to be 
attached to specific grids left needed to be done manually.

![updates](rvt_void_dimensioning.gif)

