## high level overview

nestedRvtDurchbruchCutter is very configurable, 
but also assumes a configured setup:

* main model
    * linked container
        * nested void models content from 3rd-party to be regularly updated.
          ordered per model or floor or discipline

It further requires in the main model sets of levels that match the script 
configuration, e.g. set of raw levels set of finish levels.

To run it under [RevitPythonShell](https://github.com/architecture-building-systems/revitpythonshell) 
or [pyRevit](https://github.com/eirannejad/pyRevit) , select a nested void model.
The script will create a temporary mirror of the link and search through the model
for the configured voids (remote voids). It will represent the void in the main 
model (local void), so they can be safely annotated (tags/dimensions) and the 
annotation for stable voids can be kept.

When searching through the link for voids, there are three standard cases:

* **obsolete**: uncut and comment local void as obsolete
* **existing**: update local void data only, re-cut if needed
* **new void**: create, cut local void, timestamp & write data


#### obsolete
If a local void only exists in the main model but not as remote void in 
the external link, the local void gets commented as obsolete and receives 
a timestamp, when it got considered as obsolete. 
It will also be uncut from elements it was cutting.


#### existing
If local void exists at the exact location as a remote void in link model, 
the void stays geometrically untouched, but receives a data update: 
Parameter values are re-written according to the linked model.
A re-cut evaluation is performed, to see if the void needs to re-cut. 
e.g. new elements added to its surrounding after the void's creation.


#### new void
If a remote void is found in the nested link model, but not as local void in 
main model, a new accordingly sized instance of the configured void family 
is created, with a new size type if needed. configured parameter values are 
mirrored from the link to the new instance. 

If the linked void has rotation, a note is written to the void comments 
parameter. (as currently there is not rotation performed by the script, 
also see known limitations)

A collision check with preconfigured sets (via type name white and black lists) 
of wall, floor, roof and structural-framing elements is performed. 
If the local void does not intersect with the solid of any existing element, 
a warning is written to the void comments parameter.
The intersecting elements will be cut, if it matches the white listing for 
its category. 
If the void cuts elements of multiple categories (e.g. walls **and** floors) 
a note is written to the void comments parameter. 
Depending on the element being cut by the local void, the visibility parameters 
of the void family are set accordingly: wall cut symbol, floor cut symbol, etc..

Distances to above, below or insert levels are written to the void as configured.

3d model lines ("antennas") of the void family are extended, to ensure 
visibility in floor plans and reflected ceiling plans.
Antennas for voids cutting walls are extended over the level height. 
Antennas for voids cutting floors or roofs are extended until raw floors 
above and below. The latter will additionally receive a parameter value 
of its closest raw floor, so voids can easily be filtered.


#### known limitations
* There need to be at least a raw and finish level below the lowest void, 
so it can be properly placed.
* Void rotation is not implemented (yet).


-----

nestedRvtDurchbruchCutter sent out to * 
20200422 10:34 as base for collaboration.
