# solar
Tiny, pure python solar module for handling the variety of calculations that are commonly required 
to make sense of the complex relationship between a point on the earths surface and the sun

**There are several other libraries that do the same things and more**
One of them, [pysolar](https://pysolar.readthedocs.io/en/latest/), does a good job comparing
itself to other modules. This module is most like pysolar, and that module 
offers more long term support.


## What it can do.

Started with replicating astronomical algorithms by Jean Meeus. 
You can read about and see some alternate ways of computing these at
NOAA's ESRL [here](https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html)

Parts have been expanded by more recent work by NREL, namely
[this document](https://www.nrel.gov/docs/fy08osti/34302.pdf)

The core calculator class manages intermediate computations necessary
for returning whatever value is ultimately desired. It supports numpy vectorized computations
for 

* A single point in space, at a single time (scalar).
* A single point in space, at many times (1d).
* A meshgrid of 2d points, at a single time (2d).
* A meshgrid of 2d points at many times (3d).

One of the use cases served is to use some reference raster data source
and compute solar values for every nominal pixel location at a corresponding time.


## Examples

todo