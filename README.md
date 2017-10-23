TINRAD & TINGEOFIX readme file
==============================================================================

When using any of the scripts described, please cite our article:

> Moreno, H.A., Ogden, F.L., and Alvarez, L.V. XXXX. Unstructured-Mesh Terrain 
> Analysis and Incident Solar Radiation for Continuous Hydrologic Modeling in
> Mountain Watersheds. Computers and Geosciences, Issue X, Vol X, pp XX--XX.

--------


This is the readme file for the **TINGEOFIX** and **TINRAD** R scripts that
produce Triangular Irregular Networks (TIN), analyze their topographic
properties, and estimate the total incident radiation accounting for direct,
diffuse, atmospheric backscattered, and albedo-reflected components. The
repository also contains demo data for the Green River basin ready to run
TINGEOFIX.



TINGEOFIX
==============================================================================

If you are starting from scratch and have set of `X,Y,Z` coordinates in a text
file it is perfect for running the **TINGEOFIX** script. Before running, the
following R libraries need to be installed in your R version:

 * rgeos
 * sp
 * maptools
 * spatstat
 * fields
 * colorRamps
 * plotrix
 * moments
 
You can install them all using the following command:

    install.packages(c('rgeos','sp','maptools','spatstat','fields','colorRamps','plotrix','moments'))
 
Additionally, you will need to provide a shapefile (i.e. `.dbf`, `.shp`, `.shx`,
and `.prj`) containing the watershed divide or polygon of interest to constrain
the tessellation and avoid getting weird TIN edges.

The inputs to TINGEOFIX are:

 * **WD**:         A string to set working directory

 * **divisoria**:  A string with the name of the watershed divide shapefile to 
                   constrain the tessellation.

 * **vertices**:   A string with the name of a textfile with the TIN triangle 
                   vertices containing `X,Y,Z` (`Z` is elevation). Coordinates 
                   in projected system in meters.

 * **maxres**:     Is a numeric with the length units of largest expected 
                   triangle edge within the TIN to find nearest three nodes.
                   Rule of thumb is if triangle coordinates were derived from a 
                   30m DEM then use 30m to start.

 * **tolerancia**: Is a numeric tolerance value for distance when simplifying
                   the tessellation achieved with Delaunay. Begin with tolerance 
                   of 0 and if processing is taking long then increase 
                   accordingly.

 * **slices**:     Is a numeric argument to divide the visible horizon from each 
                   triangle (e.g. 4, 8, 16, 32, 64, 128). 
                   The more slices the more accurate the calculations.

 * **RE**:         Is numeric value for the Earth's radius at the equator GRS80 
                   ellipsoid in meters.

 * **makeplots**:  Logic value. `1` to make plots, `0` to not make plots.
                   Not making the plots may save time.

The outputs from to TINGEOFIX are:

 * **clipTIN**:     A Delaunay triangulation array.

 * **ordersky**:    A vector of sky view fractions according to `n` slice 
                    directions

 * **Shading**:     A vector of elements (rows) slices (cols) according to `n` 
                    shading directions. 

 * **TIN_Centers**: A vector of `X,Y,Z` with the central coordinates of each TIN 
                    element.

 * **vertices2**:   A vector of TIN vertices ready to ingest in TINRAD.

 * **nu**:          A vector containing the `nux`, `nuy`, and `nuz` components 
                    of the normal unit vector to each TIN element.



TINRAD
==============================================================================

After you have run **TINGEOFIX**, the outputs are used as inputs for **TINRAD**.
This script creates an incident radiation matrix per hour (cols) and TIN element
(rows) and its components.

The Inputs to TINRAD are in GMT standard time:

 * **iyear**:    Numeric. Initial year of simulation 
 * **imonth**:   Numeric. Initial month of simulation
 * **iday**:     Numeric. Initial day of simulation
 * **ihour**:    Numeric. Initial hour of simulation
 * **imin**:     Numeric. Initial minute of simulation
 * **isec**:     Numeric. Initial second of simulation
 * **fyear**:    Numeric. Final year of simulation 
 * **fmonth**:   Numeric. Final month of simulation
 * **fday**:     Numeric. Final day of simulation
 * **fhour**:    Numeric. Final hour of simulation
 * **fmin**:     Numeric. Final minute of simulation
 * **fsec**:     Numeric. Final second of simulation
 * *sconstant*:  Numeric. solar constant in w/m2
 * *RE*:         Numeric. Earth's radius at equator in meters GRS80 ellipsoid
 * *Lambda00*:   Numeric. Reference meridian in degrees. 
 * *FalseEast*:  Numeric. False East of desired coordinate system
 * *FalseNorth*: NUmeric. False North of  desired coordinate system

 * *Ta*:         Vector. Time step air temperatures in Kelvin. This can be a 
                 vector with  values at each triangle center of TINCENTERS.

 * *RH*:         Vector. Time step air relative humidities in percent. This can 
                 be a vector with values at each triangle center of TINCENTERS.

 * *albedo*:     Vector. Time step albedo values. It can be one basin-average 
                 value for all the simulation or one for each element and time 
                 step.

The outputs from to TINRAD are:

 * **finalrad**: A numeric matrix of TIN element (rows) x time step (cols) with 
                 the final estimated incident radiation values including direct, 
                 self and remote shading, sky view fraction, diffuse,
                 atmospheric backscattered and albedo-reflected.
 
 * **Kdif**:     A numeric matrix of TIN element (rows) x time step (cols) with 
                 the final estimated diffuse radiation.

 * **Kbs**:      A numeric matrix of TIN element (rows) x time step (cols) with 
                 the final estimated atmospheric backscattered radiation.

 * **Kal**:      A numeric matrix of TIN element (rows) x time step (cols) with 
                 the final estimated albedo-reflected radiation.
