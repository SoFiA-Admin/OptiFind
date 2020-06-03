# OptiFind

Python wrapper script for catalogue-based source finding with SoFiA 2

## Usage

```
optifind.py <par_file> <source_list> <r_spat> <r_spec> [<sofia_exe>]
```

 ## Arguments:
 
 * `<par_file>`     Name of the SoFiA 2 control parameter file.
 * `<source_list>`  Name of the source catalogue file.
 * `<r_spat>`       Spatial radius of sub-region in pixels.
 * `<r_spec>`       Spectral radius of sub-region in channels.
 * `<sofia_exe>`    Optional name of the SoFiA 2 executable. Default: sofia.

## Description

 OptiFind serves as a wrapper script around the SoFiA 2 source finding pipeline
 to allow source finding  on multiple sub-regions of a data cube centred on the
 positions from a user-supplied source catalogue. The user will need to provide
 a template SoFiA 2 parameter file that will be used in each source finding run
 spawned by OptiFind.

 The template file must specify the input data cube to be searched and can also
 define an output file name which will be used as the base name for all output.
 In addition,  the user must specify a source catalogue containing the position
 of each source to be searched. The catalogue must specify the world coordinate
 position  of each source  on a separate line.  The following,  comma-separated
 parameters must be supplied with each source:

 `id, coord_1, coord_2, coord_3, ...`

 Here, id is a unique source ID  that will be used as an identifier for output
 products, while coord_n denotes the coordinates in each dimension of the data
 cube.  The coordinate must be specified in the raw units of the data cube  as
 specified in the FITS header.  A coordinate value must be given for each axis
 of the cube in the correct order.

 For example, if a cube has four axes (right ascension, declination, frequency
 and Stokes), then four coordinate values need to be provided for each source,
 and they must  be given  in the native units  specified  in the header,  e.g.
 degrees for right ascension  and Hz for frequency.  An example  for a 3D data
 cube with RA, declination and velocity axis might look like this:
 
 ```
 # Example catalogue
 Source 1, 180.7, 62.0, 1300000.0
 Source 2, 180.9, 62.3, 1200000.0
 ```
 
 Here, RA and declination are in degrees,  while velocity is specified in m/s,
 which are the default units defined by the FITS standard.

 Separate output catalogues and products will be created for each SoFiA 2 run.
 They will be named either  `optifind` + suffix  or  output.filename + suffix,
 depending on whether an output file name  was defined  in the parameter file.
 The suffix  will be an underscore  followed by the source ID  provided in the
 catalogue file.

 In addition to the individual output catalogues from each run,  OptiFind will
 also create a single, merged catalogue called `optifind_merged_catalogue.txt`
 in the  same  output directory.   Note that  this feature  is currently  only
 available for plain-text ASCII catalogues, but not for XML or SQL catalogues.
 
 ## Copyright and licence

 Copyright (C) 2020 Tobias Westmeier

 OptiFind is free software: you can redistribute it and/or modify it under the
 terms of the GNU General Public License as published by the Free Software
 Foundation, either version 3 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License  along with
this program. If not, see http://www.gnu.org/licenses/.
