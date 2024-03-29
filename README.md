# Abstract

This program converts the administrative boundary data (GeoJSON file) downloaded from http://nlftp.mlit.go.jp to meshed data.
The output file is numpy npz file containing following fields, 
- lons: 2-dimensional array presents the grid points for longitude axis 
- lats: 2-dimensional array presents the grid points for latitude axis 
- data: 2-dimensional array presents the administrative by the index number (0, 1, 2, ...) 
- name: the administrative name included at its corresponding index number of data

In order to run this program, your system is necessary to support displaying Japanese on the console.

# General Usage

```
python mesh_N03.py -j {GeoJSON} -m {Sea/Land mask} -l {AdminLevel} -o {output} [options...]
```

## Mandatory arguments

  ### -h, --help            
  
  show this help message and exit
  
  ### --version             
  
  show program's version number and exit
  
  ### -j GeoJSON, --json GeoJSON
  
  File name of GeoJSON file
  
  ### -m Sea/Land mask, --mask Sea/Land mask
  
  File name of land/sea mask data as numpy npz file 
  
  The file should contain following fields at least, 
  - lons: 2-dimensional array presents the grid points for longitude axis 
  - lats: 2-dimensional array presents the grid points for latitude axis 
  - data: 2-dimensional array present land (=1) and sea (=0) 
  
  ### -l AdminLevel, --level AdminLevel 
  
  the administrative level, 
  
  
  - N03_001: prefectures
  - N03_002: branches of prefecture (only Hokkaido)
  - N03_003: metropolis and city designated by ordinance 
  - N03_004: cities, towns and villages 
  - N03_007: code number of administrative division 
  
  See http://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N03-v2_2.html
  
  ### -o output, --output output
  
  Output file name as numpy npz file
  
  ## Optional arguments
  
  ### --nearest-distance
  
  A threshold of the nearest neighbor algorithm.
  The nearest neighbor algorithm will work with this threshold more than zero.
  Default is inf. Unit is km.
  The algorithm will be used when a grid point is not located in the administrative region.
  This threshold should be set to zero if the mask data fully filled with 1.

  ### --prio-ordcity

  An ordinance-designated city name is used if the administrative range is a part of the ordinance-designated city. 
  This option works only with '-l N03_004'

  ### --inc-metropolis

  Metropolis names will include in the town and village names.
  This option works only with '-l N03_004'

  ### --avoid-empty003

  avoid enpty name (e.g. city name)
  Because city name without ordinance-designated city has no information in N03_003 level,
  the output names without this option with '-l N03_003' will be empty for the normal cities.
  Using this option is recommended with '-l N03_003'.
  This option works only with '-l N03_003'

# GeoJSON data

The GeoJSON data can be downloaded at http://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N03-v2_2.html.
After the downloading and unziping, you can find the geojson file (.geojson).

# Sea/Land mask data

Sea/land mask data is necessary for this program as a numpy npz file.
This mask data will also provide the grid point spacing for the output file.

The mask data must contain following fiels,
 - lons: 2-dimensional array presents the grid points for longitude axis 
 - lats: 2-dimensional array presents the grid points for latitude axis 
 - data: 2-dimensional array present land (=1) and sea (=0) 
  
 If you have no mask data, you can use an array of which all the lements are filled with 1, meaning land.
 An example to create such mask data is shown as follow.
 
 ```
 >>> import numpy as np
 >>> lons = np.array([110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120], dtype=np.float)
 >>> lats = np.array([30, 31, 32, 33, 34, 35], dtype=np.float)
 >>> lons, lats = np.meshgrid(lons, lats)
 >>> data = np.ones([lats.size, lons.size], dtype=np.int)
 >>> np.savez("mask.npz", lons=lons, lats=lats, data=data)
 ```
 
