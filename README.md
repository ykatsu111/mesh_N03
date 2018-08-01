# Abstract

This program converts the administrative boundary data (GeoJSON file) downloaded from http://nlftp.mlit.go.jp to meshed data.
The output file is numpy npz file containing following fields, 
- lons: 1-dimensional array presents the grid points for longitude axis 
- lats: 1-dimensional array presents the grid points for latitude axis 
- data: 2-dimensional array presents the administrative by the index number (0, 1, 2, ...) 
- name: the administrative name included at its corresponding index number of data

# General Usage

```
python mesh_N03.py -j {GeoJSON} -m {Sea/Land mask} -l {AdminLevel} -o {output} [--disable-nearest]
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
  - lons: 1-dimensional array presents the grid points for longitude axis 
  - lats: 1-dimensional array presents the grid points for latitude axis 
  - data: 2-dimensional array present land (=1) and sea (=0) 
  
  If you do not have the mask data, an array filled with 1 can be assigned as the data field. 
  If so, note that you should not use the nearest neighbor algorithm (set --disable-nearest).
  
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
  
  ### --disable-nearest     
  
  Use nearest neighbor algorithm or not 
  
  This algorithm will be used when a grid point is not located in the administrative region. 
  You should use this option if the mask data filled with 1 is assigned.

# GeoJSON data

The GeoJSON data can be downloaded at http://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N03-v2_2.html.
After the downloading and unziping, you can find the geojson file (.geojson).

# Sea/Land mask data

Sea/land mask data is necessary for this program as a numpy npz file.
This mask data will also provide the grid point spacing for the output file.

The mask data must contain following fiels,
 - lons: 1-dimensional array presents the grid points for longitude axis 
 - lats: 1-dimensional array presents the grid points for latitude axis 
 - data: 2-dimensional array present land (=1) and sea (=0) 
  
 If you have no mask data, you can use an array of which all the lements are filled with 1, meaning land.
 An example to create such mask data is shown as follow.
 
 ```
 >>> import numpy as np
 >>> lons = np.array([110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120], dtype=np.float)
 >>> lats = np.array([30, 31, 32, 33, 34, 35], dtype=np.float)
 >>> data = np.ones([lats.size, lons.size], dtype=np.int)
 >>> np.savez("mask.npz", lons=lons, lats=lats, data=data)
 ```
 
