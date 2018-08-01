# mesh_N03

This program converts the administrative boundary data (GeoJSON file) downloaded from http://nlftp.mlit.go.jp to meshed data.
The output file is numpy npz file containing following fields, 
- lons: 1-dimensional array presents the grid points for longitude axis 
- lats: 1-dimensional array presents the grid points for latitude axis 
- data: 2-dimensional array presents the administrative by the index number (0, 1, 2, ...) 
- name: the administrative name included at its corresponding index number of data

# Usage

```
python mesh_N03.py -j {GeoJSON} -m {Seal/Land mask} -l {AdminLevel} -o {output} [--disable-nearest]
```

## optional arguments

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
  
  ### --disable-nearest     
  
  Use nearest neighbor algorithm or not 
  
  This algorithm will be used when a grid point is not located in the administrative region. 
  You should use this option if the mask data filled with 1 is assigned.
