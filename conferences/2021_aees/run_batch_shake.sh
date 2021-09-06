#!/bin/bash

#conda activate shakemap

for folder in 201*; do
    echo $folder
    
    shake $folder assemble <<< "Scenario assemble 2021-09-03"
    shake $folder model gridxml
    #shake $folder kml gridxml
    #shake $folder model contour
done

for folder in 202*; do
    echo $folder

    shake $folder assemble <<< "Scenario assemble 2021-09-03"
    shake $folder model gridxml

done

