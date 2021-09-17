#!/bin/bash

#conda activate shakemap

for folder in 197*; do
    echo $folder
    
    shake $folder assemble <<< "Scenario assemble 2021-09-09"
    shake $folder model gridxml
    #shake $folder kml gridxml
    #shake $folder model contour
done

for folder in 198*; do
    echo $folder

    shake $folder assemble <<< "Scenario assemble 2021-09-09"
    shake $folder model gridxml

done

