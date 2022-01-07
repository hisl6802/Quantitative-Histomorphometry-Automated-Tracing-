# Quantitative-Histomorphometry-Automated-Tracing

# NOTE: Future updates of this script will incorporate a color wheel, input for scale factor and number of different colors.

## Overview

This script allows users to automatically trace and calculate distances of thick line boundaries. Additionally, this script incorporates manual selection of non-enclosed thick lines using a geodesic distance approach between selected point on the line. Currently this script is built for quantitative histomorphometry analysis using a white thick line for tracing full boundaries (on a boundary image), a white for Label 1 (Labels image), and a yellow for Label 2. Additionally, any locations where Label 1 and Label 2 overlap a blue color is used to trace. The blue trace is incorporated into the distance calculations for both Labels (1 and 2). **After analysis you will be prompted to save a .csv file, which you can assign a name to.**

# IMPORTANT NAMING CONVENTION NOTE


Boundary file should have the following appended to it: "_#_BonePerimeter_#.png"
Labels file should have the following appended to it: "_#_Labels.png"

**If you would like to update this convention you can update lines 63-65 of the script.**


### RECOMMENDED UPDATES

#### Scale factor: update on line 60

#### Output Labels: update lines 72-78

#### Boundary trace color: update on line 100

#### Label 1 color: Update line 518

#### Label 2 color: Update line 883

#### Dual Label color: Update line 525 & line 890

