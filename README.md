# aircraft-wing-analysis

A MATLAB code which analyses an aerofoil section using 2D potential flow panel method and integral boundary layer equation solver. This produces aerofoil performance metrics, most notably lift-to-drag ratio. A wing analysis surface generator (WASG) was created to design our own aerofoil design, optimised for both high and low Reynold's numbers.

RUN:
- main code in foil.m
- choose which section to run when asked, see Geometry folder for names of sections
- plot streamlines in streamline.m
- create new wing section with WASG.m
- change number of panels, reynolds number, range of angle incidences to iterarate, in Parfiles folder
