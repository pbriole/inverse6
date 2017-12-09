# inverse6
Modelling of earthquake slip by inversion of GPS and InSAR data assuming homogenous elastic medium
<a href="https://doi.org/10.5281/zenodo.1098391"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1098391.svg" alt="DOI"></a>




### <a name="Objectives"></a>Objectives

This programme is an update of inverse3.f that was not able to ingest InSAR LOS (line of sight) observations
The associated files relate to the 2017 Kos earthquake (see https://www.emsc-csem.org/Files/event/606346/Kos_report_30-7-2017.pdf)


### <a name="installation"></a>Installation

The programmes are in Fortran 77. 


### <a name="Usage"></a>Usage

Using inverse6
1. Compile the programme (with gfortran or f77). Nota: the dimensions can be modified to fit with a given case
2. Three input files needed:
  - kosn14.txt = input file for the case of north dipping fault, koss.txt = input file for the case of south dipping fault
  - kosh.txt = horizontal displacements at the GPS stations
  - kosvda.txt = vertical motion at the GPS stations + LOS (line of sight motion) observed in descending and ascending InSAR
    the number of observations is defined at line 785 of the source code and can be changed there
3. Three output files are produced: DIFF, FAIL, DEPL

Using direct4a (this code produces a synthetic interferogram)
1. Compile the programme
2. Two input files are needed:
  - the file parameters (kosn14.txt or kosss.txt depending on the selected fault plane)
  - a file containing the parameters of the future synthetic interferogram and the LOS (line of sight) angle:
    kosa.grd (for ascending), kosd.grd (for descending)
3. The programme asks for the name of a output image file: you can put 'out.raw'. This output file will be in .raw binaty format
  (1 bype per pixel)


### <a name="authors"></a>Authors

* Pierre Briole (briole@ens.fr)


### <a name="references"></a>References

*  Okada 1985 for the model, Tarantola and Valette 1982 for the inversion, Briole et al., 1986 for the initial version of the software https://www.scopus.com/record/display.uri?eid=2-s2.0-0022879463&origin=inward&txGid=23cc1f7dfe0c13e5a489090f02957406


### <a name="license"></a>License

Copyright 2017 Pierre Briole - CNRS/ENS

Licensed under the GPL v3.0
