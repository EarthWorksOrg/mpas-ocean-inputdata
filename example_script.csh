#

# This is a sample script for ocean and seaice generation. In this example
#   the inputfiles from /glade are put in the 'inputdata' subdirectory.
#   A 'workdata' subdirectory is used for intermediate files, and a
#   'outputdata' subdirectory is used for the files the MPAS ocean and
#   seaice may read in EarthWorks.

# This example will build files for the 60k quasi-uniform grid. Adapt for your needs.

# general glade files
setenv BATHYMETRY_FILE_15s inputdata/SRTM15_plus_earth_relief_15s_original.nc
setenv REFBOTTOMDEPTH inputdata/refBottomDepth
setenv LEVITUS_TEMPERATURE_FILE inputdata/PotentialTemperature.100levels.Levitus.EN4_1900estimate.200813.nc
setenv LEVITUS_SALINITY_FILE inputdata/Salinity.100levels.Levitus.EN4_1900estimate.200813.nc
setenv TRANSECT_TEMPLATE_FILE inputdata/QU60_mocBasinsAndTransects20210623.nc
setenv SALINITY_PATH inputdata/salinity'

# target grid specific input files
setenv GLOBAL_BASE_MESH_FILE inputdata/cami_01-01-2000_00Z_mpasa60_L58_CFSR_c240905.nc
setenv OCEAN_INPUT_RST_FILE inputdata/omip_120.mpaso.rst.0311-01-01_00000.nc
setenv ICE_INPUT_RST_FILE inputdata/omip_120.mpassi.rst.0311-01-01_00000.nc

# intermediate files
setenv BATHYMETRY_FILE_WORK workdata/bathy.nc
setenv INTERPOLATION_FILE_WORK workdata/interp_weights.nc

# output files for EarthWorks execution
setenv OCEAN_MESH_FILE outputdata/ocean.QU.060km.251222.nc
setenv OCEAN_ESMFMESH_FILE outputdata/oQU060_ESMFmesh.251222.nc
setenv OCEAN_GRAPH_FILE outputdata/mpas-o.graph.info.QU060.251222

# optional output files
setenv ICE_RST_FILE outputdata/omipv2.5_060.mpassi.rst.0001-01-01_0000.nc
setenv OCEAN_RST_FILE outputdata/omipv2.5_060.mpaso.rst.0001-01-01_0000.nc
setenv OCEAN_RSTIC_FILE outputdata/oceanIC.QU.060km.251222-2000-01-01.nc
setenv SALINITY_FORCING_FILE outputdata/sss.monthlyClimatology_251222.oQU060.nc
setenv TRANSECT_BASIN_FILE outputdata/QU60_moc_BasinsAndTransects20251222.nc


# compilation commands (these here are for gfortran on a Mas PRO)
setenv F90 gfortran -ffree-line-length-none
setenv F90_net gfortran -ffree-line-length-none -I/opt/homebrew/include/ -L/opt/homebrew/lib -lnetcdff -lnetcdf -lhdf5setenv 
setenv F90_net_pnet mpif90 -ffree-line-length-none -I/opt/homebrew/include/ -I${HOME}/include  -L/opt/homebrew/lib -L${HOME}/lib -lnetcdff -lnetcdf -lhdf5 -lpnetcdf
setenv OMP_FLAG -fopenmp