**Creating the files needed for MPAS-Ocean/Seaice in EarthWorks**
Input files needed by EarthWorks:
  Required:
    Mesh and Input file - this has all the mpas metrics plus initial conditions for the ocean. Can be a restart file.
    ESMF mesh file.
    Graph decomposition file.
  Optional:
    Runoff remapping file - needed to properly couple to runoff models
    Surface salinity forcing file - for OMIP style spin up runs.
    Transect and basin file - defines regions for ocean diatnostics.

Initial conditions if not starting from a restart

Ocean temperature and salinity are provided by the mesh and input file. Velocities are zero. The sea surface height is flat.
The model initializes the seaice to be disks of uniform thickness surrounding each pole. This is placed atop the flat ocean.

**Overview of steps to creating these files:**

Start with a global base mesh file - this can be an atmospheric initial condition file or one of the mesh files found at https://mpas-dev.github.io/atmosphere/atmosphere_meshes.html
