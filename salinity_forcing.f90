PROGRAM salinity_forcing

! compile with netcdf, pnetcdf

! This program will interpolate WOA13v2 values of salinity in the top layer
!  to the MPAS grid and create a salinity forcing file with monthly climatology
! Adapted from woa13v2_sub.f90

   USE netcdf, only : nf90_open, nf90_inquire, nf90_close, nf90_write, nf90_nowrite,   &
                      nf90_get_var, nf90_inq_varid, nf90_inq_dimid,                     &
                      nf90_inquire_dimension, nf90_put_var!, nf90_strerror
   USE pnetcdf, only : nf90mpi_create, nf90mpi_close, nf90mpi_unlimited => nf90_unlimited,          &
                       nf90mpi_write=>nf90_write, nf90mpi_inq_dimid, nf90mpi_inquire_dimension, &
                       nf90mpi_inq_varid, nf90mpi_get_var_all, nf90mpi_put_var_all,   &
                       nf90mpi_def_var, NF90mpi_DOUBLE=>NF90_DOUBLE, nf90mpi_char=>nf90_char,   &
                       NF90mpi_GLOBAL=>NF90_GLOBAL, nf90mpi_enddef, nf90mpi_put_att,    &
                       nf90mpi_def_dim, NF90MPI_CLOBBER=>NF90_CLOBBER,     &
                       NF90MPI_64BIT_DATA=>NF90_64BIT_DATA
   USE mpi

   IMPLICIT NONE

   ! VARIABLE DECLARATIONS

   ! MPAS salinity_forcing file to be created
   CHARACTER(LEN=256) :: filename

   ! MPAS file with grid information to be read
   CHARACTER(LEN=256) :: gridfile, spath
   
   CHARACTER(LEN=64) :: xtime
   INTEGER :: dimid_ncells, dimid_time, dimid_strlen

   ! NetCDF-related variables for MPAS file
   INTEGER :: rc, ncid
   INTEGER(KIND=MPI_OFFSET_KIND) :: pcount(2), pstart(2)

   ! mpi variables for pnetcdf
   integer :: err, rank, nprocs, info

   ! MPAS grid, s variables
   INTEGER :: ncells
   INTEGER*8 :: pdim
   REAL*8, DIMENSION(:), ALLOCATABLE :: lonCell, latCell
   REAL*8, DIMENSION(:), ALLOCATABLE :: layertop, layerbase                        
   ! The following is just the surface layer, all months done one point at a time
   REAL*8, DIMENSION(12) :: salinity
   INTEGER :: mi  ! the character location of month
   CHARACTER(LEN=2) :: mon

   ! WOA13-related variables
   INTEGER :: nlon, nlat
   REAL*4, DIMENSION(:), ALLOCATABLE :: lon_in, lat_in
   REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: s_in
   INTEGER :: res_in
   CHARACTER(len=256) :: sfilename
   REAL :: delta_ll

   ! other netcdf variables
   INTEGER :: varid, dimid

   ! loop indices
   INTEGER :: k, m
   
   ! cartesian grid locations on the unit sphere
   REAL*8, DIMENSION(:,:), ALLOCATABLE :: xyz_m
   REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: xyz_in
   REAL*8 :: dtor

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! BEGIN EXECUTABLE CODE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! set an MPI-IO hint to disable file offset alignment for
          ! fixed-size variables
          call MPI_Info_create(info, err)
          call MPI_Info_set(info, "nc_var_align_size", "1", err)

      dtor = 4.*atan(1.) / 180.


   ! inquire for the path for the salinity input data
   PRINT *,'enter the salinity input data path'
   READ(*,'(A256)') spath

   ! inquire for the input resolution to be used
   PRINT *,'enter 1 for 1 degree WOA data, 4 for 0.25 degree'
   READ(*,*)res_in

   ! read in the MPAS filenames
   PRINT *,"Enter the name of the MPAS grid file "
   READ(*,'(A256)') gridfile
   PRINT *,"Enter the name of the file to create "
   READ(*,'(A256)') filename


   ! read in the WOA data
   SELECT CASE(res_in)
   CASE(1)
      sfilename = trim(spath)//'/res_1.0/woa18_decav_sxx_01.nc'
      mi = 67
      delta_ll = 1.0
   CASE(4)
      sfilename = trim(spath)//'/res_0.25/woa18_decav_sxx_04.nc'
      mi = 68
      delta_ll= 0.25
   CASE DEFAULT
      PRINT *,'invalid input resolution selected - stopping'
      STOP
   END SELECT

   DO m = 1,12
      WRITE(UNIT=mon,FMT='(i2.2)')m
      sfilename(mi:mi+1) = mon(1:2)
      print*,'sfilename',trim(sfilename)
      rc = nf90_open(TRIM(sfilename),nf90_nowrite,ncid)
      print*,'open',rc,m
       
      IF(m == 1) THEN
         ! get dimensions, allocate input salinity
         rc = nf90_inq_dimid(ncid,'lon',dimid)
print*,'inqdimid lon ',rc
         rc = nf90_inquire_dimension(ncid,dimid,len=nlon)
print*,'inquiredimension lon ',rc
         rc = nf90_inq_dimid(ncid,'lat',dimid)
print*,'inqdimid lat ',rc
         rc = nf90_inquire_dimension(ncid,dimid,len=nlat)
print*,'inquiredimension lat ',rc

         ALLOCATE(lon_in(nlon))
         ALLOCATE(lat_in(nlat))
         ALLOCATE(s_in(nlon,nlat,12))
         
         rc = nf90_inq_varid(ncid,'lon',varid)
print*,'inqvarid lon ',rc
         rc = nf90_get_var(ncid,varid,lon_in)
print*,'getvar lon ',rc
         rc = nf90_inq_varid(ncid,'lat',varid)
print*,'inqvarid lat ',rc
         rc = nf90_get_var(ncid,varid,lat_in)
print*,'getvar lat ',rc
      ENDIF
   
      rc = nf90_inq_varid(ncid,'s_an',varid)
print*,'inqvarid s ',rc
      rc = nf90_get_var(ncid,varid,s_in(:,:,m),start=(/1,1,1,1/),count=(/nlon,nlat,1,1/))
print*,'getvar s ',rc
      rc = nf90_close(ncid)
   ENDDO


   ! get the grid variables
   CALL get_grid()
   print*,'get_grid done'
      
   ! create the output file
   rc = nf90mpi_create(mpi_comm_world,'salinity_forcing/'//TRIM(filename),IOR(NF90MPI_CLOBBER, NF90MPI_64BIT_DATA),info,ncid)
print*,'nf90mpi_create',rc
   rc = nf90mpi_def_dim(ncid,'nCells',ncells*1_mpi_offset_kind,dimid_ncells)
print*,'nf90mpi_defdim ncells',rc
   rc = nf90mpi_def_dim(ncid,'Time',nf90mpi_unlimited*1_mpi_offset_kind,dimid_time)
print*,'nf90mpi_defdim time',rc
   rc = nf90mpi_def_dim(ncid,'StrLen',64_mpi_offset_kind,dimid_strlen)
print*,'nf90mpi_defdim strlen',rc
   rc = nf90mpi_def_var(ncid,'xtime',nf90mpi_char, (/dimid_strlen,dimid_time/), varid)
print*,'nf90mpi_defvar xtime',rc
   rc = nf90mpi_def_var(ncid,'surfaceSalinityMonthlyClimatologyValue',nf90mpi_double, &
                          (/dimid_ncells,dimid_time/), varid)
print*,'nf90mpi_defvar salinity',rc
   rc = nf90mpi_put_att(ncid, varid, 'long_name',    &
    "Salinity is restored toward this field at a rate controlled by salinityPistonVelocity.")
print*,'nf90mpi_putatt longname',rc
   rc = nf90mpi_put_att(ncid, varid, 'units', "PSU")
print*,'nf90mpi_putatt units',rc
   rc = nf90mpi_enddef(ncid)
print*,'nf90mpi_enddef units',rc
   rc = nf90mpi_inq_varid(ncid,'xtime',varid)
print*,'nf90mpi_inqvarid xtime',rc
   pcount = (/64,1/)
   xtime = '0000-mm-15_00:00:00                                             '
   DO m = 1,12
      pstart = (/1,m/)
      WRITE(UNIT=xtime(6:7),FMT='(i2.2)')m
      rc = nf90mpi_put_var_all(ncid,varid,xtime,start=pstart,count=pcount)
print*,'nf90mpi_putvarll xtime',rc
   ENDDO
   
print *,'start process_level_k'
   ! for the surface layer, interpolate WOA data to the MPAS grid and write 
   CALL process_level_k()
print *,'done process_level_k'

   ! close the output file
   rc = nf90mpi_close(ncid)
print *,'done nf90mpi_close',rc

   call mpi_finalize(err)


   CONTAINS

   !-------------------------------

   SUBROUTINE get_grid()
   ! This routine will determine the dimensions of the the MPAS file,
   !  allocate MPAS variables, and read the grid variables

      REAL :: rtod
      INTEGER :: n, i, j

      ! get grid dimensions
      rc = nf90_open(TRIM(gridfile),nf90_nowrite,ncid)
      rc = nf90_inq_dimid(ncid,'nCells',dimid)
      rc = nf90_inquire_dimension(ncid,dimid,len=nCells)
      PRINT *,'MPAS dimensions nCells ', ncells

      ! allocate grid variables
      ALLOCATE(loncell(ncells))
      ALLOCATE(latcell(ncells))
      ALLOCATE(xyz_m(3,ncells))
      ALLOCATE(xyz_in(3,nlon,nlat))


      ! read in grid variables
      rc = nf90_inq_varid(ncid,'lonCell',varid)
      rc = nf90_get_var(ncid,varid,loncell)
      rc = nf90_inq_varid(ncid,'latCell',varid)
      rc = nf90_get_var(ncid,varid,latcell)

      rc = nf90_close(ncid)
      
      ! compute cartesian locations of mpas grid
      DO n = 1,ncells
         CALL ll_xyz_mpas(loncell(n),latcell(n),xyz_m(:,n))
      ENDDO

      ! compute cartesian locations of input grid
      DO j = 1,nlat
         DO i = 1,nlon
            CALL ll_xyz_in(lon_in(i),lat_in(j),xyz_in(:,i,j))
         ENDDO
      ENDDO

      ! convert latcell and loncell from radians to degrees
      rtod = 180. / (4.*atan(1.))
      latcell = latcell * rtod
      loncell = loncell * rtod

      ! put loncell in the range [-180 180]
      DO n = 1,ncells
         IF(loncell(n) > 180.) loncell(n) = loncell(n) - 360.
      ENDDO

   END SUBROUTINE get_grid

   !-------------------------------

   SUBROUTINE process_level_k()
   ! This routine will interpolate the WOA top layer salinity
   !  data to the MPAS grid.

      ! Argument list variables : none

      ! local variables
      INTEGER(KIND=MPI_OFFSET_KIND) :: pcount(2), pstart(2)
      INTEGER :: l, n, i, j, m
      INTEGER, DIMENSION(4) :: i_x, j_y       ! bounding woa lon and lat indices
      INTEGER, DIMENSION(4) :: mask
      REAL :: wgt, wsum, dist
      INTEGER :: nbad, mcount
      REAL :: d_w, d_e, d_n, d_s, dsum, dmin
      LOGICAL :: cell_exists
      REAL*8 :: dotloc, dotmax
      INTEGER :: imax, jmax

      ! BEGIN EXECUTABLE CODE

print*,'ncells ',ncells
      nbad = 0
      ! interpolate
      DO n = 1,ncells
         salinity = 0.0
         pcount = (/1,12/)
         pstart = (/n,1/)
      
         ! find the lat, lon indices of the bounding WOA lat/lon box
         j_y(1) = (latcell(n) + 90 + 0.5 * delta_ll) / delta_ll
         j_y(3) = j_y(1) + 1
         j_y(1) = MAX(1,j_y(1))
         j_y(3) = MIN(nlat, j_y(3))
         i_x(1) = (loncell(n) + 180 + 0.5 * delta_ll) / delta_ll
         i_x(3) = i_x(1) + 1
         IF(i_x(1) < 1) i_x(1) = nlon
         IF(i_x(3) > nlon) i_x(3) = 1
         j_y(2) = j_y(3)
         j_y(4) = j_y(1)
         i_x(2) = i_x(1)
         i_x(4) = i_x(3)
            
if(n==1) then
   print*,'cell loc',loncell(n),latcell(n)
   print*,'cell bounds',i_x(1),j_y(1),s_in(i_x(1),j_y(1),1)
   print*,'cell bounds',i_x(2),j_y(2),s_in(i_x(2),j_y(2),1)
   print*,'cell bounds',i_x(3),j_y(3),s_in(i_x(3),j_y(3),1)
   print*,'cell bounds',i_x(4),j_y(4),s_in(i_x(4),j_y(4),1)
endif

               ! eliminate box corners where there is no data defined
               mask = 1
               DO m = 1,4
                  IF (s_in(i_x(m),j_y(m),1) > 999.) mask(m) = 0
               ENDDO
               mcount = SUM(mask)
                              
               IF(mcount > 0) THEN                  
                  ! do the interpolation with the existing defined bounding box points
                  SELECT CASE(mcount)
                  ! using sin(delta) as an approximation to delta so no longitude correction
                  !  is needed at domain edge
                  CASE(4)   ! bilinear interpolation
                     ! d_* is the approximate distance to n, s, w, and e box edges respectively
                     d_w = SIN(dtor*(loncell(n)-lon_in(i_x(1))))
                     d_e = SIN(dtor*(lon_in(i_x(3))-loncell(n)))
                     dsum = d_e + d_w
                     dsum = 1./dsum
                     d_e = d_e * dsum
                     d_w = d_w * dsum
                     d_s = SIN(dtor*(latcell(n)-lat_in(j_y(1))))
                     d_n = SIN(dtor*(lat_in(j_y(3))-latcell(n)))
                     dsum = d_n + d_s
                     IF(j_y(1) .ne. j_y(3)) THEN
                        dsum = 1./dsum
                        d_n = d_n * dsum
                        d_s = d_s * dsum
                     ELSE 
                        d_n = 0.5
                        d_s = 0.5
                     ENDIF
                     salinity = (s_in(i_x(1),j_y(1),:)*d_e*d_n +    &
                                 s_in(i_x(2),j_y(2),:)*d_e*d_s +    &
                                 s_in(i_x(3),j_y(3),:)*d_w*d_s +    &
                                 s_in(i_x(4),j_y(4),:)*d_w*d_n )
                                 
                  CASE(1,2,3) ! distance-weighted interpolation
                     wsum = 0.0
                     DO m = 1,4
                        IF(mask(m)==1) THEN
                           dist = SQRT( SIN( dtor*(loncell(n)-lon_in(i_x(m))) )**2 +  &
                                        SIN( dtor*(latcell(n)-lat_in(j_y(m))) )**2 )
                           wgt = 1. / (1.e-7 + dist)
                           wsum = wsum + wgt
                           salinity = salinity + s_in(i_x(m),j_y(m),:) * wgt
                        ENDIF
                     ENDDO
                     salinity = salinity / wsum
                  END SELECT
               
if(n==1) then
   print*,'salinity bilin',salinity
endif
            ELSE
               ! find the nearest neighbor in WOA layer ktop out one row/column
               cell_exists = .FALSE.
               dmin = 1.
               
               ! West edge
               i = i_x(1) - 1
               IF(i < 1) i = i + nlon
               DO j = MAX(1,j_y(1)),MIN(nlat,j_y(3))
                  IF(s_in(i,j,1) < 999.) THEN
                     dist = SQRT( SIN( dtor*(loncell(n)-lon_in(i)) )**2 +  &
                                  SIN( dtor*(latcell(n)-lat_in(j)) )**2 )
                     IF (dist < dmin) THEN
                           salinity = s_in(i,j,:)
                        dmin = dist
                        cell_exists = .TRUE.
                     ENDIF
                  ENDIF
               ENDDO
               
               ! East edge
               i = i_x(1) + 1
               IF(i > nlon) i = i - nlon
               DO j = MAX(1,j_y(1)),MIN(nlat,j_y(3))
                  IF(s_in(i,j,1) < 999.) THEN
                     dist = SQRT( SIN( dtor*(loncell(n)-lon_in(i)) )**2 +  &
                                  SIN( dtor*(latcell(n)-lat_in(j)) )**2 )
                     IF (dist < dmin) THEN
                           salinity = s_in(i,j,:)
                        dmin = dist
                        cell_exists = .TRUE.
                     ENDIF
                  ENDIF
               ENDDO
               
               ! North edge
               j = j_y(3) + 1
               IF(j <= nlat) THEN
                  DO m = 1,3,2
                  i = i_x(m)
                     IF(s_in(i,j,1) < 999.) THEN
                        dist = SQRT( SIN( dtor*(loncell(n)-lon_in(i)) )**2 +  &
                                     SIN( dtor*(latcell(n)-lat_in(j)) )**2 )
                        IF (dist < dmin) THEN
                           salinity = s_in(i,j,:)
                           dmin = dist
                           cell_exists = .TRUE.
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
               
               ! South edge
               j = j_y(1) - 1
               IF(j > 0) THEN
                  DO m = 1,3,2
                  i = i_x(m)
                     IF(s_in(i,j,1) < 999.) THEN
                        dist = SQRT( SIN( dtor*(loncell(n)-lon_in(i)) )**2 +  &
                                     SIN( dtor*(latcell(n)-lat_in(j)) )**2 )
                        IF (dist < dmin) THEN
                           salinity = s_in(i,j,:)
                           dmin = dist
                           cell_exists = .TRUE.
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
               
               IF(.not.cell_exists) THEN
                  ! do a global nearest neighbor search
                  ! choose the point that produces the maximum dot product with the MPAS cell
                  nbad = nbad + 1
                  dotmax = -1.
                  DO j = 1,nlat
                     DO i = 1,nlon
                        IF ( s_in(i,j,1) < 999.) THEN
                           dotloc = DOT_PRODUCT(xyz_m(:,n),xyz_in(:,i,j))
                           IF ( dotloc > dotmax) THEN
                              dotmax = dotloc
                              imax = i
                              jmax = j
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
                  salinity = s_in(imax,jmax,:)
               ENDIF
if(n==1) then
   print*,'salinity nn',salinity
endif
            ENDIF
         ! write the variables  to the MPAS file

         rc = nf90mpi_inq_varid(ncid,'surfaceSalinityMonthlyClimatologyValue',varid)
         rc = nf90mpi_put_var_all(ncid,varid,salinity,start=pstart,count=pcount)

      ENDDO

      IF(nbad > 0) THEN
         PRINT*,'Global nearest neighbor search for ',nbad,' cells.'
      ENDIF

   END SUBROUTINE process_level_k
   
   !-------------------------------

   SUBROUTINE ll_xyz_mpas(lon,lat,xyz)
   ! compute unit sphere cartesian location of an MPAS grid cell 
   ! lat and lon input are radians and 8-byte reals
   REAL*8, INTENT(IN) :: lon,lat
   REAL*8, INTENT(OUT), DIMENSION(3) :: xyz
      
   xyz(1) = cos(lon) * cos(lat)
   xyz(2) = sin(lon) * cos(lat)
   xyz(3) = sin(lat)
   
   END SUBROUTINE ll_xyz_mpas
   
   !-------------------------------

   SUBROUTINE ll_xyz_in(lon,lat,xyz)
   ! compute unit sphere cartesian location of an input WOA grid cell 
   ! lat and lon input are degrees and 4-byte reals
   REAL*4, INTENT(IN) :: lon,lat
   REAL*8, INTENT(OUT), DIMENSION(3) :: xyz
   
   REAL*8 :: lonloc, latloc
   
   lonloc = lon  * dtor
   latloc = lat  * dtor
   xyz(1) = cos(lonloc) * cos(latloc)
   xyz(2) = sin(lonloc) * cos(latloc)
   xyz(3) = sin(latloc)
   
   END SUBROUTINE ll_xyz_in
   
   !-------------------------------

END PROGRAM salinity_forcing
