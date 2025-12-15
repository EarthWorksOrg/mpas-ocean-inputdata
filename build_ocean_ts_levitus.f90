program build_ocean_ts_levitus

! compile with mpi, hdf5, netcdf and pnetcdf libraries


use netcdf
use pnetcdf, only : nf90mpi_open, nf90mpi_close, nf90mpi_inq_dimid, &
                    nf90mpi_put_var_all, nf90mpi_inq_varid, nf90mpi_inquire_dimension, &
                    nf90mpi_put_var, nf90mpi_write=>nf90_write, nf90mpi_sync, nf90mpi_get_var_all
use mpi

implicit none


!!!
! Declarations

   ! master ocean data files
   character(len=128) :: ocnfile, sfile, tfile
   !resolution specific variables
   integer caseres

   ! netcdf-related variables
   integer :: rc, ncid, ncidf
   integer :: dimid, varid, attid, dimidcf, dimidvf, dimidef, &
              dimidt, dimidme, dimid2, dimidvd, dimidme2, dimidsl,   &
              dimid1, dimidcat, dimidily, dimidsly, varid2, dimidl
   integer :: natts, ncids, ncidt
   character(len=80) :: attname

   ! ocean data arrays
   real*4, dimension(0:361,0:181,100,1:2) :: ocdat
   integer, parameter :: indxs = 1, indxt = 2

   ! grid variables
   ! vertices
   integer :: ncellsc, ncellsf, maxlinks, nverticesc, nverticesf, nedgesc, nedgesf
   
   integer, allocatable, dimension(:) :: indexcella2o, indexverta2o, indexedgea2o
   integer, allocatable, dimension(:) :: indexcello2a, indexverto2a, indexedgeo2a
   integer, allocatable, dimension(:) :: nedgesonCell, indexedgeid, indexedgeid_a
   integer, allocatable, dimension(:) :: indexcellid
   integer :: temcells(6), ncellstem
   
   real*8 r8min, r8max, wrk1, wrk2, areatem, sphere_radius
   integer :: ne, nlen
   

   ! generic output arrays
   real*8, allocatable, dimension(:) :: real8_1d
   real*8, dimension(100) :: refbottomdepth
   real*8, allocatable, dimension(:) :: bottomdepth, loncello, latcello

   
   
   ! loop indices
   integer :: n, na, m, id, k, iii, kk, indxo, l

   
   ! ocean data arrays
   real*8 pi, d2r, coslat, latloc, lonloc, lon_bb, z_b, lat_bb, dist
   integer imax, imax1, imin, imin1, jmax, jmin, ncount, ntot, nedges, item
   real*8 latmax, latmin, lonmax, lonmin, lon, x_a, y_a, z_a
   integer ncorr, dimid3, dimid4, dimid5
   
   
   ! mpi variables
   integer :: err, rank, nprocs, info, nctem
           integer(kind=MPI_OFFSET_KIND) G_NY, myOff, block_start, &
                                        global_nx, global_ny
          integer(kind=MPI_OFFSET_KIND) start(2), count(2)
          integer(kind=MPI_OFFSET_KIND) malloc_size, sum_size, ntem
  

!!!
! Executable code begins
          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! set an MPI-IO hint to disable file offset alignment for
          ! fixed-size variables
          call MPI_Info_create(info, err)
          call MPI_Info_set(info, "nc_var_align_size", "1", err)


   pi = 4._8 * atan(1._8)
   d2r = pi/180.
!!!
! preliminaries

   print*,'enter the ocean mesh file name'
   read(5,'(a128)') ocnfile

   ! read ocean data
   tfile = '/Volumes/disk11/mpas-ocean/inputdata/PotentialTemperature.100levels.Levitus.EN4_1900estimate.200813.nc'
   sfile = '/Volumes/disk11/mpas-ocean/inputdata/Salinity.100levels.Levitus.EN4_1900estimate.200813.nc'
   rc = nf90_open(sfile,nf90_nowrite,ncids)
   rc = nf90_open(tfile,nf90_nowrite,ncidt)
   rc = nf90_inq_varid(ncids,'SALT',varid)
      rc = nf90_get_var(ncids,varid,ocdat(1:360,1:180,1:100,indxs))
   rc = nf90_inq_varid(ncidt,'TEMP',varid)
      rc = nf90_get_var(ncidt,varid,ocdat(1:360,1:180,1:100,indxt))
   ocdat(0,:,:,:) = ocdat(360,:,:,:)
   ocdat(361,:,:,:) = ocdat(1,:,:,:)
   ocdat(:,0,:,:) = ocdat(:,1,:,:)
   ocdat(:,181,:,:) = ocdat(:,180,:,:)
   rc = nf90_close(ncids)
   rc = nf90_close(ncidt)
   

          call MPI_Info_free(info, err)
   rc = nf90mpi_open(mpi_comm_world,trim(ocnfile),nf90mpi_write,info,ncidf)
print*,'dbg create ',rc
   if(rc .ne. 0) stop 'troubles opening the ocean mesh file - stopping'

!!!!!! define and copy the dimensions   
   rc = nf90mpi_inq_dimid(ncidf,'nCells',dimid)
   print*,'inqdimidncells',rc, dimid
   rc = nf90mpi_inquire_dimension(ncidf,dimid,len=ntem)
      ncellsf = ntem
   print*,'nCells',rc,nCellsf


print*,'some allocates 1'
!!!!!! allocate necessary work variables
   allocate(real8_1d(ncellsf))
   allocate(latcello(ncellsf))
   allocate(loncello(ncellsf))
   allocate(bottomdepth(ncellsf))

   !lonCell
   rc = nf90mpi_inq_varid(ncidf,'lonCell',varid)
   rc = nf90mpi_get_var_all(ncidf,varid,loncello)
   loncello = loncello * 180./pi
   where (loncello < 0.) loncello = loncello + 360.
   where (loncello > 360.) loncello = loncello - 360.

   !latCell
   rc = nf90mpi_inq_varid(ncidf,'latCell',varid)
   rc = nf90mpi_get_var_all(ncidf,varid,latcello)
   latcello = latcello * 180./pi
   
   rc = nf90mpi_inq_varid(ncidf,'refBottomDepth',varid)
   rc = nf90mpi_get_var_all(ncidf,varid,refbottomdepth)
print*,'refbottomdepth',rc
   rc = nf90mpi_inq_varid(ncidf,'bottomDepth',varid)
   rc = nf90mpi_get_var_all(ncidf,varid,bottomdepth)
print*,'bottomdepth',rc
      
print*,'start temperature'
   !double temperature(Time, nCells, nVertLevels) ;
   !        temperature:long_name = "temperature" ;
   !        temperature:units = "C" ;
   !METHOD: Levitus 
   rc = nf90mpi_inq_varid(ncidf,'temperature',varid)
   do k = 1,100
      do n = 1,ncellsf
         real8_1d(n) = 0.0
         if(k==1) then
            real8_1d(n) = get_ts(loncello(n),latcello(n),k,indxt)
         else
            if(bottomDepth(n) > refbottomdepth(k-1) ) then
               real8_1d(n) = get_ts(loncello(n),latcello(n),k,indxt)
            else
               real8_1d(n) = -1.e34
            endif
         endif
      enddo
       rc = nf90mpi_put_var_all(ncidf,varid,real8_1d, (/k,1,1/)*1_mpi_offset_kind, (/1,ncellsf,1/)*1_mpi_offset_kind)
print*,'temperature ',k,rc
   enddo
   

print*,'start salinity'
   !double salinity(Time, nCells, nVertLevels) ;
   !        salinity:long_name = "salinity" ;
   !        salinity:units = "grams salt per kilogram seawater" ;
   !METHOD: Levitus 
   rc = nf90mpi_inq_varid(ncidf,'salinity',varid)
   do k = 1,100
      do n = 1,ncellsf
         real8_1d(n) = 0.0
         if(k==1) then
            real8_1d(n) = get_ts(loncello(n),latcello(n),k,indxs)
         else
            if(bottomDepth(n) > refbottomdepth(k-1) ) then
               real8_1d(n) = get_ts(loncello(n),latcello(n),k,indxs)
            else
               real8_1d(n) = -1.e34
            endif
         endif
      enddo
       rc = nf90mpi_put_var_all(ncidf,varid,real8_1d, (/k,1,1/)*1_mpi_offset_kind, (/1,ncellsf,1/)*1_mpi_offset_kind)
print*,'salinity ',k,rc
   enddo


! close up the output file
   rc = nf90mpi_close(ncidf)
   
   call mpi_finalize(err)

contains

   real*8 function get_ts(lon,lat,k,index)
      real*8, intent(in) :: lon, lat
      integer, intent(in) :: k, index
      
      real*8 dx, dy, tlon, tlat, dxr, dyr
      integer i, j
      
      i = lon + 0.5
      j = lat + 90.5
      tlon = real(i)-0.5
      tlat = real(j)-90.5
      
      dx = lon - tlon
      dy = lat - tlat
      dxr = 1.-dx
      dyr = 1.-dy
      get_ts = ocdat(i  ,j  ,k,index) * dxr * dyr + &
               ocdat(i+1,j  ,k,index) * dx  * dyr + &
               ocdat(i  ,j+1,k,index) * dxr * dy  + &
               ocdat(i+1,j+1,k,index) * dx  * dy
      
   end function get_ts   
   
end program
