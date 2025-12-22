program build_ice_restart_file

! compile with netcdf, pnetcdf, and, optionally openmp

! setenv OMP_NUM_THREADS 16

use netcdf
use pnetcdf, only : nf90mpi_create, nf90mpi_close, nf90mpi_enddef, nf90mpi_def_dim, &
                    nf90mpi_put_att, nf90mpi_def_var, nf90mpi_put_var_all, nf90mpi_inq_varid, &
                    nf90mpi_put_var, nf90mpi_redef, &
                    nf90mpi_clobber => NF90_CLOBBER, NF90mpi_64BIT_DATA=>NF90_64BIT_DATA, &
                    nf90mpi_unlimited => nf90_unlimited, NF90mpi_GLOBAL=>NF90_GLOBAL, &
                    NF90mpi_DOUBLE=>NF90_DOUBLE, NF90mpi_INT=>NF90_INT, nf90mpi_char=>nf90_char
use mpi

implicit none


!!!
! Declarations

   ! master ocean data files
   character(len=128) :: infile, terpfile, coarsefile, finefile
   !resolution specific variables
   integer caseres

   ! netcdf-related variables
   !   'c' suffix is coarse grid
   !   'f' suffix is fine grid
   integer :: rc, ncid, ncidc, ncidf
   integer :: dimid, varid, attid, dimidcf, dimidvf, dimidef, &
              dimidt, dimidme, dimid2, dimidvd, dimidme2, dimidsl,   &
              dimid1, dimidcat, dimidily, dimidsly, varid2
   integer :: natts
   character(len=80) :: attname

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
   
   ! generic input arrays
   integer, allocatable, dimension(:) :: integer1d_in
   integer, allocatable, dimension(:,:) :: integer2d_in
   real*8, allocatable, dimension(:) :: real8_1d_in
   real*8, allocatable, dimension(:,:) :: real8_2d_in

   ! generic atm arrays
   integer, allocatable, dimension(:) :: integer1d_a
   integer, allocatable, dimension(:) :: integer1d_edge_a
   integer, allocatable, dimension(:) :: integer1d_vert_a
   integer, allocatable, dimension(:,:) :: integer2d_a
   real*8, allocatable, dimension(:) :: real8_1d_a
   real*8, allocatable, dimension(:,:) :: real8_2d_a
   real*8, allocatable, dimension(:) :: real8_vert_a
   real*8, allocatable, dimension(:) :: real8_edge_a

   ! generic output arrays
   integer, allocatable, dimension(:) :: integer1d
   integer, allocatable, dimension(:) :: integer1d_edge
   integer, allocatable, dimension(:) :: integer1d_vert
   integer, allocatable, dimension(:,:) :: integer2d
   real*8, allocatable, dimension(:) :: real8_1d
   real*8, allocatable, dimension(:,:) :: real8_2d
   real*8, allocatable, dimension(:) :: real8_vert
   real*8, allocatable, dimension(:) :: real8_edge
   real*8, allocatable, dimension(:,:) :: real8_2d_edge
   
   
   ! output arrays
   character*64 :: time
   character*16 :: tstring
   
   ! loop indices
   integer :: n, na, m, id, k, iii, kk, indxo, l
   
   ! ocean data arrays
   real*8 pi, d2r, coslat, latloc, lonloc, lon_bb, z_b, lat_bb, dist
   integer imax, imax1, imin, imin1, jmax, jmin, ncount, ntot, nedges, item
   real*8 latmax, latmin, lonmax, lonmin, lon, x_a, y_a, z_a
   integer ncorr, dimid3, dimid4, dimid5
   
   integer start2(2), start3(3), count2(2), count3(3)
   
   ! remapping arrays
   real*4, dimension(:,:), allocatable :: remap_wgts
   integer, dimension(:,:), allocatable :: remap_indices
   real*4, dimension(:,:), allocatable :: remap_wgtsv
   integer, dimension(:,:), allocatable :: remap_indicesv

   
   ! mpi variables
   integer :: err, rank, nprocs, info, nctem
           integer(kind=MPI_OFFSET_KIND) G_NY, myOff, block_start, &
                                        global_nx, global_ny
          integer(kind=MPI_OFFSET_KIND) start(2), count(2)
          integer(kind=MPI_OFFSET_KIND) malloc_size, sum_size
  
        integer, dimension(:), allocatable :: intwork1d
        integer, dimension(:,:), allocatable :: intwork2d
        real*8, dimension(:), allocatable :: r8work1d
        real*8, dimension(:,:), allocatable :: r8work2d

!!!
! Executable code begins
print*,'start'
          call MPI_Init(err)
print*,'start1'
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
print*,'start1.1'
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)
print*,'start2'

          ! set an MPI-IO hint to disable file offset alignment for
          ! fixed-size variables
          call MPI_Info_create(info, err)
          call MPI_Info_set(info, "nc_var_align_size", "1", err)
print*,'start3'


   pi = 4._8 * atan(1._8)
   d2r = pi/180.
!!!
! preliminaries

   !choose the case and define the case-specific variables
   print*,'enter the ocean mesh file name'
   read(5,'(a128)')infile
   print*,'enter the interpolation weights file name'
   read(5,'(a128)')terpfile
   print*,'enter the input seaice restart file name'
   read(5,'(a128)')coarsefile
   print*,'enter the seaice restart file name to be created'
   read(5,'(a128)')finefile
   
   rc = nf90mpi_create(mpi_comm_world,trim(finefile),IOR(NF90_CLOBBER, NF90_64BIT_DATA),info,ncidf)

          call MPI_Info_free(info, err)

   ! open the necessary netcdf files
   rc = nf90_open(trim(infile),nf90_nowrite,ncid)
   rc = nf90_open(trim(coarsefile),nf90_nowrite,ncidc)
   print*,'open coarse restart',rc
      ! read the coarse grid dimensions
      rc = nf90_inq_dimid(ncidc,'nCells',dimid)
      rc = nf90_inquire_dimension(ncidc,dimid,len=ncellsc)
      rc = nf90_inq_dimid(ncidc,'nVertices',dimid)
      rc = nf90_inquire_dimension(ncidc,dimid, len=nverticesc)
      rc = nf90_inq_dimid(ncidc,'nEdges',dimid)
      rc = nf90_inquire_dimension(ncidc,dimid, len=nedgesc)

   allocate(real8_1d_in(ncellsc))

!!!!!! define and copy the dimensions   
   rc = nf90_inq_dimid(ncid,'nCells',dimid)
   rc = nf90_inquire_dimension(ncid,dimid,len=ncellsf)
   rc = nf90_inq_dimid(ncidc,'nCells',dimid)
   rc = nf90_inquire_dimension(ncidc,dimid,len=ncellsc)
   rc = nf90mpi_def_dim(ncidf,'nCells',ncellsf*1_mpi_offset_kind,dimidcf)

   rc = nf90_inq_dimid(ncid,'nEdges',dimid)
   rc = nf90_inquire_dimension(ncid,dimid,len=nedgesf)
   rc = nf90_inq_dimid(ncidc,'nEdges',dimid)
   rc = nf90_inquire_dimension(ncidc,dimid,len=nedgesc)
   rc = nf90mpi_def_dim(ncidf,'nEdges',nedgesf*1_mpi_offset_kind,dimidef)

   rc = nf90_inq_dimid(ncid,'nVertices',dimid)
   rc = nf90_inquire_dimension(ncid,dimid,len=nverticesf)
   rc = nf90_inq_dimid(ncidc,'nVertices',dimid)
   rc = nf90_inquire_dimension(ncidc,dimid,len=nverticesc)
   rc = nf90mpi_def_dim(ncidf,'nVertices',nverticesf*1_mpi_offset_kind,dimidvf)

!ice files
   !!! define and add additional required dimensions
   
   rc = nf90mpi_def_dim(ncidf,'Time',nf90mpi_unlimited*1_mpi_offset_kind,dimidt)
   rc = nf90mpi_def_dim(ncidf,'maxEdges',6*1_mpi_offset_kind,dimidme)
   rc = nf90mpi_def_dim(ncidf,'TWO',2*1_mpi_offset_kind,dimid2)
   rc = nf90mpi_def_dim(ncidf,'vertexDegree',3*1_mpi_offset_kind,dimidvd)
   rc = nf90mpi_def_dim(ncidf,'maxEdges2',12*1_mpi_offset_kind,dimidme2)
   rc = nf90mpi_def_dim(ncidf,'StrLen',64*1_mpi_offset_kind,dimidsl)
   rc = nf90mpi_def_dim(ncidf,'ONE',1*1_mpi_offset_kind,dimid1)
   rc = nf90mpi_def_dim(ncidf,'nCategories',5*1_mpi_offset_kind,dimidcat)
   rc = nf90mpi_def_dim(ncidf,'nIceLayers',7*1_mpi_offset_kind,dimidily)
   rc = nf90mpi_def_dim(ncidf,'nSnowLayers',5*1_mpi_offset_kind,dimidsly)

   rc = nf90_get_att(ncid,NF90_GLOBAL,'sphere_radius',sphere_radius)
   if(sphere_radius < 10.) then
      ! unit sphere - reset sphere_radius to real earth
      sphere_radius = 6371229.
      rc = nf90mpi_put_att(ncidf, NF90mpi_GLOBAL, 'sphere_radius', sphere_radius)
   else
      ! this is the real earth - write then reset to 1 so grid metrics don't change
      rc = nf90mpi_put_att(ncidf, NF90mpi_GLOBAL, 'sphere_radius', sphere_radius)
      sphere_radius= 1.0
   endif

print*,'some allocates 1'
!!!!!! allocate necessary work variables
   allocate(real8_1d(ncellsf))
   allocate(real8_edge(nedgesf))
   allocate(integer1d(ncellsf))
   allocate(integer1d_edge(nedgesf))
   allocate(integer1d_vert(nverticesf))
   allocate(real8_vert(nverticesf))
   
!!!!!! define and copy the existing mesh variables from the culled mesh
   rc = nf90_inq_varid(ncid,'lonCell',varid)
   rc = nf90_get_var  (ncid,varid,real8_1d)
print*,'mpi def loncell',rc
   rc = nf90mpi_def_var(ncidf,'lonCell',nf90mpi_DOUBLE, (/dimidcf/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,real8_1d)
print*,'mpi put loncell',rc
   rc = nf90mpi_redef(ncidf)
print*,'done loncell'      

   rc = nf90_inq_varid(ncid,'latCell',varid)
   rc = nf90_get_var  (ncid,varid,real8_1d)
   rc = nf90mpi_def_var(ncidf,'latCell',nf90mpi_DOUBLE, (/dimidcf/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,real8_1d)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'meshDensity',varid)
   rc = nf90_get_var  (ncid,varid,real8_1d)
   rc = nf90mpi_def_var(ncidf,'meshDensity',nf90mpi_DOUBLE, (/dimidcf/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,real8_1d)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'xCell',varid)
   rc = nf90_get_var  (ncid,varid,real8_1d)
   rc = nf90mpi_def_var(ncidf,'xCell',nf90mpi_DOUBLE, (/dimidcf/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_1d = real8_1d * sphere_radius
   rc = nf90mpi_put_var_all(ncidf,varid,real8_1d)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'yCell',varid)
   rc = nf90_get_var  (ncid,varid,real8_1d)
   rc = nf90mpi_def_var(ncidf,'yCell',nf90mpi_DOUBLE, (/dimidcf/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_1d = real8_1d * sphere_radius
   rc = nf90mpi_put_var_all(ncidf,varid,real8_1d)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'zCell',varid)
   rc = nf90_get_var  (ncid,varid,real8_1d)
   rc = nf90mpi_def_var(ncidf,'zCell',nf90mpi_DOUBLE, (/dimidcf/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_1d = real8_1d * sphere_radius
   rc = nf90mpi_put_var_all(ncidf,varid,real8_1d)
   rc = nf90mpi_redef(ncidf)
print*,'done xyzcell'      

   rc = nf90_inq_varid(ncid,'indexToCellID',varid)
   rc = nf90_get_var  (ncid,varid,integer1d)
   rc = nf90mpi_def_var(ncidf,'indexToCellID',nf90mpi_int, (/dimidcf/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,integer1d)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'areaCell',varid)
   rc = nf90_get_var  (ncid,varid,real8_1d)
   rc = nf90mpi_def_var(ncidf,'areaCell',nf90mpi_DOUBLE, (/dimidcf/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_1d = real8_1d * sphere_radius**2
   rc = nf90mpi_put_var_all(ncidf,varid,real8_1d)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'nEdgesOnCell',varid)
   rc = nf90_get_var  (ncid,varid,integer1d)
   rc = nf90mpi_def_var(ncidf,'nEdgesOnCell',nf90mpi_int, (/dimidcf/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,integer1d)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'lonEdge',varid)
   rc = nf90_get_var  (ncid,varid,real8_edge)
   rc = nf90mpi_def_var(ncidf,'lonEdge',nf90mpi_DOUBLE, (/dimidef/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,real8_edge)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'latEdge',varid)
   rc = nf90_get_var  (ncid,varid,real8_edge)
   rc = nf90mpi_def_var(ncidf,'latEdge',nf90mpi_DOUBLE, (/dimidef/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,real8_edge)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'xEdge',varid)
   rc = nf90_get_var  (ncid,varid,real8_edge)
   rc = nf90mpi_def_var(ncidf,'xEdge',nf90mpi_DOUBLE, (/dimidef/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_edge = real8_edge * sphere_radius
   rc = nf90mpi_put_var_all(ncidf,varid,real8_edge)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'yEdge',varid)
   rc = nf90_get_var  (ncid,varid,real8_edge)
   rc = nf90mpi_def_var(ncidf,'yEdge',nf90mpi_DOUBLE, (/dimidef/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_edge = real8_edge * sphere_radius
   rc = nf90mpi_put_var_all(ncidf,varid,real8_edge)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'zEdge',varid)
   rc = nf90_get_var  (ncid,varid,real8_edge)
   rc = nf90mpi_def_var(ncidf,'zEdge',nf90mpi_DOUBLE, (/dimidef/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_edge = real8_edge * sphere_radius
   rc = nf90mpi_put_var_all(ncidf,varid,real8_edge)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'indexToEdgeID',varid)
   rc = nf90_get_var  (ncid,varid,integer1d_edge)
   rc = nf90mpi_def_var(ncidf,'indexToEdgeID',nf90mpi_int, (/dimidef/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,integer1d_edge)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'nEdgesOnEdge',varid)
   rc = nf90_get_var  (ncid,varid,integer1d_edge)
   rc = nf90mpi_def_var(ncidf,'nEdgesOnEdge',nf90mpi_int, (/dimidef/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,integer1d_edge)
   rc = nf90mpi_redef(ncidf)

print*,'done xyzedge'      

   rc = nf90_inq_varid(ncid,'dvEdge',varid)
   rc = nf90_get_var  (ncid,varid,real8_edge)
   rc = nf90mpi_def_var(ncidf,'dvEdge',nf90mpi_DOUBLE, (/dimidef/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_edge = real8_edge * sphere_radius
   rc = nf90mpi_put_var_all(ncidf,varid,real8_edge)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'dcEdge',varid)
   rc = nf90_get_var  (ncid,varid,real8_edge)
   rc = nf90mpi_def_var(ncidf,'dcEdge',nf90mpi_DOUBLE, (/dimidef/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_edge = real8_edge * sphere_radius
   rc = nf90mpi_put_var_all(ncidf,varid,real8_edge)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'angleEdge',varid)
   rc = nf90_get_var  (ncid,varid,real8_edge)
   rc = nf90mpi_def_var(ncidf,'angleEdge',nf90mpi_DOUBLE, (/dimidef/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,real8_edge)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'lonVertex',varid)
   rc = nf90_get_var  (ncid,varid,real8_vert)
   rc = nf90mpi_def_var(ncidf,'lonVertex',nf90mpi_DOUBLE, (/dimidvf/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,real8_vert)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'latVertex',varid)
   rc = nf90_get_var  (ncid,varid,real8_vert)
   rc = nf90mpi_def_var(ncidf,'latVertex',nf90mpi_DOUBLE, (/dimidvf/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,real8_vert)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'xVertex',varid)
   rc = nf90_get_var  (ncid,varid,real8_vert)
   rc = nf90mpi_def_var(ncidf,'xVertex',nf90mpi_DOUBLE, (/dimidvf/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_vert = real8_vert * sphere_radius
   rc = nf90mpi_put_var_all(ncidf,varid,real8_vert)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'yVertex',varid)
   rc = nf90_get_var  (ncid,varid,real8_vert)
   rc = nf90mpi_def_var(ncidf,'yVertex',nf90mpi_DOUBLE, (/dimidvf/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_vert = real8_vert * sphere_radius
   rc = nf90mpi_put_var_all(ncidf,varid,real8_vert)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'zVertex',varid)
   rc = nf90_get_var  (ncid,varid,real8_vert)
   rc = nf90mpi_def_var(ncidf,'zVertex',nf90mpi_DOUBLE, (/dimidvf/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_vert = real8_vert * sphere_radius
   rc = nf90mpi_put_var_all(ncidf,varid,real8_vert)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'areaTriangle',varid)
   rc = nf90_get_var  (ncid,varid,real8_vert)
   rc = nf90mpi_def_var(ncidf,'areaTriangle',nf90mpi_DOUBLE, (/dimidvf/), varid)
   rc = nf90mpi_enddef(ncidf)
   real8_vert = real8_vert * sphere_radius**2
   rc = nf90mpi_put_var_all(ncidf,varid,real8_vert)
   rc = nf90mpi_redef(ncidf)

print*,'done xyzvertex'      

   rc = nf90_inq_varid(ncid,'indexToVertexID',varid)
   rc = nf90_get_var  (ncid,varid,integer1d_vert)
   rc = nf90mpi_def_var(ncidf,'indexToVertexID',nf90mpi_int, (/dimidvf/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,integer1d_vert)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'cellsOnCell',varid)
   rc = nf90mpi_def_var(ncidf,'cellsOnCell',nf90mpi_INT, (/dimidme,dimidcf/), varid2)
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,6
      rc = nf90_get_var  (ncid,varid,integer1d,(/k,1/),(/1,ncellsf/))
      rc = nf90mpi_put_var_all(ncidf,varid2,integer1d,(/k,1/)*1_mpi_offset_kind,(/1,ncellsf/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)
print*,'done cellsoncell'      

   rc = nf90_inq_varid(ncid,'cellsOnEdge',varid)
   rc = nf90mpi_def_var(ncidf,'cellsOnEdge',nf90mpi_INT, (/dimid2,dimidef/), varid2)
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,2
      rc = nf90_get_var  (ncid,varid,integer1d_edge,(/k,1/),(/1,nedgesf/))
      rc = nf90mpi_put_var_all(ncidf,varid2,integer1d_edge,(/k,1/)*1_mpi_offset_kind,(/1,nedgesf/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)
print*,'done cellsonedge'      


   rc = nf90_inq_varid(ncid,'weightsOnEdge',varid)
   rc = nf90mpi_def_var(ncidf,'weightsOnEdge',NF90mpi_DOUBLE, (/dimidme2,dimidef/), varid2)
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,12
      rc = nf90_get_var  (ncid,varid,real8_edge,(/k,1/),(/1,nedgesf/))
      rc = nf90mpi_put_var_all(ncidf,varid2,real8_edge,(/k,1/)*1_mpi_offset_kind,(/1,nedgesf/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)

print*,'done weightsonedge'      

   rc = nf90_inq_varid(ncid,'edgesOnEdge',varid)
   rc = nf90mpi_def_var(ncidf,'edgesOnEdge',nf90mpi_INT, (/dimidme2,dimidef/), varid2)
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,12
      rc = nf90_get_var  (ncid,varid,integer1d_edge,(/k,1/),(/1,nedgesf/))
      rc = nf90mpi_put_var_all(ncidf,varid2,integer1d_edge,(/k,1/)*1_mpi_offset_kind,(/1,nedgesf/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)

print*,'done edgesonedge'      

   rc = nf90_inq_varid(ncid,'verticesOnEdge',varid)
   rc = nf90mpi_def_var(ncidf,'verticesOnEdge',nf90mpi_INT, (/dimid2,dimidef/), varid2)
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,2
      rc = nf90_get_var  (ncid,varid,integer1d_edge,(/k,1/),(/1,nedgesf/))
      rc = nf90mpi_put_var_all(ncidf,varid2,integer1d_edge,(/k,1/)*1_mpi_offset_kind,(/1,nedgesf/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)
print*,'done verticesonedge'      

   rc = nf90_inq_varid(ncid,'kiteAreasOnVertex',varid)
   rc = nf90mpi_def_var(ncidf,'kiteAreasOnVertex',NF90mpi_DOUBLE, (/dimidvd,dimidvf/), varid2)
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,3
      rc = nf90_get_var  (ncid,varid,real8_vert,(/k,1/),(/1,nverticesf/))
      real8_vert = real8_vert * sphere_radius ** 2
      rc = nf90mpi_put_var_all(ncidf,varid2,real8_vert,(/k,1/)*1_mpi_offset_kind,(/1,nverticesf/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'edgesOnVertex',varid)
   rc = nf90mpi_def_var(ncidf,'edgesOnVertex',nf90mpi_INT, (/dimidvd,dimidvf/), varid2)
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,3
      rc = nf90_get_var  (ncid,varid,integer1d_vert,(/k,1/),(/1,nverticesf/))
      rc = nf90mpi_put_var_all(ncidf,varid2,integer1d_vert,(/k,1/)*1_mpi_offset_kind,(/1,nverticesf/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'cellsOnVertex',varid)
   rc = nf90mpi_def_var(ncidf,'cellsOnVertex',nf90mpi_INT, (/dimidvd,dimidvf/), varid2)
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,3
      rc = nf90_get_var  (ncid,varid,integer1d_vert,(/k,1/),(/1,nverticesf/))
      rc = nf90mpi_put_var_all(ncidf,varid2,integer1d_vert,(/k,1/)*1_mpi_offset_kind,(/1,nverticesf/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'edgesOnCell',varid)
   rc = nf90mpi_def_var(ncidf,'edgesOnCell',nf90mpi_INT, (/dimidme,dimidcf/), varid2)
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,6
      rc = nf90_get_var  (ncid,varid,integer1d,(/k,1/),(/1,ncellsf/))
      rc = nf90mpi_put_var_all(ncidf,varid2,integer1d,(/k,1/)*1_mpi_offset_kind,(/1,ncellsf/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)

print*,'done edgesoncell'      


   rc = nf90_inq_varid(ncid,'verticesOnCell',varid)
   rc = nf90mpi_def_var(ncidf,'verticesOnCell',nf90mpi_INT, (/dimidme,dimidcf/), varid2)
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,6
      rc = nf90_get_var  (ncid,varid,integer1d,(/k,1/),(/1,ncellsf/))
      rc = nf90mpi_put_var_all(ncidf,varid2,integer1d,(/k,1/)*1_mpi_offset_kind,(/1,ncellsf/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)

print*,'mesh variables copied'

   rc = nf90_inq_varid(ncid,'fCell',varid)
   rc = nf90_get_var  (ncid,varid,real8_1d)
   rc = nf90mpi_def_var(ncidf,'fCell',nf90mpi_DOUBLE, (/dimidcf/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,real8_1d)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'fEdge',varid)
   rc = nf90_get_var  (ncid,varid,real8_edge)
   rc = nf90mpi_def_var(ncidf,'fEdge',nf90mpi_DOUBLE, (/dimidef/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,real8_edge)
   rc = nf90mpi_redef(ncidf)

   rc = nf90_inq_varid(ncid,'fVertex',varid)
   rc = nf90_get_var  (ncid,varid,real8_vert)
   rc = nf90mpi_def_var(ncidf,'fVertex',nf90mpi_DOUBLE, (/dimidvf/), varid)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid,real8_vert)
   rc = nf90mpi_redef(ncidf)
print*,'coriolis copied'

print*,'start initialize_interpolate_variable'
      ALLOCATE(remap_indices(3,ncellsf))  !
      ALLOCATE(remap_wgts(3,ncellsf))     !
      ALLOCATE(remap_indicesv(4,nverticesf))  !
      ALLOCATE(remap_wgtsv(4,nverticesf))     !
   rc = nf90_open(trim(terpfile),nf90_nowrite,nctem)
      rc = nf90_inq_varid(nctem,'remap_indices',varid)
print*,'remap_indices inqvarid ',rc
      rc = nf90_get_var(nctem,varid,remap_indices,(/1,1,1/),(/3,ncellsf,1/))
print*,'remap_indices getvar ',rc
      rc = nf90_inq_varid(nctem,'remap_wgts',varid)
print*,'remap_wgts inqvarid ',rc
      rc = nf90_get_var(nctem,varid,remap_wgts,(/1,1,1/),(/3,ncellsf,1/))
print*,'remap_wgts getvar ',rc
      rc = nf90_inq_varid(nctem,'remap_indicesv',varid)
      rc = nf90_get_var(nctem,varid,remap_indicesv)
      rc = nf90_inq_varid(nctem,'remap_wgtsv',varid)
      rc = nf90_get_var(nctem,varid,remap_wgtsv)
   rc = nf90_close(nctem)
print*,'done initialize_interpolate_variable'
   
print*,'start ice variables'

! ice restart variables to be interpolated
   !char xtime(Time, StrLen) ;
   !METHOD: copy
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'xtime',  &
                            pdimids=(/dimidsl,dimidt/), pxtype=nf90mpi_char  )
   rc = nf90_inq_varid(ncidc,'xtime',varid)
   rc = nf90_get_var  (ncidc,varid,time)
   !time(1:4)='0001'
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid2,time, &
        start=(/1_mpi_offset_kind,1_mpi_offset_kind/),count=(/64_mpi_offset_kind,1_mpi_offset_kind/))
   rc = nf90mpi_redef(ncidf)

   !char simulationStartTime(StrLen) ;
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'simulationStartTime', &
                            pdimids=(/dimidsl/), pxtype=nf90mpi_char  )
   rc = nf90_inq_varid(ncidc,'simulationStartTime',varid)
   rc = nf90_get_var  (ncidc,varid,time)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid2,time, &
        start=(/1_mpi_offset_kind,1_mpi_offset_kind/),count=(/64_mpi_offset_kind,1_mpi_offset_kind/))
   rc = nf90mpi_redef(ncidf)

   !double iceAreaCategory(Time, nCells, nCategories, ONE) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'iceAreaCategory', debug=.true., &
                            pdimids=(/dimid1, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  ) 
   print*,'define iceareacategory done'
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'iceAreaCategory',varid)
   !DO k = 1,5
   print*,'begin interpolate_and_write_icearea'
      CALL interpolate_and_write_icearea(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1,1,1/),         &
                                          countin=(/1,1,ncellsc,1/), &
                                          countout=(/1,5,ncellsf,1/), debug=.true. )
   !ENDDO
   rc = nf90mpi_redef(ncidf)
print*,'iceareacategory done'
        
   !double iceVolumeCategory(Time, nCells, nCategories, ONE) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'iceVolumeCategory', &
                            pdimids=(/dimid1, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'iceVolumeCategory',varid)
   DO k = 1,5
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,k,1,1/),         &
                                          countin=(/1,1,ncellsc,1/), &
                                          countout=(/1,1,ncellsf,1/) )
   ENDDO
   rc = nf90mpi_redef(ncidf)
        
   !double snowVolumeCategory(Time, nCells, nCategories, ONE) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'snowVolumeCategory', &
                            pdimids=(/dimid1, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'snowVolumeCategory',varid)
   DO k = 1,5
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,k,1,1/),         &
                                          countin=(/1,1,ncellsc,1/), &
                                          countout=(/1,1,ncellsf,1/) )
   ENDDO
   rc = nf90mpi_redef(ncidf)
        
   !double surfaceTemperature(Time, nCells, nCategories, ONE) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'surfaceTemperature', &
                            pdimids=(/dimid1, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'surfaceTemperature',varid)
   DO k = 1,5
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,k,1,1/),         &
                                          countin=(/1,1,ncellsc,1/), &
                                          countout=(/1,1,ncellsf,1/) )
   ENDDO
   rc = nf90mpi_redef(ncidf)
        
   !double iceEnthalpy(Time, nCells, nCategories, nIceLayers) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'iceEnthalpy', &
                            pdimids=(/dimidily, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'iceEnthalpy',varid)
   DO l = 1,7
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                              ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/), zlim = 1.0_8 )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)
print*,'iceenthalpy done'
        
   !double iceSalinity(Time, nCells, nCategories, nIceLayers) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'iceSalinity', &
                            pdimids=(/dimidily, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'iceSalinity',varid)
   DO l = 1,7
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                              ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)

   !double snowEnthalpy(Time, nCells, nCategories, nSnowLayers) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'snowEnthalpy',  &
                            pdimids=(/dimidsly, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'snowEnthalpy',varid)
   DO l = 1,5
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/), zlim = 1.0_8 )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)

   !double snowIceMass(Time, nCells, nCategories, nSnowLayers) ;
   !        snowIceMass:units = "kg/m^3" ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'snowIceMass', &
                            units="kg/m^3", &
                            pdimids=(/dimidsly, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'snowIceMass',varid)
   DO l = 1,5
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)
print*,'snowicemass done'

   !double snowLiquidMass(Time, nCells, nCategories, nSnowLayers) ;
   !        snowLiquidMass:units = "kg/m^3" ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'snowLiquidMass', &
                            units="kg/m^3", &
                            pdimids=(/dimidsly, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'snowLiquidMass',varid)
   DO l = 1,5
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)
   
   !double snowGrainRadius(Time, nCells, nCategories, nSnowLayers) ;
   !        snowGrainRadius:units = "um" ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'snowGrainRadius', &
                            units="um", &
                            pdimids=(/dimidsly, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'snowGrainRadius',varid)
   DO l = 1,5
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)
   
   !double snowDensity(Time, nCells, nCategories, nSnowLayers) ;
   !        snowDensity:units = "kg/m^3" ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'snowDensity', &
                            units="kg/m^3", &
                            pdimids=(/dimidsly, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'snowDensity',varid)
   DO l = 1,5
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/), zlim = 1.0_8 )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)
print*,'snowdensity done'
   
   !double iceAge(Time, nCells, nCategories, ONE) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'iceAge', &
                            pdimids=(/dimid1, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'iceAge',varid)
   DO l = 1,1
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)

   !double firstYearIceArea(Time, nCells, nCategories, ONE) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'firstYearIceArea', &
                            pdimids=(/dimid1, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'firstYearIceArea',varid)
   DO l = 1,1
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)

   !double levelIceArea(Time, nCells, nCategories, ONE) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'levelIceArea', &
                            pdimids=(/dimid1, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'levelIceArea',varid)
   DO l = 1,1
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)
print*,'levelicearea done'

   !double levelIceVolume(Time, nCells, nCategories, ONE) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'levelIceVolume', &
                            pdimids=(/dimid1, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'levelIceVolume',varid)
   DO l = 1,1
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)

   !double pondArea(Time, nCells, nCategories, ONE) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'pondArea', &
                            pdimids=(/dimid1, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'pondArea',varid)
   DO l = 1,1
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)

   !double pondDepth(Time, nCells, nCategories, ONE) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'pondDepth' , &
                            pdimids=(/dimid1, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'pondDepth',varid)
   DO l = 1,1
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)

   !double pondLidThickness(Time, nCells, nCategories, ONE) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'pondLidThickness', &
                            pdimids=(/dimid1, dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'pondLidThickness',varid)
   DO l = 1,1
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/l,k,1,1/),         &
                                             countin=(/1,1,ncellsc,1/), &
                                             countout=(/1,1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)
print*,'pondlidthickness done'

   !double uVelocity(Time, nVertices) ;
   ! METHOD: set to zero for now and test spinup
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'uVelocity', &
                            pdimids=(/dimidvf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   real8_vert = 0.0
   rc = nf90mpi_put_var_all(ncidf,varid2,real8_vert,start=(/1,1/)*1_mpi_offset_kind,count=(/nverticesf,1/)*1_mpi_offset_kind)
   rc = nf90mpi_redef(ncidf)

   !double vVelocity(Time, nVertices) ;
   ! METHOD: set to zero for now and test spinup
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'vVelocity', &
                            pdimids=(/dimidvf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   real8_vert = 0.0
   rc = nf90mpi_put_var_all(ncidf,varid2,real8_vert,start=(/1,1/)*1_mpi_offset_kind,count=(/nverticesf,1/)*1_mpi_offset_kind)
   rc = nf90mpi_redef(ncidf)

   !double stress11var(Time, nCells, maxEdges) ;
   ! METHOD: set to zero for now and test spinup
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'stress11var', &
                            pdimids=(/dimidme, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   real8_1d = 0.0
   do k = 1,6
     rc = nf90mpi_put_var_all(ncidf,varid2,real8_vert,start=(/k,1,1/)*1_mpi_offset_kind,count=(/1,ncellsf,1/)*1_mpi_offset_kind)
   enddo
   rc = nf90mpi_redef(ncidf)

   !double stress22var(Time, nCells, maxEdges) ;
   ! METHOD: set to zero for now and test spinup
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'stress22var', &
                            pdimids=(/dimidme, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   real8_1d = 0.0
   do k = 1,6
      rc = nf90mpi_put_var_all(ncidf,varid2,real8_vert,start=(/k,1,1/)*1_mpi_offset_kind,count=(/1,ncellsf,1/)*1_mpi_offset_kind)
   enddo
   rc = nf90mpi_redef(ncidf)

   !double stress12var(Time, nCells, maxEdges) ;
   ! METHOD: set to zero for now and test spinup
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'stress12var', &
                            pdimids=(/dimidme, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   real8_1d = 0.0
   do k = 1,6
      rc = nf90mpi_put_var_all(ncidf,varid2,real8_vert,start=(/k,1,1/)*1_mpi_offset_kind,count=(/1,ncellsf,1/)*1_mpi_offset_kind)
   enddo
   rc = nf90mpi_redef(ncidf)

   !int solveVelocityPrevious(Time, nVertices) ;
   ! METHOD: set to 0 for now and test spinup
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'solveVelocityPrevious', &
                            pdimids=(/dimidvf, dimidt/), pxtype=nf90mpi_int  )
   rc = nf90mpi_enddef(ncidf)
   integer1d_vert = 0.0
   rc = nf90mpi_put_var_all(ncidf,varid2,real8_vert,start=(/1,1/)*1_mpi_offset_kind,count=(/nverticesf,1/)*1_mpi_offset_kind)
   rc = nf90mpi_redef(ncidf)

   !double freezeOnset(Time, nCells) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'freezeOnset', &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'freezeOnset',varid)
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)
print*,'freezeonset done'

   !double snowfallRate(Time, nCells) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'snowfallRate', &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'snowfallRate',varid)
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double pondSnowDepthDifference(Time, nCells, nCategories) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'pondSnowDepthDifference', &
                            pdimids=(/dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'pondSnowDepthDifference',varid)
   DO l = 1,1
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/k,1,1/),         &
                                             countin=(/1,ncellsc,1/), &
                                             countout=(/1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)

   !double pondLidMeltFluxFraction(Time, nCells, nCategories) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'pondLidMeltFluxFraction', &
                            pdimids=(/dimidcat, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'pondLidMeltFluxFraction',varid)
   DO l = 1,1
      DO k = 1,5
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/k,1,1/),         &
                                             countin=(/1,ncellsc,1/), &
                                             countout=(/1,ncellsf,1/) )
      ENDDO
   ENDDO
   rc = nf90mpi_redef(ncidf)

   !double solarZenithAngleCosine(Time, nCells) ;
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'solarZenithAngleCosine', &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'solarZenithAngleCosine',varid)
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double shortwaveScalingFactor(Time, nCells) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'shortwaveScalingFactor', &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'shortwaveScalingFactor',varid)
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)
print*,'shortwaveScalingFactor done'

   !double shortwaveVisibleDirectDown(Time, nCells) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'shortwaveVisibleDirectDown', &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'shortwaveVisibleDirectDown',varid)
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double shortwaveVisibleDiffuseDown(Time, nCells) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'shortwaveVisibleDiffuseDown', &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'shortwaveVisibleDiffuseDown',varid)
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double shortwaveIRDirectDown(Time, nCells) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'shortwaveIRDirectDown', &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'shortwaveIRDirectDown',varid)
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double shortwaveIRDiffuseDown(Time, nCells) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'shortwaveIRDiffuseDown', &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'shortwaveIRDiffuseDown',varid)
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)
print*,'shortwaveIRDiffuseDown done'

   !double oceanStressCellU(Time, nCells) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'oceanStressCellU', &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'oceanStressCellU',varid)
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double oceanStressCellV(Time, nCells) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'oceanStressCellV', &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'oceanStressCellV',varid)
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double seaSurfaceTemperature(Time, nCells) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'seaSurfaceTemperature', &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double , debug=.true.  )
print*,'varid,varid2',varid,varid2
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'seaSurfaceTemperature',varid)
print*,'varid,varid2',varid,varid2
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) , debug=.true. )
   rc = nf90mpi_redef(ncidf)

   !double freezingMeltingPotential(Time, nCells) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'freezingMeltingPotential', &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double , debug=.true.  )
print*,'varid,varid2',varid,varid2
   rc = nf90mpi_enddef(ncidf)
print*,'varid,varid2',varid,varid2
   rc = nf90_inq_varid(ncidc,'freezingMeltingPotential',varid)
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) , debug=.true.)
   rc = nf90mpi_redef(ncidf)
print*,'freezingMeltingPotential done'

   !double airOceanDragCoefficientRatio(Time, nCells) ;
   !METHOD: interpolate 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'airOceanDragCoefficientRatio',  &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double , debug=.true.  )
print*,'varid,varid2',varid,varid2
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid(ncidc,'airOceanDragCoefficientRatio',varid)
print*,'varid,varid2',varid,varid2
         CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                             ncidf,           &
                                             varid2,        &
                                             start=(/1,1/),         &
                                             countin=(/ncellsc,1/), &
                                             countout=(/ncellsf,1/) , debug=.true.)
   rc = nf90mpi_redef(ncidf)
print*,'airOceanDragCoefficientRatio done'

   !int newlyFormedIce(nCells, nCategories) ;
   !METHOD: assign constant 0 
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'newlyFormedIce', &
                            pdimids=(/dimidcat,dimidcf/), pxtype=nf90mpi_int  )
print*,'newlyFormedIce defined'
   rc = nf90mpi_enddef(ncidf)
print*,'newlyFormedIce alloc'
   integer1d = 0
print*,'newlyFormedIce assigned'
   do k = 1,5
      rc = nf90mpi_put_var_all(ncidf,varid2,integer1d,(/k,1/)*1_mpi_offset_kind,(/1,ncellsf/)*1_mpi_offset_kind)
   enddo
   rc = nf90mpi_redef(ncidf)
print*,'done ice variables'


! close up the output file
   rc = nf90mpi_close(ncidf)
   
          call MPI_finalize(err)

   contains
      
 
!------------------------------------------

      SUBROUTINE define_rst_variable(varid, pvarid,                   &
                                     ncidin, pncidout,       &
                                     varname,          &
                                     units, long_name ,               &
                                     pdimids, pxtype, debug)
      ! input arguments
      INTEGER, INTENT(IN) :: ncidin,   & ! input file id
                             pncidout    ! pnetcdf output file id
      INTEGER, INTENT(IN), DIMENSION(:) :: pdimids    ! pnetcdf output file dimension ids
      INTEGER, INTENT(IN), OPTIONAL :: pxtype !  pnetcdf variables types
      CHARACTER(LEN=*), INTENT(IN) :: varname ! variable name
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: units ! units attribute
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: long_name ! long_name attribute
      LOGICAL, INTENT(IN), OPTIONAL :: debug
      
      ! output arguments
      INTEGER, INTENT(OUT) :: varid, pvarid ! netcdf and pnetcdf output file variable ids
      
      ! local variables
      INTEGER :: rc, vid
      
      rc = nf90_inq_varid(ncidin ,TRIM(varname),varid  )
      rc = nf90mpi_def_var(pncidout,TRIM(varname),pxtype, pdimids, pvarid)
      if(PRESENT(debug))PRINT*,'rc: pnet define_rst_variable def_var',trim(varname),rc

      IF(PRESENT(long_name)) THEN
         rc = nf90mpi_put_att(pncidout, pvarid, 'long_name', trim(long_name))
         if(PRESENT(debug))PRINT*,'rc: pnet define_rst_variable put_att1',trim(varname),rc
      ENDIF
                                     
      IF(PRESENT(units)) THEN
         rc = nf90mpi_put_att(pncidout, pvarid, 'units', trim(units))
         if(PRESENT(debug))PRINT*,'rc: pnet define_rst_variable put_att2',trim(varname),rc
      ENDIF
                                     
      END SUBROUTINE define_rst_variable

!------------------------------------------

      SUBROUTINE interpolate_and_write_variable(       &
                             lev,                      &
                             ncidr, varid,             &
                             pncid,              &
                             pvarid,            &
                             start, countin, countout, &
                             spval, debug, zlim )
      ! This routine will perform a horizontal remapping from the coarser to finer grid.
      ! It will read in and write out one horizontal slice of the variable.
      
      ! Argument list
      INTEGER, INTENT(IN) :: lev  ! the vertical level of interpolation.
                                  ! call this routine with lev=k for variables with an oceaninc vertical dimension
                                  !   otherwise call with lev=1
      INTEGER, INTENT(IN) :: ncidr, varid  ! file and variable id for the input coarse grid file
      INTEGER, INTENT(IN) :: pncid, pvarid ! file and variable id for the output fine grid pnetcdf file
      INTEGER, INTENT(IN), DIMENSION(:) :: start, countin, countout
      REAL*8, INTENT(IN), OPTIONAL :: spval
      REAL*8, INTENT(IN), OPTIONAL :: zlim
      LOGICAL, INTENT(IN), OPTIONAL :: debug
      
      ! local variables
      INTEGER :: n, l  ! loop indices
      INTEGER :: rc

      REAL*8 :: spval_tem, zlim_tem
      real*8 :: wgtsum
      
      if(lev.ne.1)stop 'unavailable lev'

      rc = nf90_get_var(ncidr, varid, real8_1d_in, start, countin)
      if(PRESENT(debug))PRINT*,'rc: interpolate_and_write_variable get_var',rc
!$OMP parallel num_threads(1)
      if(present(spval)) then
         !$OMP do
         DO n = 1,ncellsf
            real8_1d(n) = spval
         enddo
         spval_tem = spval
      else
         !$OMP do
         DO n = 1,ncellsf
            real8_1d(n) = 0.0
         enddo
         spval_tem = 0.0
      endif
      if(present(zlim)) then
         zlim_tem = zlim
      else
         zlim_tem = -1.0
      endif
      !print*,'lev, spval',lev,spval_tem
      !print*,'start, ',start
      !print*,'countin, ',countin
      !print*,'indices',remap_indices(:,1563,lev)
      !print*,'weights',remap_wgts(:,1563,lev)
      !print*,'input  ',real8_1d_in(remap_indices(1,1563,lev)),real8_1d_in(remap_indices(2,1563,lev)),real8_1d_in(remap_indices(3,1563,lev))
      !$OMP do private(l, wgtsum)
      DO n = 1,ncellsf
         wgtsum = 0.0
         DO l = 1,3
            IF(remap_indices(l,n) > 0) then
               if(real8_1d(n) == spval_tem) real8_1d(n) = 0.0
               if( abs(real8_1d_in(remap_indices(l,n))) >= zlim_tem ) then
                  real8_1d(n) = real8_1d(n) +     &
                                real8_1d_in(remap_indices(l,n)) * remap_wgts(l,n)
                  wgtsum = wgtsum + remap_wgts(l,n)
               endif
            endif
         ENDDO
         if(wgtsum > 0.0)   &
            real8_1d(n) = real8_1d(n) / wgtsum
         if(PRESENT(debug) .and. n < 10) then
         print*,'n',real8_1d(n) ,wgtsum
         print*,'input',real8_1d_in(remap_indices(1,n)),real8_1d_in(remap_indices(2,n)),real8_1d_in(remap_indices(3,n))
         endif
      ENDDO
!$OMP end parallel

      !print*,'output ',real8_1d(1563)
      rc = nf90mpi_put_var_all(pncid,pvarid,real8_1d, start*1_mpi_offset_kind, countout*1_mpi_offset_kind)
      if(PRESENT(debug))PRINT*,'rc: interpolate_and_write_variable nfmpi_put_var',rc
                             
      END SUBROUTINE interpolate_and_write_variable

!------------------------------------------

      SUBROUTINE interpolate_and_write_icearea(       &
                             lev,                      &
                             ncidr, varid,             &
                             pncid,              &
                             pvarid,            &
                             start, countin, countout, &
                             spval, debug )
      ! This routine will perform a horizontal remapping from the coarser to finer grid.
      ! It will read in and write out one horizontal slice of the variable.
      
      ! Argument list
      INTEGER, INTENT(IN) :: lev  ! the vertical level of interpolation.
                                  ! call this routine with lev=k for variables with an oceaninc vertical dimension
                                  !   otherwise call with lev=1
      INTEGER, INTENT(IN) :: ncidr, varid  ! file and variable id for the input coarse grid file
      INTEGER, INTENT(IN) :: pncid, pvarid ! file and variable id for the output fine grid pnetcdf file
      INTEGER, INTENT(IN), DIMENSION(:) :: start, countin, countout
      REAL*8, INTENT(IN), OPTIONAL :: spval
      LOGICAL, INTENT(IN), OPTIONAL :: debug
      
      ! local variables
      INTEGER :: n, l, k  ! loop indices
      INTEGER :: rc

      REAL*8 :: spval_tem
      real*8 :: icearea(5,ncellsf)

      if(lev.ne.1)stop 'unavailable icearea lev'
      
      do k = 1,5
      rc = nf90_get_var(ncidr, varid, real8_1d_in, (/1,k,1,1/), countin)
      print*,'1 min icearea k ',k,minval(real8_1d_in),minloc(real8_1d_in)
      print*,'1 max icearea k ',k,maxval(real8_1d_in),maxloc(real8_1d_in)
      if(PRESENT(debug))PRINT*,'rc: interpolate_and_write_variable get_var',rc
!$OMP parallel  
      if(present(spval)) then
!$OMP  do 
         DO n = 1,ncellsf
            icearea(k,n) = spval
         enddo
         spval_tem = spval
      else
!$OMP  do 
         DO n = 1,ncellsf
            icearea(k,n) = 0.0
         enddo
         spval_tem = 0.0
      endif
      !print*,'1 min icearea k ',k,minval(icearea(k,:))
      !print*,'lev, spval',lev,spval_tem
      !print*,'start, ',start
      !print*,'countin, ',countin
      !print*,'indices',remap_indices(:,1563,lev)
      !print*,'weights',remap_wgts(:,1563,lev)
      !print*,'input  ',real8_1d_in(remap_indices(1,1563,lev)),real8_1d_in(remap_indices(2,1563,lev)),real8_1d_in(remap_indices(3,1563,lev))

!$OMP  do private(l)
      DO n = 1,ncellsf
         DO l = 1,3
            IF(remap_indices(l,n) > 0) then
               if(icearea(k,n) == spval_tem) icearea(k,n) = 0.0
               icearea(k,n) = icearea(k,n) +     &
                             real8_1d_in(remap_indices(l,n)) * remap_wgts(l,n)
            endif
         !if(n==25863 .and. k==5) then
         !   print*,'l,index,weight',l,remap_indices(l,n,lev),remap_wgts(l,n,lev)
         !   print*,'input',real8_1d_in(remap_indices(l,n,lev))
         !   print*,'output',icearea(k,n)
         !endif
         ENDDO
         
         ! keep the sum < 1.0
         icearea(k,n) = icearea(k,n) * 0.999999
      
      ENDDO
!$OMP end parallel  
      
      print*,'2 min icearea k ',k,minval(icearea(k,:)),minloc(icearea(k,:))
      print*,'2 max icearea k ',k,maxval(icearea(k,:)),maxloc(icearea(k,:))
      enddo
      
      
      !print*,'output ',real8_1d(1563)
      rc = nf90mpi_put_var_all(pncid,pvarid,icearea, start*1_mpi_offset_kind, countout*1_mpi_offset_kind)
      if(PRESENT(debug))PRINT*,'rc: pnet interpolate_and_write_variable put_var',rc
                             
      END SUBROUTINE interpolate_and_write_icearea

!------------------------------------------
   
end program build_ice_restart_file
