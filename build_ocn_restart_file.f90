program build_ocn_restart_file

! compile with netcdf, pnetcdf, and, optionally openmp

! setenv OMP_NUM_THREADS 16

use netcdf
use pnetcdf, only : nf90mpi_create, nf90mpi_close, nf90mpi_enddef, nf90mpi_def_dim, &
                    nf90mpi_put_att, nf90mpi_def_var, nf90mpi_put_var_all, nf90mpi_inq_varid, &
                    nf90mpi_put_var, nf90mpi_redef, nf90mpi_open,nf90mpi_inq_dimid,  &
                    nf90mpi_inquire_dimension, nf90mpi_get_var_all,   &
                    nf90mpi_clobber => NF90_CLOBBER, NF90mpi_64BIT_DATA=>NF90_64BIT_DATA, &
                    nf90mpi_unlimited => nf90_unlimited, NF90mpi_GLOBAL=>NF90_GLOBAL, &
                    NF90mpi_DOUBLE=>NF90_DOUBLE, NF90mpi_INT=>NF90_INT, nf90mpi_char=>nf90_char, &
                    nf90mpi_nowrite=>nf90_nowrite, nf90mpi_sync
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
              dimid1, dimidcat, dimidily, dimidsly, varid2, dimidl, &
              dimidv1, dimidfg
   integer :: natts, maxlev, maxiter
   character(len=80) :: attname

   ! grid variables
   ! vertices
   integer :: ncellsc, ncellsf, maxlinks, nverticesc, nverticesf, nedgesc, nedgesf
   
   
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
   real*4, dimension(:,:), allocatable :: remap_wgts, remap_wgts3d
   integer, dimension(:,:), allocatable :: remap_indices, remap_indices3d
   real*4, dimension(:,:), allocatable :: remap_wgtsv
   integer, dimension(:,:), allocatable :: remap_indicesv

   
   ! mpi variables
   integer :: err, rank, nprocs, info, nctem
           integer(kind=MPI_OFFSET_KIND) G_NY, myOff, block_start, &
                                        global_nx, global_ny
          integer(kind=MPI_OFFSET_KIND) start(2), count(2), dimtem
          integer(kind=MPI_OFFSET_KIND) malloc_size, sum_size
  
        integer, dimension(:), allocatable :: intwork1d
        integer, dimension(:,:), allocatable :: intwork2d
        real*8, dimension(:), allocatable :: r8work1d, bottomdepth_in, ssh_in, areacell, wrkarr, restingthickness
        real*8, dimension(:,:), allocatable :: r8work2d, layerthickness_in

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

   !choose the case and define the case-specific variables
   print*,'enter the ocean mesh file name'
   read(5,'(a128)')infile
   print*,'enter the interpolation weights file name'
   read(5,'(a128)')terpfile
   print*,'enter the input ocean restart file name'
   read(5,'(a128)')coarsefile
   print*,'enter the ocean restart file name to be created'
   read(5,'(a128)')finefile
   
   maxlev = 100
   rc = nf90mpi_create(mpi_comm_world,trim(finefile),IOR(NF90_CLOBBER, NF90_64BIT_DATA),info,ncidf)
   print*,'create fine restart',rc
   if(rc .ne. 0)print*,trim(finefile)

          call MPI_Info_free(info, err)

   ! open the necessary netcdf files
   rc = nf90_open(trim(infile),nf90_nowrite,ncid)
   print*,'open input meshfile',rc
   if(rc .ne. 0)print*,trim(infile)
   rc = nf90mpi_open(mpi_comm_world,trim(coarsefile),nf90mpi_nowrite,info,ncidc)
   print*,'open coarse restart',rc
   if(rc .ne. 0)print*,trim(coarsefile)
      ! read the coarse grid dimensions
      rc = nf90mpi_inq_dimid(ncidc,'nCells',dimid)
      rc = nf90mpi_inquire_dimension(ncidc,dimid,len=dimtem)
      ncellsc = dimtem
      rc = nf90mpi_inq_dimid(ncidc,'nVertices',dimid)
      rc = nf90mpi_inquire_dimension(ncidc,dimid, len=dimtem)
      nverticesc = dimtem
      rc = nf90mpi_inq_dimid(ncidc,'nEdges',dimid)
      rc = nf90mpi_inquire_dimension(ncidc,dimid, len=dimtem)
      nedgesc = dimtem

   allocate(real8_1d_in(ncellsc))

!!!!!! define and copy the dimensions   
   rc = nf90_inq_dimid(ncid,'nCells',dimid)
   rc = nf90_inquire_dimension(ncid,dimid,len=ncellsf)
   !rc = nf90_inq_dimid(ncidc,'nCells',dimid)
   !rc = nf90_inquire_dimension(ncidc,dimid,len=ncellsc)
   rc = nf90mpi_def_dim(ncidf,'nCells',ncellsf*1_mpi_offset_kind,dimidcf)

   rc = nf90_inq_dimid(ncid,'nEdges',dimid)
   rc = nf90_inquire_dimension(ncid,dimid,len=nedgesf)
   !rc = nf90_inq_dimid(ncidc,'nEdges',dimid)
   !rc = nf90_inquire_dimension(ncidc,dimid,len=nedgesc)
   rc = nf90mpi_def_dim(ncidf,'nEdges',nedgesf*1_mpi_offset_kind,dimidef)

   rc = nf90_inq_dimid(ncid,'nVertices',dimid)
   rc = nf90_inquire_dimension(ncid,dimid,len=nverticesf)
   !rc = nf90_inq_dimid(ncidc,'nVertices',dimid)
   !rc = nf90_inquire_dimension(ncidc,dimid,len=nverticesc)
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
   rc = nf90mpi_def_dim(ncidf,'nVertLevels',100*1_mpi_offset_kind,dimidl)
   rc = nf90mpi_def_dim(ncidf,'nVertLevelsP1',101*1_mpi_offset_kind,dimidv1)
   rc = nf90mpi_def_dim(ncidf,'nForcingGroupsMax',4*1_mpi_offset_kind,dimidfg)

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
   rc = nf90mpi_def_var(ncidf,'weightsOnEdge',nf90mpi_double, (/dimidme2,dimidef/), varid2)
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,12
      rc = nf90_get_var  (ncid,varid,real8_edge,(/k,1/),(/1,nedgesf/))
      rc = nf90mpi_put_var_all(ncidf,varid2,real8_edge,(/k,1/)*1_mpi_offset_kind,(/1,nedgesf/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)

print*,'done weightsonedge'      

   rc = nf90_inq_varid(ncid,'edgesOnEdge',varid)
   rc = nf90mpi_def_var(ncidf,'edgesOnEdge',nf90mpi_int, (/dimidme2,dimidef/), varid2)
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
   rc = nf90mpi_def_var(ncidf,'kiteAreasOnVertex',nf90mpi_double, (/dimidvd,dimidvf/), varid2)
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

   rc = nf90_inq_varid(ncid,'cellsOnCell',varid)
   rc = nf90mpi_def_var(ncidf,'cellsOnCell',nf90mpi_INT, (/dimidme,dimidcf/), varid2)
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,6
      rc = nf90_get_var  (ncid,varid,integer1d,(/k,1/),(/1,ncellsf/))
      rc = nf90mpi_put_var_all(ncidf,varid2,integer1d,(/k,1/)*1_mpi_offset_kind,(/1,ncellsf/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)
print*,'done cellsoncell'      

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
      ALLOCATE(remap_indices3d(3,ncellsf))  !
      ALLOCATE(remap_wgts3d(3,ncellsf))     !
      ALLOCATE(remap_indicesv(4,nverticesf))  !
      ALLOCATE(remap_wgtsv(4,nverticesf))     !
   rc = nf90_open(trim(terpfile),nf90_nowrite,nctem)
      rc = nf90_inq_varid(nctem,'remap_indices',varid)
      rc = nf90_get_var(nctem,varid,remap_indices,(/1,1,1/),(/3,ncellsf,1/))
      rc = nf90_inq_varid(nctem,'remap_wgts',varid)
      rc = nf90_get_var(nctem,varid,remap_wgts,(/1,1,1/),(/3,ncellsf,1/))
      rc = nf90_inq_varid(nctem,'remap_indicesv',varid)
      rc = nf90_get_var(nctem,varid,remap_indicesv)
      rc = nf90_inq_varid(nctem,'remap_wgtsv',varid)
      rc = nf90_get_var(nctem,varid,remap_wgtsv)
   !rc = nf90_close(nctem)
print*,'done initialize_interpolate_variable'
   
print*,'start ocn variables'

! ice restart variables to be interpolated
   !char xtime(Time, StrLen) ;
   !METHOD: copy
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'xtime',  &
                            pdimids=(/dimidsl,dimidt/), pxtype=nf90mpi_char  )
   rc = nf90mpi_get_var_all  (ncidc,varid,time)
   print*,'time ',rc, trim(time)
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
   rc = nf90mpi_get_var_all  (ncidc,varid,time)
   print*,'simstarttime ',rc, trim(time)
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid2,time, &
        start=(/1_mpi_offset_kind,1_mpi_offset_kind/),count=(/64_mpi_offset_kind,1_mpi_offset_kind/))
   rc = nf90mpi_redef(ncidf)

   !double refBottomDepth(nVertLevels) ;
   !        refBottomDepth:units = "m" ;
   !        refBottomDepth:long_name = "Reference depth of ocean for each vertical level. Used in \'z-level\' type runs." ;
   !METHOD: copy 
   print*,'refbottomdepth'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'refBottomDepth', &
                            units='m',                    &
                            long_name="Reference depth of ocean for each vertical level. Used in \'z-level\' type runs.", &
                            pdimids=(/dimidl/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_get_var_all  (ncidc,varid,real8_1d(1:100))
   rc = nf90mpi_put_var_all(ncidf,varid2,real8_1d(1:100))
   rc = nf90mpi_redef(ncidf)
   
   !double bottomDepth(nCells) ;
   !        bottomDepth:units = "m" ;
   !        bottomDepth:long_name = "Depth of the bottom of the ocean. Given as a positive distance from sea level." ;
   !METHOD: read from file produced by offline script 
   print*,'bottomdepth'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'bottomDepth', &
                            units='m',                    &
                            long_name="Depth of the bottom of the ocean. Given as a positive distance from sea level.", &
                            pdimids=(/dimidcf/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid  (ncid,'bottomDepth',varid)
   rc = nf90_get_var  (ncid,varid,real8_1d)
   rc = nf90mpi_put_var_all(ncidf,varid2,real8_1d)
   rc = nf90mpi_redef(ncidf)
   
   !double restingThickness(nCells, nVertLevels) ;
   !        restingThickness:units = "m" ;
   !        restingThickness:long_name = "Layer thickness when the ocean is at rest, i.e. without SSH or internal perturbations." ;
   !METHOD: compute (from bottomdepth and refbottomdepth) 
   print*,'restingthickness'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'restingThickness', &
                            units='m',                    &
                            long_name="Layer thickness when the ocean is at rest, i.e. without SSH or internal perturbations.", &
                            pdimids=(/dimidl, dimidcf/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid  (ncid,'restingThickness',varid)
   do k = 1,100
     rc = nf90_get_var  (ncid,varid,real8_1d,(/k,1/),(/1,ncellsf/))
     rc = nf90mpi_put_var_all(ncidf,varid2,real8_1d,(/k,1/)*1_mpi_offset_kind,(/1,ncellsf/)*1_mpi_offset_kind)
   enddo
   rc = nf90mpi_redef(ncidf)
   
   !int maxLevelCell(nCells) ;
   !        maxLevelCell:units = "unitless" ;
   !        maxLevelCell:long_name = "Index to the last active ocean cell in each column." ;
   !METHOD: compute (from restingthickness) 
   print*,'maxlevelcell'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'maxLevelCell',  &
                            units='unitless',                    &
                            long_name="Index to the last active ocean cell in each column.", &
                            pdimids=(/dimidcf/), pxtype=nf90mpi_int  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid  (ncid,'maxLevelCell',varid)
   rc = nf90_get_var  (ncid,varid,integer1d)
   rc = nf90mpi_put_var_all(ncidf,varid2,integer1d)
   rc = nf90mpi_redef(ncidf)
   

   !int minLevelCell(nCells) ;
   !        minLevelCell:units = "unitless" ;
   !        minLevelCell:long_name = "Index to the last active ocean cell in each column." ;
   !METHOD: compute (from restingthickness) 
   print*,'minlevelcell'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'minLevelCell',  &
                            units='unitless',                    &
                            long_name="Index to the first active ocean cell in each column.", &
                            pdimids=(/dimidcf/), pxtype=nf90mpi_int  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90_inq_varid  (ncid,'minLevelCell',varid)
   rc = nf90_get_var  (ncid,varid,integer1d)
   rc = nf90mpi_put_var_all(ncidf,varid2,integer1d)
   rc = nf90mpi_redef(ncidf)
   
   !double vertCoordMovementWeights(nVertLevels) ;
   !        vertCoordMovementWeights:units = "unitless" ;
   !        vertCoordMovementWeights:long_name = "Weights used for distribution of sea surface height perturbations through multiple vertical levels." ;
   !METHOD: copy 
   print*,'vertCoordMovementWeights'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'vertCoordMovementWeights', &
                            units='unitless',                    &
                            long_name="Weights used for distribution of sea surface height perturbations' &
                            //' through multiple vertical levels.", &
                            pdimids=(/dimidl/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_get_var_all  (ncidc,varid,real8_1d(1:100))
   rc = nf90mpi_put_var_all(ncidf,varid2,real8_1d(1:100))
   rc = nf90mpi_redef(ncidf)
 

   !double temperatureSurfaceValue(Time, nCells) ;
   !        temperatureSurfaceValue:long_name = "potential temperature extrapolated to ocean surface" ;
   !        temperatureSurfaceValue:units = "degrees Celsius" ;
   !METHOD: interpolate 
   print*,'temperatureSurfaceValue'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'temperatureSurfaceValue',  &
                            units='degrees Celsius',                    &
                            long_name="potential temperature extrapolated to ocean surface", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/),debug=.true. )
   rc = nf90mpi_redef(ncidf)

      
   !double salinitySurfaceValue(Time, nCells) ;
   !        salinitySurfaceValue:long_name = "salinity extrapolated to ocean surface" ;
   !        salinitySurfaceValue:units = "PSU" ;
   !METHOD: interpolate 
   print*,'salinitySurfaceValue'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'salinitySurfaceValue',  &
                            units='PSU',                    &
                            long_name="salinity extrapolated to ocean surface", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double surfaceVelocityZonal(Time, nCells) ;
   !        surfaceVelocityZonal:long_name = "Zonal surface velocity reconstructed at cell centers" ;
   !        surfaceVelocityZonal:units = "m s^{-1}" ;
   !METHOD: interpolate 
   print*,'surfaceVelocityZonal'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'surfaceVelocityZonal',  &
                            units='m s^{-1}',                    &
                            long_name="Zonal surface velocity reconstructed at cell centers", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double surfaceVelocityMeridional(Time, nCells) ;
   !        surfaceVelocityMeridional:long_name = "Meridional surface velocity reconstructed at cell centers" ;
   !        surfaceVelocityMeridional:units = "m s^{-1}" ;
   !METHOD: interpolate 
   print*,'surfaceVelocityMeridional'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'surfaceVelocityMeridional',  &
                            units='m s^{-1}',                    &
                            long_name="Meridional surface velocity reconstructed at cell centers", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double SSHGradientZonal(Time, nCells) ;
   !        SSHGradientZonal:long_name = "Zonal gradient of SSH reconstructed at cell centers" ;
   !        SSHGradientZonal:units = "m m^{-1}" ;
   !METHOD: interpolate 
   print*,'SSHGradientZonal'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'SSHGradientZonal',  &
                            units='m m^{-1}',                    &
                            long_name="Zonal gradient of SSH reconstructed at cell centers", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double SSHGradientMeridional(Time, nCells) ;
   !        SSHGradientMeridional:long_name = "Meridional gradient of SSH reconstructed at cell centers" ;
   !        SSHGradientMeridional:units = "m m^{-1}" ;
   !METHOD: interpolate 
   print*,'SSHGradientMeridional'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'SSHGradientMeridional',  &
                            units='m m^{-1}',                    &
                            long_name="Meridional gradient of SSH reconstructed at cell centers", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double vertNonLocalFluxTemp(Time, nCells, nVertLevelsP1) ;
   !        vertNonLocalFluxTemp:long_name = "CVMix/KPP: nonlocal boundary layer mixing term for temperature" ;
   !        vertNonLocalFluxTemp:units = "nondimensional" ;
   !METHOD: assign constant 0 
   print*,'vertNonLocalFluxTemp'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'vertNonLocalFluxTemp', &
                            units='nondimensional',                    &
                            long_name="CVMix/KPP: nonlocal boundary layer mixing term for temperature", &
                            pdimids=(/dimidv1,dimidcf,dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   real8_1d = 0.0
   rc = nf90mpi_put_var_all(ncidf,varid2,real8_1d)
   rc = nf90mpi_redef(ncidf)
   
   !int indMLD(Time, nCells) ;
   !        indMLD:units = "unitless" ;
   !        indMLD:long_name = "index of model where mixed layer depth occurs (always one past)" ;
   !METHOD: assign constant 3 
   print*,'indMLD'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'indMLD',  &
                            units='unitless',                    &
                            long_name="index of model where mixed layer depth occurs (always one past)", &
                            pdimids=(/dimidcf,dimidt/), pxtype=nf90mpi_int  )
   rc = nf90mpi_enddef(ncidf)
   integer1d = 3
   rc = nf90mpi_put_var_all(ncidf,varid2,integer1d)
   rc = nf90mpi_redef(ncidf)
   
   !double landIcePressure(Time, nCells) ;
   !        landIcePressure:units = "Pa" ;
   !        landIcePressure:long_name = "Pressure defined at the sea surface due to land ice." ;
   !METHOD: interpolate 
   print*,'landIcePressure'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'landIcePressure',  &
                            units='Pa',                    &
                            long_name="Pressure defined at the sea surface due to land ice.", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double landIceDraft(Time, nCells) ;
   !        landIceDraft:units = "m" ;
   !        landIceDraft:long_name = "The elevation of the interface between land ice and the ocean." ;
   !METHOD: interpolate 
   print*,'landIceDraft'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'landIceDraft',  &
                            units='m',                    &
                            long_name="The elevation of the interface between land ice and the ocean.", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double landIceFraction(Time, nCells) ;
   !        landIceFraction:units = "unitless" ;
   !        landIceFraction:long_name = "The fraction of each cell covered by land ice" ;
   !METHOD: interpolate 
   print*,'landIceFraction'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'landIceFraction',  &
                            units='unitless',                    &
                            long_name="The fraction of each cell covered by land ice.", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double seaIcePressure(Time, nCells) ;
   !        seaIcePressure:units = "Pa" ;
   !        seaIcePressure:long_name = "Pressure at the sea surface due to sea ice." ;
   !METHOD: interpolate 
   print*,'seaIcePressure'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'seaIcePressure',  &
                            units='Pa',                    &
                            long_name="Pressure at the sea surface due to sea ice.", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double atmosphericPressure(Time, nCells) ;
   !        atmosphericPressure:units = "Pa" ;
   !        atmosphericPressure:long_name = "Pressure at the sea surface due to the atmosphere." ;
   !METHOD: interpolate 
   print*,'atmosphericPressure'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'atmosphericPressure',  &
                            units='Pa',                    &
                            long_name="Pressure at the sea surface due to the atmosphere.", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double accumulatedFrazilIceMass(Time, nCells) ;
   !        accumulatedFrazilIceMass:units = "kg m^{-2}" ;
   !        accumulatedFrazilIceMass:long_name = "Mass per unit area of frazil ice produced. Reset to zero at each coupling interval" ;
   !METHOD: interpolate 
   print*,'accumulatedFrazilIceMass'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'accumulatedFrazilIceMass',  &
                            units='kg m^{-2}',                    &
                            long_name="Mass per unit area of frazil ice produced. Reset to zero at each coupling interval", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double accumulatedFrazilIceSalinity(Time, nCells) ;
   !        accumulatedFrazilIceSalinity:units = "kg m^{-2}" ;
   !        accumulatedFrazilIceSalinity:long_name = "Salinity associated with accumulatedFrazilIceMass. Reset to zero at each coupling interval" ;
   !METHOD: interpolate 
   print*,'accumulatedFrazilIceSalinity'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'accumulatedFrazilIceSalinity',  &
                            units='kg m^{-2}',                    &
                            long_name="Salinity associated with accumulatedFrazilIceMass. ' &
                            //'Reset to zero at each coupling interval" , &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double accumulatedLandIceFrazilMass(Time, nCells) ;
   !        accumulatedLandIceFrazilMass:units = "kg m^{-2}" ;
   !        accumulatedLandIceFrazilMass:long_name = "Mass per unit area of frazil ice produced under land ice.  Only computed when not coupled to a dynamic land-ice model." ;
   !METHOD: interpolate 
   print*,'accumulatedLandIceFrazilMass'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'accumulatedLandIceFrazilMass',  &
                            units='kg m^{-2}',                    &
                            long_name="Mass per unit area of frazil ice produced under land ice. ' &
                            //' Only computed when not coupled to a dynamic land-ice model.", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double frazilSurfacePressure(Time, nCells) ;
   !        frazilSurfacePressure:units = "Pa" ;
   !        frazilSurfacePressure:long_name = "surface pressure forcing due to weight of frazil ice" ;
   !METHOD: interpolate 
   print*,'frazilSurfacePressure'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'frazilSurfacePressure',  &
                            units='Pa',                    &
                            long_name="surface pressure forcing due to weight of frazil ice", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !double filteredSSHGradientZonal(Time, nCells) ;
   !        filteredSSHGradientZonal:units = "m m^{-1}" ;
   !        filteredSSHGradientZonal:long_name = "Time filtered zonal gradient of SSH" ;
   !METHOD: interpolate 
   print*,'filteredSSHGradientZonal'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'filteredSSHGradientZonal',  &
                            units='m m^{-1}',                    &
                            long_name="Time filtered zonal gradient of SSH", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)
        
   !double filteredSSHGradientMeridional(Time, nCells) ;
   !        filteredSSHGradientMeridional:units = "m m^{-1}" ;
   !        filteredSSHGradientMeridional:long_name = "Time filtered meridional gradient of SSH" ;
   !METHOD: interpolate 
   print*,'filteredSSHGradientMeridional'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'filteredSSHGradientMeridional',  &
                            units='m m^{-1}',                    &
                            long_name="Time filtered meridional gradient of SSH", &
                            pdimids=(/dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
      CALL interpolate_and_write_variable(1, ncidc, varid,        &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/1,1/),         &
                                          countin=(/ncellsc,1/), &
                                          countout=(/ncellsf,1/) )
   rc = nf90mpi_redef(ncidf)

   !char forcingGroupNames(Time, nForcingGroupsMax, StrLen) ;
   !METHOD: copy
   print*,'forcingGroupNames'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'forcingGroupNames',  &
                            pdimids=(/dimidsl, dimidfg, dimidt/), pxtype=nf90mpi_char  )
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,5
      rc = nf90mpi_get_var_all  (ncidc,varid,time,start=(/1,k,1/)*1_mpi_offset_kind,count=(/64,1,1/)*1_mpi_offset_kind)
      rc = nf90mpi_put_var_all(ncidf,varid2,time,start=(/1,k,1/)*1_mpi_offset_kind,count=(/64,1,1/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)
        
   !char forcingGroupRestartTimes(Time, nForcingGroupsMax, StrLen) ;
   !METHOD: copy
   print*,'forcingGroupRestartTimes'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'forcingGroupRestartTimes',  &
                            pdimids=(/dimidsl, dimidfg, dimidt/), pxtype=nf90mpi_char  )
   rc = nf90mpi_enddef(ncidf)
   DO k = 1,5
      rc = nf90mpi_get_var_all  (ncidc,varid,time,start=(/1,k,1/)*1_mpi_offset_kind,count=(/64,1,1/)*1_mpi_offset_kind)
      rc = nf90mpi_put_var_all(ncidf,varid2,time,start=(/1,k,1/)*1_mpi_offset_kind,count=(/64,1,1/)*1_mpi_offset_kind)
   ENDDO
   rc = nf90mpi_redef(ncidf)
        
print*,'done copied ocn variables'


! SSH related variables (layerthickness)
   !double layerThickness(Time, nCells, nVertLevels) ;
   !        layerThickness:units = "m" ;
   !        layerThickness:long_name = "layer thickness" ;
   !METHOD: compute SSH from coarse grid layer thickness and bottom depth
   !        interpolate SSH, then use it to adjust fine grid resting thickness   
   print*,'layerThickness'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'layerThickness',   &
                            units='m',                    &
                            long_name="layer thickness", &
                            pdimids=(/dimidl, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_sync(ncidf)
   rc = nf90mpi_enddef(ncidf)
   
   allocate(areacell(ncellsf))
   allocate(bottomdepth_in(ncellsc))
   allocate(ssh_in(ncellsc))
   allocate(layerthickness_in(100,ncellsc))
   rc = nf90mpi_inq_varid  (ncidc,'bottomDepth',varid)
   rc = nf90mpi_get_var_all(ncidc,varid,bottomdepth_in)
      print*,'ssh_in ',minval(bottomdepth_in),maxval(bottomdepth_in)
   rc = nf90mpi_inq_varid  (ncidc,'layerThickness',varid)
   rc = nf90mpi_get_var_all(ncidc,varid,layerthickness_in)
   rc = nf90_inq_varid  (ncid,'areaCell',varid)
   rc = nf90_get_var(ncid,varid,areaCell)
   ! compute coarse grid ssh
   do n = 1,ncellsc
      ssh_in(n) = sum(layerthickness_in(:,n)) - bottomdepth_in(n) 
   enddo
      print*,'ssh_in ',minval(ssh_in),maxval(ssh_in)
      !copied interpolation code to interpolate ssh to fine grid
      real8_1d = 0.0
      DO n = 1,ncellsf
         DO l = 1,3
            IF(remap_indices(l,n) > 0)     &
               real8_1d(n) = real8_1d(n) +     &
                             ssh_in(remap_indices(l,n)) * remap_wgts(l,n)
         ENDDO
      ENDDO
      
      print*,'ssh_terp ',minval(real8_1d),maxval(real8_1d)
      ! reset global mean of fine grid ssh to zero
      wrk1 = sum(areacell)
      wrk2 = 0.0
      do n = 1,ncellsf
         wrk2 = wrk2 + real8_1d(n)*areacell(n)
      enddo
      wrk2 = wrk2 / wrk1
         print*,'wrk1,wrk2',wrk1,wrk2
      real8_1d = real8_1d - wrk2
   deallocate(ssh_in, bottomdepth_in, layerthickness_in, areacell)

   allocate(wrkarr(ncellsf))      
   allocate(restingthickness(ncellsf))      
   rc = nf90_inq_varid  (ncid,'restingThickness',varid)
   wrkarr = 0.0
   do k = 1,100
     rc = nf90_get_var(ncid,varid,restingthickness,(/k,1/),(/1,ncellsf/))
     wrkarr = wrkarr + restingthickness
   enddo
   wrkarr = (wrkarr + real8_1d )/ wrkarr
   do k = 1,100
     rc = nf90_get_var(ncid,varid,restingthickness,(/k,1/),(/1,ncellsf/))
     restingthickness = restingthickness * wrkarr
     rc = nf90mpi_put_var_all(ncidf,varid2,restingthickness,(/k,1,1/)*1_mpi_offset_kind,(/1,ncellsf,1/)*1_mpi_offset_kind)
   enddo
   deallocate(wrkarr,restingthickness)
   rc = nf90mpi_redef(ncidf)
   
   rc = nf90mpi_sync(ncidf)
! velocities
   print*,'solve_for_velocity_and_write'
   call solve_for_velocity_and_write0()

   !rc = nf90mpi_close(ncidc)
   rc = nf90mpi_sync(ncidf)


   !rc = nf90mpi_open(mpi_comm_world,trim(coarsefile),nf90mpi_nowrite,info,ncidc)
   !print*,'open coarse restart',rc
   if(rc .ne. 0)print*,trim(coarsefile)

   !double temperature(Time, nCells, nVertLevels) ;
   !        temperature:long_name = "potential temperature" ;
   !        temperature:units = "degrees Celsius" ;
   !METHOD: interpolate 
   print*,'temperature'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'temperature',  &
                            units='degrees Celsius',                    &
                            long_name="potential temperature", &
                            pdimids=(/dimidl, dimidcf, dimidt/), pxtype=nf90mpi_double, debug=.true.  )
   rc = nf90mpi_enddef(ncidf)
   real8_1d(1:100) = 0.
   do k = 1,100
     rc = nf90_get_var  (ncid,varid,real8_1d,(/k,1/),(/1,ncellsf/))
     rc = nf90mpi_put_var_all(ncidf,varid2,real8_1d,(/k,1/)*1_mpi_offset_kind,(/1,ncellsf/)*1_mpi_offset_kind)
      CALL interpolate_and_write_variable3d(k, ncidc, varid,       &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/k,1,1/),         &
                                          countin=(/1,ncellsc,1/), &
                                          countout=(/1,ncellsf,1/), &
                                          spval = -9.99999979021477d+33 )
      !rc = nf90mpi_put_var_all(ncidf,varid2,real8_1d(1:100), (/1,k,1/)*1_mpi_offset_kind, (/100,1,1/)*1_mpi_offset_kind)
   enddo
   rc = nf90mpi_sync(ncidf)
   rc = nf90mpi_redef(ncidf)
  

   !double salinity(Time, nCells, nVertLevels) ;
   !        salinity:long_name = "salinity" ;
   !        salinity:units = "grams salt per kilogram seawater" ;
   !METHOD: interpolate 
   print*,'salinity'
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'salinity',  &
                            units='grams salt per kilogram seawater',                    &
                            long_name="salinity", &
                            pdimids=(/dimidl, dimidcf, dimidt/), pxtype=nf90mpi_double  )
   rc = nf90mpi_enddef(ncidf)
   do k = 1,100
     rc = nf90_get_var  (ncid,varid,real8_1d,(/k,1/),(/1,ncellsf/))
     rc = nf90mpi_put_var_all(ncidf,varid2,real8_1d,(/k,1/)*1_mpi_offset_kind,(/1,ncellsf/)*1_mpi_offset_kind)
      CALL interpolate_and_write_variable3d(k, ncidc, varid,       &
                                          ncidf,           &
                                          varid2,        &
                                          start=(/k,1,1/),         &
                                          countin=(/1,ncellsc,1/), &
                                          countout=(/1,ncellsf,1/), &
                                          spval = -9.99999979021477d+33)
      !rc = nf90mpi_put_var_all(ncidf,varid2,real8_1d(1:100), (/1,k,1/)*1_mpi_offset_kind, (/100,1,1/)*1_mpi_offset_kind)
   enddo
   rc = nf90mpi_redef(ncidf)
   print*,'close'
   
! close up the output file
   rc = nf90mpi_close(ncidf)
    print*,'closed'
  
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
      
      rc = nf90mpi_inq_varid(ncidin ,TRIM(varname),varid  )
      if(PRESENT(debug))PRINT*,'rc: input inq_varid',trim(varname),rc
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
      
      if(present(debug))print*,'lev ',lev
      if(present(debug))print*,'ncidr ',ncidr
      if(present(debug))print*,'varid ',varid
      if(present(debug))print*,'pncid ',pncid
      if(present(debug))print*,'pvarid ',pvarid
      if(present(debug))print*,'start ',start
      if(present(debug))print*,'countin ',countin
      if(present(debug))print*,'countout ',countout
      
      
      
      if(lev.ne.1)stop 'unavailable lev'
      rc = nf90mpi_get_var_all(ncidr, varid, real8_1d_in, start*1_mpi_offset_kind, countin*1_mpi_offset_kind)
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
!print*,'remap_wghts ',remap_wgts(:,111:113)
!print*,'remap_indices ',remap_indices(:,111:113)
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
      ENDDO
!$OMP end parallel

      !print*,'output ',real8_1d(1563)
      rc = nf90mpi_put_var_all(pncid,pvarid,real8_1d, start*1_mpi_offset_kind, countout*1_mpi_offset_kind)
      if(PRESENT(debug))PRINT*,'rc: interpolate_and_write_variable nfmpi_put_var',rc
                             
      END SUBROUTINE interpolate_and_write_variable


!------------------------------------------   

   subroutine solve_for_velocity_and_write0()
   
   integer nb, ne
   real*8, dimension(:,:), allocatable :: velocity
   real*8, dimension(:), allocatable :: velocity0
   allocate(velocity(1:100,10000))
   velocity = 0.
print*,'write velocity'
   ! write
   !double normalVelocity(Time, nEdges, nVertLevels) ;
   !        normalVelocity:units = "m s^{-1}" ;
   !        normalVelocity:long_name = "horizontal velocity, normal component to an edge" ;
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'normalVelocity',  &
                            units='m s^{-1}',                    &
                            long_name="horizontal velocity, normal component to an edge", &
                            pdimids=(/dimidl, dimidef, dimidt/), pxtype=nf90mpi_double, debug=.true.  )

 print*,'after define_rst_variable'
  rc = nf90mpi_enddef(ncidf)
   do nb =1,nedgesf,10000
      ne = min(nb+9999,nedgesf)
      print*,'nb,ne',nb,ne,(ne+1-nb)
      rc = nf90mpi_put_var_all(ncidf,varid2,velocity(1:100,1:(ne+1-nb)))
   enddo
   rc = nf90mpi_sync(ncidf)
   rc = nf90mpi_redef(ncidf)

 print*,'deallocate-allocate'
   deallocate(velocity)
   allocate(velocity0(nedgesf))
   velocity0 = 0.
   !double normalBarotropicVelocity(Time, nEdges) ;
   !        normalBarotropicVelocity:units = "m s^{-1}" ;
   !        normalBarotropicVelocity:long_name = "barotropic velocity, used in split-explicit time-stepping" ;
   call define_rst_variable(varid, varid2,                    &
                            ncidc, ncidf,              &
                            'normalBarotropicVelocity',  &
                            units='m s^{-1}',                    &
                            long_name="barotropic velocity, used in split-explicit time-stepping", &
                            pdimids=(/dimidef, dimidt/), pxtype=nf90mpi_double, debug=.true.  )
 print*,'after define_rst_variable'
   rc = nf90mpi_enddef(ncidf)
   rc = nf90mpi_put_var_all(ncidf,varid2,velocity0(:))
   rc = nf90mpi_redef(ncidf)
   
   deallocate(velocity0)
 print*,'end of solve_for_velocity_and_write0'

   end subroutine solve_for_velocity_and_write0

      SUBROUTINE interpolate_and_write_variable3d(       &
                             lev,                     &
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
      !integer, intent(in) :: nctem
      INTEGER, INTENT(IN), DIMENSION(:) :: start, countin, countout
      REAL*8, INTENT(IN), OPTIONAL :: spval
      REAL*8, INTENT(IN), OPTIONAL :: zlim
      LOGICAL, INTENT(IN), OPTIONAL :: debug
      
      ! local variables
      INTEGER :: n, l  ! loop indices
      INTEGER :: rc, varidt

      REAL*8 :: spval_tem, zlim_tem
      real*8 :: wgtsum
      
      rc = nf90mpi_get_var_all(ncidr, varid, real8_1d_in, start*1_mpi_offset_kind, countin*1_mpi_offset_kind)
      if(Present(debug))print*,'input ',real8_1d_in(1:10)
      if(PRESENT(debug))PRINT*,'rc: interpolate_and_write_variable get_var',rc
!!$OMP parallel num_threads(1)
      if(present(spval)) then
         !!$OMP do
         !DO n = 1,ncellsf
         !   real8_1d(n) = spval
         !enddo
         spval_tem = spval
      else
         !!$OMP do
         !DO n = 1,ncellsf
         !   real8_1d(n) = 0.0
         !enddo
         spval_tem = 0.0
      endif
         !!$OMP do
         DO n = 1,ncellsf
            real8_1d(n) = 0.0
         enddo
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
      if(PRESENT(debug))PRINT*,'rc: interpolate_and_write_variable get remap',rc
      rc = nf90_inq_varid(nctem,'remap_indices',varidt)
      rc = nf90_get_var(nctem,varidt,remap_indices3d,(/1,1,lev/),(/3,ncellsf,1/) )
      rc = nf90_inq_varid(nctem,'remap_wgts',varidt)
      rc = nf90_get_var(nctem,varidt,remap_wgts3d,(/1,1,lev/),(/3,ncellsf,1/) )
      if(PRESENT(debug))PRINT*,'rc: top of loop',rc
      !$OMP do private(l, wgtsum)
      DO n = 1,ncellsf
         wgtsum = 0.0
         DO l = 1,3
            IF(remap_indices3d(l,n) > 0) then
               !if(real8_1d(n) == spval_tem) real8_1d(n) = 0.0
               if( abs(real8_1d_in(remap_indices3d(l,n))) >= zlim_tem ) then
                  real8_1d(n) = real8_1d(n) +     &
                                real8_1d_in(remap_indices3d(l,n)) * remap_wgts3d(l,n)
                  wgtsum = wgtsum + remap_wgts3d(l,n)
               endif
            endif
         ENDDO
         if(wgtsum > 0.0) then
            real8_1d(n) = real8_1d(n) / wgtsum
         else
            real8_1d(n) = spval_tem
         endif
      ENDDO
!$OMP end parallel
      if(PRESENT(debug))PRINT*,'rc: end of loop',rc

      !print*,'output ',real8_1d(1563)
      rc = nf90mpi_put_var_all(pncid,pvarid,real8_1d, start*1_mpi_offset_kind, countout*1_mpi_offset_kind)
      if(PRESENT(debug))PRINT*,'rc: interpolate_and_write_variable nfmpi_put_var',rc
                             
      END SUBROUTINE interpolate_and_write_variable3d

!------------------------------------------
   

end program build_ocn_restart_file
