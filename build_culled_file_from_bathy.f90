program build_culled_file_from_bathy

! compile with mpi, hdf5, netcdf, and pnetcdf libraries

use netcdf
use pnetcdf, only : nf90mpi_create, nf90mpi_close, nf90mpi_enddef, nf90mpi_def_dim, &
                    nf90mpi_put_att, nf90mpi_def_var, nf90mpi_put_var_all, nf90mpi_inq_varid, &
                    nf90mpi_put_var, nf90mpi_redef, &
                    nf90mpi_clobber => NF90_CLOBBER, NF90mpi_64BIT_DATA=>NF90_64BIT_DATA, &
                    nf90mpi_unlimited => nf90_unlimited, NF90mpi_GLOBAL=>NF90_GLOBAL, &
                    NF90mpi_DOUBLE=>NF90_DOUBLE, NF90mpi_INT=>NF90_INT, nf90mpi_char=>nf90_char, &
                    nf90mpi_begin_indep_data, nf90mpi_end_indep_data, nf90mpi_sync
use mpi

implicit none

!!!
! Declarations

   ! master ocean data files
   character(len=128) :: infile, tfile, sfile, rfile
   !resolution specific variables
   integer caseres
   character(len=128) :: atmfile, outfile, poutfile

   ! netcdf-related variables
   !   'a' suffix is atmfile
   !   'o' suffix is ocean outfile
   !   no  suffix is input ocean master file
   integer :: rc, ncida, ncido, ncid, ncidp, ncids, ncidt
   integer :: dimid, varid, attid, varid2, varid3,pvarid, ncidi
   integer :: dimide, dimidc, dimidv, dimidt, dimidme, dimid2, dimidvd, dimidl, dimidme2
   integer :: pdimide, pdimidc, pdimidv, pdimidt, pdimidme, pdimid2, pdimidvd, pdimidl, pdimidme2
   integer :: natts
   character(len=80) :: attname

   ! grid variables
   ! vertices
   integer :: ncellsa, ncells, ncellso, maxlinks, nverticesa, nverticeso, nedgesa, nedgeso
   real*8, allocatable, dimension(:) :: latcello, loncello
   real*8, allocatable, dimension(:) :: latcella, loncella, latcell, loncell, &
                                        areacell, areacella, areaatm
   integer, allocatable, dimension(:) :: numlinks, sumarea         
   integer, allocatable, dimension(:,:) :: linktoatm         
   real*8, allocatable, dimension(:,:) :: wgts         
   real*8, allocatable, dimension(:,:) :: xyzcella, xyzcell
   real*8 :: dotmax, dotloc, ofrac_glob, globareainv, arealim
   integer :: nnear
   integer, allocatable, dimension(:) :: cellexistsa, vertexexistsa, edgeexistsa
   integer, allocatable, dimension(:,:) :: cellsonedgea, cellsonvertexa, cellsonedge, cellsonvertex
   integer, allocatable, dimension(:) :: indexcella2o, indexverta2o, indexedgea2o
   integer, allocatable, dimension(:) :: indexcello2a, indexverto2a, indexedgeo2a
   integer, allocatable, dimension(:,:) :: cellsoncella
   logical, allocatable, dimension(:) :: angleswap
   integer :: temcells(6), ncellstem, count1, countloc
   
   real*8 r8min, r8max, wrk1, wrk2, areatem, sphere_radius
   real*8, dimension(100) :: refbottomdepth, refzmid
   integer :: ne, nlen, nb
   
   ! generic input arrays
   integer, allocatable, dimension(:) :: integer1d_in
   integer, allocatable, dimension(:,:) :: integer2d_in
   real*8, allocatable, dimension(:) :: real8_1d_in
   real*8, allocatable, dimension(:,:) :: real8_2d_in

   ! generic atm arrays
   integer, allocatable, dimension(:) :: integer1d_a, nedgesoncella
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
   real*8, allocatable, dimension(:) :: real8_edge, angleedge
   real*8, allocatable, dimension(:,:) :: real8_2d_edge
   
   ! other input
   real*8, allocatable, dimension(:) :: bathy
   
   ! output arrays
   real*8, allocatable, dimension(:) :: bottomdepth
   real*8, allocatable, dimension(:,:) :: restingthickness
   real*8, allocatable, dimension(:) :: areatriangle
   real*8, allocatable, dimension(:,:) :: resting

   ! loop indices
   integer :: n, na, m, id, k, iii, kk, indxo
      
   real*8 pi
   
   integer start2(2), start3(3), count2(2), count3(3)
   
   ! mpi variables
   integer :: err, rank, nproc, info, nctem
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

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nproc, err)

          ! set an MPI-IO hint to disable file offset alignment for
          ! fixed-size variables
          call MPI_Info_create(info, err)
          call MPI_Info_set(info, "nc_var_align_size", "1", err)

   pi = 4._8 * atan(1._8)
!!!
! preliminaries

   ! specify input and output files
   print*,'enter the global base mesh file name'
   read(5,'(a128)') atmfile
   print*,'enter the global base mesh bathymetry file name'
   read(5,'(a128)') infile
   print*,'enter the ocean mesh file name'
   read(5,'(a128)') outfile

   rc = nf90mpi_create(mpi_comm_world,trim(outfile),IOR(NF90_CLOBBER, NF90_64BIT_DATA),info,ncido)
print*,'dbg create ',rc
   if(rc .ne. 0) stop 'troubles creating the ocean mesh file - stopping'
          call MPI_Info_free(info, err)

   ! open the necessary netcdf files
   rc = nf90_open(trim(infile),nf90_nowrite,ncid)
   if(rc .ne. 0) stop 'troubles opening the global base mesh bathymetry file - stopping'
   rc = nf90_open(trim(atmfile),nf90_nowrite,ncida)
   if(rc .ne. 0) stop 'troubles opening the global base mesh file - stopping'

!!!
! We are going to do nearest neighbor interpolation.

   ! allocate and read ocean and atm cell locations
   rc = nf90_inq_dimid(ncida,'nCells',dimid)
   rc = nf90_inquire_dimension(ncida,dimid,len=ncellsa)
   rc = nf90_inq_dimid(ncida,'nEdges',dimid)
   rc = nf90_inquire_dimension(ncida,dimid,len=nedgesa)
   rc = nf90_inq_dimid(ncida,'nVertices',dimid)
   rc = nf90_inquire_dimension(ncida,dimid,len=nverticesa)
   allocate(latcella(ncellsa))
   allocate(loncella(ncellsa))
   allocate(areacella(ncellsa))
   allocate(cellsoncella(6,ncellsa))
   allocate(areaatm (ncellsa))
   allocate(bathy (ncellsa))
   rc = nf90_inq_varid(ncida,'latCell',varid)
   rc = nf90_get_var  (ncida,varid,latcella)
   rc = nf90_inq_varid(ncida,'lonCell',varid)
   rc = nf90_get_var  (ncida,varid,loncella)
   rc = nf90_inq_varid(ncid,'bathymetry',varid)
   rc = nf90_get_var  (ncid,varid,bathy)
   rc = nf90_inq_varid(ncida,'areaCell',varid)
   rc = nf90_get_var  (ncida,varid,areacella)
   rc = nf90_inq_varid(ncida,'cellsOnCell',varid)

! convert land values of bathy to 1, and define south of 75S as land.
do n = 1,ncellsa
   rc = nf90_get_var  (ncida,varid,cellsoncella(:,n),(/1,n/),(/6,1/))
   if(latcella(n) < -75.*pi/180. .and. bathy(n) < -5990. ) bathy(n) = 1.
enddo

!! 15km adjustment - this needs explanation
!   bathy(2341120) = 1.0
   
   allocate(cellexistsa(ncellsa))
   allocate(edgeexistsa(nedgesa))
   allocate(vertexexistsa(nverticesa))

   allocate(integer1d_a(ncellsa))
   allocate(nedgesoncella(ncellsa))
   rc = nf90_inq_varid(ncida,'nEdgesOnCell',varid)
   do n = 1,ncellsa,1000
      ne = min(ncellsa,n+999)
      nlen = min(1000,ncellsa+1-n)
      rc = nf90_get_var  (ncida,varid,nedgesoncella(n:ne),(/n/),(/nlen/))
   enddo

!!!!!! iteration to fill one edge cells
   count1 = 99999
   do while (count1 > 0)

      vertexexistsa = 0
      cellexistsa = 0
      edgeexistsa = 0
      print*,'loop through bathymetry'
 
   
   ! loop through the atm grid and identify the ocean cells as bathy < -1
      ncellso = 0
      do na = 1,ncellsa
         if(bathy(na) < -1.) then
            ncellso = ncellso + 1
            cellexistsa(na) = 1
         endif
      enddo
   
      print*,'ncellsa/o',ncellsa,ncellso

      !!!
      ! determine number of vertices and edges, and their mapping from the atmfile
      allocate(cellsonedgea(2,nedgesa))
      allocate(cellsonvertexa(3,nverticesa))
      rc = nf90_inq_varid(ncida,'cellsOnEdge',varid)
      do n = 1,nedgesa
         rc = nf90_get_var  (ncida,varid,cellsonedgea(:,n),(/1,n/),(/2,1/))
      enddo
      rc = nf90_inq_varid(ncida,'cellsOnVertex',varid)
      do n = 1,nverticesa
         rc = nf90_get_var  (ncida,varid,cellsonvertexa(:,n),(/1,n/),(/3,1/))
      enddo
      do na = 1,nverticesa
         do m = 1,3
            id = cellsonvertexa(m,na)
            if(id>0) then
               if(cellexistsa(id)==1) vertexexistsa(na) = 1
            endif
         enddo
      enddo
      nverticeso = sum(vertexexistsa)

      do na = 1,nedgesa
         do m = 1,2
            id = cellsonedgea(m,na)
            if(id>0) then
               if(cellexistsa(id)==1) edgeexistsa(na) = 1
            endif
         enddo
      enddo
      nedgeso = sum(edgeexistsa)

      allocate(indexcella2o(0:ncellsa))
      allocate(indexedgea2o(0:nedgesa))
      allocate(indexverta2o(0:nverticesa))
      indexcella2o = 0
      indexedgea2o = 0
      indexverta2o = 0
      allocate(indexcello2a(0:ncellso))
      allocate(indexedgeo2a(0:nedgeso))
      allocate(indexverto2a(0:nverticeso))
      indexcello2a = 0
      indexedgeo2a = 0
      indexverto2a = 0

      n = 0
      do na = 1,ncellsa
         if(cellexistsa(na)==1) then
            n = n+1
            indexcella2o(na) = n
            indexcello2a(n) = na
         endif
      enddo
      n = 0
      do na = 1,nedgesa
         if(edgeexistsa(na)==1) then
            n = n+1
            indexedgea2o(na) = n
            indexedgeo2a(n) = na
         endif
      enddo
      n = 0
      do na = 1,nverticesa
         if(vertexexistsa(na)==1) then
            n = n+1
            indexverta2o(na) = n
            indexverto2a(n) = na
         endif
      enddo

      ! intermediate developmental test writes
      print*,'ncellsa',ncellsa
      print*,'nedgesa',nedgesa
      print*,'nverticesa',nverticesa
      print*,'ncellso',ncellso
      print*,'nedgeso',nedgeso
      print*,'nverticeso',nverticeso
 
      allocate(integer1d(ncellso))
      allocate(integer2d(6,ncellso))
      allocate(latcello(ncellso))

      ! redefine cellsoncell
      count2 = (/6,1/)
      do n = 1,ncellso
         start2 = (/1,n/)
         do k = 1,6
            integer1d(k) = indexcella2o(cellsoncella(k,indexcello2a(n)))
            integer2d(k,n) = integer1d(k)
         enddo
      enddo

      call compress_1d_int(ncellsa, ncellso, nedgesoncella, integer1d, cellexistsa)
      call compress_1d_r8(ncellsa, ncellso, latcella, latcello, cellexistsa)
      print*,'nedgesoncella ',nedgesoncella(1:10)
      print*,'nedgesoncell ',integer1d(1:10)
      
      count1 = 0
      do n = 1,ncellso
         countloc = 0
         do m = 1,integer1d(n)
            if(integer2d(m,n)>0) countloc = countloc + 1
         enddo
         if(countloc==1 .and. abs(latcello(n))>0.25*pi) count1 = count1 + 1
      enddo
      
      print*,count1,' cells to be filled in'
      
      if(count1 > 0) then
         do n = 1,ncellso
         countloc = 0
            do m = 1,integer1d(n)
               if(integer2d(m,n)>0) countloc = countloc + 1
            enddo
            if(countloc==1 .and. abs(latcello(n))>0.25*pi) then
               na = indexcello2a(n)
               bathy(na) = 1.
            endif
         enddo
      endif
      
      deallocate(integer1d,integer2d,latcello)
      if(count1 > 0) then
         deallocate(indexcella2o)
         deallocate(indexedgea2o)
         deallocate(indexverta2o)
		 deallocate(indexcello2a)
		 deallocate(indexedgeo2a)
		 deallocate(indexverto2a)
		 deallocate(cellsonedgea)
		 deallocate(cellsonvertexa)
      endif
   enddo
   
!!!!!!!end fill iteration 

write(caseres)indexcello2a
write(caseres)indexcella2o
write(caseres)indexedgeo2a
write(caseres)indexedgea2o
write(caseres)indexverto2a
write(caseres)indexverta2o

   allocate(latcello(ncellso))
   allocate(loncello(ncellso))

   allocate(integer1d(ncellso))
   allocate(integer1d_edge(nedgeso))
   allocate(angleswap(nedgeso))
   allocate(integer1d_vert(nverticeso))
   allocate(integer2d(100,ncellso))
   allocate(real8_1d(ncellso))
   allocate(real8_2d(100,ncellso))
   allocate(real8_edge(nedgeso))
   allocate(angleedge(nedgeso))
   !allocate(real8_2d_edge(100,nedgeso))
   allocate(real8_vert(nverticeso))

   allocate(bottomdepth(ncellso))
   allocate(areatriangle(nverticeso))
   allocate(restingthickness(100,ncellso))
   
   !allocate(integer1d_in(ncells))
   !allocate(integer2d_in(100,ncells))
   !allocate(real8_1d_in(ncells))
   !allocate(real8_2d_in(100,ncells))
   
   allocate(integer1d_edge_a(nedgesa))
   allocate(integer1d_vert_a(nverticesa))
   allocate(integer2d_a(100,ncellsa))
   allocate(real8_1d_a(ncellsa))
   allocate(real8_2d_a(100,ncellsa))
   allocate(real8_edge_a(nedgesa))
   allocate(real8_vert_a(nverticesa))
   
   
!!!
! start defining the output file
   rc = nf90mpi_def_dim(ncido,'nEdges',nedgeso*1_mpi_offset_kind,dimide)
   rc = nf90mpi_def_dim(ncido,'nCells',ncellso*1_mpi_offset_kind,dimidc)
   rc = nf90mpi_def_dim(ncido,'nVertices',nverticeso*1_mpi_offset_kind,dimidv)
   rc = nf90mpi_def_dim(ncido,'Time',nf90mpi_unlimited*1_mpi_offset_kind,dimidt)
   rc = nf90mpi_def_dim(ncido,'maxEdges',6*1_mpi_offset_kind,dimidme)
   rc = nf90mpi_def_dim(ncido,'TWO',2*1_mpi_offset_kind,dimid2)
   rc = nf90mpi_def_dim(ncido,'vertexDegree',3*1_mpi_offset_kind,dimidvd)
   rc = nf90mpi_def_dim(ncido,'nVertLevels',100*1_mpi_offset_kind,dimidl)
   rc = nf90mpi_def_dim(ncido,'maxEdges2',12*1_mpi_offset_kind,dimidme2)
   

!   rc = nf90_inquire(ncid,nattributes=natts)
!   do n = 1,natts
!      rc = nf90_inq_attname(ncid, NF90_GLOBAL, n, attname)
!      rc = nf90_copy_att(ncid, NF90_GLOBAL, trim(attname), ncido, NF90_GLOBAL)
!   enddo

   rc = nf90mpi_def_var(ncido,'cellsOnCell',NF90mpi_INT, (/dimidme,dimidc/), varid)
      rc = nf90mpi_def_var(ncido,'bottomDepth',NF90mpi_DOUBLE, (/dimidc/), varid)
      rc = nf90mpi_def_var(ncido,'lonCell',NF90mpi_DOUBLE, (/dimidc/), varid)
      rc = nf90mpi_def_var(ncido,'latCell',NF90mpi_DOUBLE, (/dimidc/), varid)
      rc = nf90mpi_def_var(ncido,'areaCell',NF90mpi_DOUBLE, (/dimidc/), varid)
    rc = nf90mpi_def_var(ncido,'cellsOnEdge',NF90mpi_INT, (/dimid2,dimide/), varid)
      rc = nf90mpi_def_var(ncido,'angleEdge',NF90mpi_DOUBLE, (/dimide/), varid)
   rc = nf90mpi_def_var(ncido,'cellsOnVertex',NF90mpi_INT, (/dimidvd,dimidv/), varid)
      rc = nf90mpi_def_var(ncido,'dcEdge',NF90mpi_DOUBLE, (/dimide/), varid)
      rc = nf90mpi_def_var(ncido,'dvEdge',NF90mpi_DOUBLE, (/dimide/), varid)
   rc = nf90mpi_def_var(ncido,'edgesOnCell',NF90mpi_INT, (/dimidme,dimidc/), varid)
   rc = nf90mpi_def_var(ncido,'edgesOnEdge',NF90mpi_INT, (/dimidme2,dimide/), varid)
   rc = nf90mpi_def_var(ncido,'edgesOnVertex',NF90mpi_INT, (/dimidvd,dimidv/), varid)
      rc = nf90mpi_def_var(ncido,'indexToCellID',NF90mpi_INT, (/dimidc/), varid)
      rc = nf90mpi_def_var(ncido,'indexToVertexID',NF90mpi_INT, (/dimidv/), varid)
      rc = nf90mpi_def_var(ncido,'indexToEdgeID',NF90mpi_INT, (/dimide/), varid)
   rc = nf90mpi_def_var(ncido,'kiteAreasOnVertex',NF90mpi_DOUBLE, (/dimidvd,dimidv/), varid)
      rc = nf90mpi_def_var(ncido,'areaTriangle',NF90mpi_DOUBLE, (/dimidv/), varid)
      rc = nf90mpi_def_var(ncido,'latEdge',NF90mpi_DOUBLE, (/dimide/), varid)
      rc = nf90mpi_def_var(ncido,'latVertex',NF90mpi_DOUBLE, (/dimidv/), varid)
      rc = nf90mpi_def_var(ncido,'lonEdge',NF90mpi_DOUBLE, (/dimide/), varid)
      rc = nf90mpi_def_var(ncido,'lonVertex',NF90mpi_DOUBLE, (/dimidv/), varid)
      rc = nf90mpi_def_var(ncido,'meshDensity',NF90mpi_DOUBLE, (/dimidc/), varid)
      rc = nf90mpi_def_var(ncido,'nEdgesOnCell',NF90mpi_INT, (/dimidc/), varid)
      rc = nf90mpi_def_var(ncido,'nEdgesOnEdge',NF90mpi_INT, (/dimide/), varid)
   rc = nf90mpi_def_var(ncido,'verticesOnCell',NF90mpi_INT, (/dimidme,dimidc/), varid)
   rc = nf90mpi_def_var(ncido,'verticesOnEdge',NF90mpi_INT, (/dimid2,dimide/), varid)
   rc = nf90mpi_def_var(ncido,'weightsOnEdge',NF90mpi_DOUBLE, (/dimidme2,dimide/), varid)
      rc = nf90mpi_def_var(ncido,'xCell',NF90mpi_DOUBLE, (/dimidc/), varid)
      rc = nf90mpi_def_var(ncido,'xEdge',NF90mpi_DOUBLE, (/dimide/), varid)
      rc = nf90mpi_def_var(ncido,'xVertex',NF90mpi_DOUBLE, (/dimidv/), varid)
      rc = nf90mpi_def_var(ncido,'yCell',NF90mpi_DOUBLE, (/dimidc/), varid)
      rc = nf90mpi_def_var(ncido,'yEdge',NF90mpi_DOUBLE, (/dimide/), varid)
      rc = nf90mpi_def_var(ncido,'yVertex',NF90mpi_DOUBLE, (/dimidv/), varid)
      rc = nf90mpi_def_var(ncido,'zCell',NF90mpi_DOUBLE, (/dimidc/), varid)
      rc = nf90mpi_def_var(ncido,'zEdge',NF90mpi_DOUBLE, (/dimide/), varid)
      rc = nf90mpi_def_var(ncido,'zVertex',NF90mpi_DOUBLE, (/dimidv/), varid)

   !rc = nf90mpi_def_var(ncido,'salinity',nf90mpi_DOUBLE, (/dimidl,dimidc,dimidt/), varid)
      rc = nf90mpi_def_var(ncido,'seaSurfacePressure',nf90mpi_DOUBLE, (/dimidc, dimidt/), varid)
      rc = nf90mpi_def_var(ncido,'ssh',nf90mpi_DOUBLE, (/dimidc, dimidt/), varid)
   !rc = nf90mpi_def_var(ncido,'temperature',nf90mpi_DOUBLE, (/dimidl,dimidc,dimidt/), varid)
      rc = nf90mpi_def_var(ncido,'normalVelocity',nf90mpi_DOUBLE, (/dimidl, dimide, dimidt/), varid)
      rc = nf90mpi_def_var(ncido,'oceanFracObserved',nf90mpi_DOUBLE, (/dimidc/), varid)
      rc = nf90mpi_def_var(ncido,'refSSH',nf90mpi_DOUBLE, (/dimidc, dimidt/), varid)
      rc = nf90mpi_def_var(ncido,'deltaSSH',nf90mpi_DOUBLE, (/dimidc, dimidt/), varid)
      rc = nf90mpi_def_var(ncido,'density',nf90mpi_DOUBLE, (/dimidl, dimidc, dimidt/), varid)
      rc = nf90mpi_def_var(ncido,'boundaryLayerDepth',nf90mpi_DOUBLE, (/dimidc, dimidt/), varid)
      rc = nf90mpi_def_var(ncido,'modifySSHMask',nf90mpi_INT, (/dimidc, dimidt/), varid)
      rc = nf90mpi_def_var(ncido,'refZMid',nf90mpi_DOUBLE, (/dimidl/), varid)

   rc = nf90mpi_def_var(ncido,'minLevelCell',nf90mpi_INT, (/dimidc/), varid)
   rc = nf90mpi_def_var(ncido,'maxLevelCell',nf90mpi_INT, (/dimidc/), varid)
   rc = nf90mpi_def_var(ncido,'restingThickness',nf90mpi_DOUBLE, (/dimidl,dimidc/), varid)
   rc = nf90mpi_def_var(ncido,'layerThickness',nf90mpi_DOUBLE, (/dimidl,dimidc,dimidt/), varid)
   rc = nf90mpi_def_var(ncido,'fCell',nf90mpi_DOUBLE, (/dimidc/), varid)
   rc = nf90mpi_def_var(ncido,'fVertex',nf90mpi_DOUBLE, (/dimidv/), varid)
   rc = nf90mpi_def_var(ncido,'fEdge',nf90mpi_DOUBLE, (/dimide/), varid)
   rc = nf90mpi_def_var(ncido,'refBottomDepth',nf90mpi_DOUBLE, (/dimidl/), varid)
   rc = nf90mpi_def_var(ncido,'vertCoordMovementWeights',nf90mpi_DOUBLE, (/dimidl/), varid)

   rc = nf90_get_att(ncida,NF90_GLOBAL,'sphere_radius',sphere_radius)
   if(sphere_radius < 10.) then
      ! unit sphere - reset sphere_radius to real earth
      sphere_radius = 6371229.
      rc = nf90mpi_put_att(ncido, NF90mpi_GLOBAL, 'sphere_radius', sphere_radius)
   else
      ! this is the real earth - write then reset to 1 so grid metrics don't change
      rc = nf90mpi_put_att(ncido, NF90mpi_GLOBAL, 'sphere_radius', sphere_radius)
      sphere_radius= 1.0
   endif

  
   rc = nf90mpi_enddef(ncido)

   ! redefine cellsoncell
   !int cellsOnCell(nCells, maxEdges) ;
   !rc = nf90_redef(ncido)
   !rc = nf90_def_var(ncido,'cellsOnCell',NF90_INT, (/dimidme,dimidc/), varid)
   rc = nf90mpi_inq_varid(ncido,'cellsOnCell',varid)
   !rc = nf90_enddef(ncido)
   do n = 1,ncellso
      do k = 1,6
         integer1d(k) = indexcella2o(cellsoncella(k,indexcello2a(n)))
         integer2d(k,n) = integer1d(k)
      enddo
   enddo
      rc = nf90mpi_put_var_all(ncido,varid,integer2d(1:6,:))
   
   allocate(cellsonvertex(3,nverticeso))

!!!

   open(unit=7,file='refBottomDepth',form='formatted')
      read(7,*)refbottomdepth
   close(7)
   
   ! get the bottom depth, then the layer thicknesses
   call compress_1d_r8(ncellsa, ncellso, bathy, bottomdepth, cellexistsa)
   bottomdepth = -bottomdepth
   bottomdepth = min(refBottomDepth(100),bottomdepth)
   bottomdepth = max(refBottomDepth(8),bottomdepth)

!   ! for now bottomdepth is rounded to the nearest refbottomdepth value
!   do n = 1,ncellso
!      kk = 1
!      do k = 2,100
!        if(bottomdepth(n) > 0.5*(refbottomdepth(k)+refbottomdepth(k-1)) )then
!          kk = k
!        endif
!      enddo
!      bottomdepth(n) = refbottomdepth(kk)
!   enddo

   refzmid(1) = 0.5 * refbottomdepth(1)
   do n = 2,100
      refzmid(n) = 0.5*(refbottomdepth(n)+refbottomdepth(n-1))
   enddo
   refzmid = -refzmid

   rc = nf90mpi_inq_varid(ncido,'refZMid',varid)
   rc = nf90mpi_put_var_all(ncido,varid,refZMid)
   
   do n = 1,ncellso
      restingthickness(1,n) = refbottomdepth(1)
      do k = 2,100
         wrk1 = refbottomdepth(k) - refbottomdepth(k-1)
         wrk2 = bottomdepth(n)-refbottomdepth(k-1)
         restingthickness(k,n) = min(wrk1, wrk2)
         restingthickness(k,n) = max(restingthickness(k,n), 0.1*wrk1)
         if(wrk2<0.0)restingthickness(k,n) = 0.0
      enddo
   enddo

   do n = 1,ncellso
      bottomdepth(n) = sum(restingthickness(:,n))
   enddo
   rc = nf90mpi_inq_varid(ncido,'bottomDepth',varid)
   rc = nf90mpi_put_var_all(ncido,varid,bottomdepth)


!!!
! define and write other variables

   !lonCell
   rc = nf90_inq_varid(ncida,'lonCell',varid)
   rc = nf90_get_var  (ncida,varid,real8_1d_a)
   call compress_1d_r8(ncellsa, ncellso, real8_1d_a, real8_1d, cellexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'lonCell',NF90_DOUBLE, (/dimidc/), varid)
   rc = nf90mpi_inq_varid(ncido,'lonCell',varid)
   !   rc = nf90_enddef(ncido)
   rc = nf90mpi_put_var_all(ncido,varid,real8_1d)
   loncello = real8_1d * 180./pi
   where (loncello < 0.) loncello = loncello + 360.
   where (loncello > 360.) loncello = loncello - 360.

   !latCell
   rc = nf90_inq_varid(ncida,'latCell',varid)
   rc = nf90_get_var  (ncida,varid,real8_1d_a)
   call compress_1d_r8(ncellsa, ncellso, real8_1d_a, real8_1d, cellexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'latCell',NF90_DOUBLE, (/dimidc/), varid)
   rc = nf90mpi_inq_varid(ncido,'latCell',varid)
   !   rc = nf90_enddef(ncido)
   rc = nf90mpi_put_var_all(ncido,varid,real8_1d)
   latcello = real8_1d * 180./pi


   !double areaCell(nCells) ;
   rc = nf90_inq_varid(ncida,'areaCell',varid)
   rc = nf90_get_var  (ncida,varid,real8_1d_a)

   call compress_1d_r8(ncellsa, ncellso, real8_1d_a, real8_1d, cellexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'areaCell',NF90_DOUBLE, (/dimidc/), varid)
   rc = nf90mpi_inq_varid(ncido,'areaCell',varid)
   !   rc = nf90_enddef(ncido)
   real8_1d = real8_1d * sphere_radius**2
   rc = nf90mpi_put_var_all(ncido,varid,real8_1d)


   
   !double angleEdge(nEdges) ;
   rc = nf90_inq_varid(ncida,'angleEdge',varid)
   rc = nf90_get_var  (ncida,varid,real8_edge_a)
   call compress_1d_r8(nedgesa, nedgeso, real8_edge_a, angleedge, edgeexistsa)

   !int cellsOnEdge(nEdges, TWO) ;
   allocate(cellsonedge (2,nedgeso))
   !rc = nf90_redef(ncido)
   !rc = nf90_def_var(ncido,'cellsOnEdge',NF90_INT, (/dimid2,dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'cellsOnEdge',varid)
   !rc = nf90_enddef(ncido)
   do k = 1,2
      do n = 1,nedgeso
         cellsonedge(k,n) = indexcella2o(cellsonedgea(k,indexedgeo2a(n)))
      enddo
   enddo
!modify angle edge too
   angleswap = .false.
   do n = 1,nedgeso
      if(cellsonedge(1,n) > cellsonedge(2,n))  then
         iii = cellsonedge(1,n)
         cellsonedge(1,n) = cellsonedge(2,n)
         cellsonedge(2,n) = iii
         angleedge(n) = angleedge(n) - sign(pi,angleedge(n))
         angleswap(n) = .true.
      endif
      if(cellsonedge(1,n) == 0)  then
         cellsonedge(1,n) = cellsonedge(2,n)
         cellsonedge(2,n) = 0
         angleswap(n) = .true.
         angleedge(n) = angleedge(n) - sign(pi,angleedge(n))
      endif
   enddo
   rc = nf90mpi_put_var_all(ncido,varid,cellsonedge)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'angleEdge',NF90_DOUBLE, (/dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'angleEdge',varid)
   !   rc = nf90_enddef(ncido)
   rc = nf90mpi_put_var_all(ncido,varid,angleedge)

   !int cellsOnVertex(nVertices, vertexDegree) ;
   rc = nf90_inq_varid(ncida,'cellsOnVertex',varid2)
!   rc = nf90_redef(ncido)
!   rc = nf90_def_var(ncido,'cellsOnVertex',NF90_INT, (/dimidvd,dimidv/), varid)
   rc = nf90mpi_inq_varid(ncido,'cellsOnVertex',varid)
!   rc = nf90_enddef(ncido)
   do n = 1,nverticeso
      rc = nf90_get_var  (ncida,varid2,integer1d_vert_a(1:3),(/1,indexverto2a(n)/),(/3,1/))
      do k = 1,3
         cellsonvertex(k,n) = indexcella2o(integer1d_vert_a(k))
      enddo
   enddo
      rc = nf90mpi_put_var_all(ncido,varid,cellsonvertex)

   !double dcEdge(nEdges) ;
   rc = nf90_inq_varid(ncida,'dcEdge',varid)
print*,'inqvarid ',rc
   rc = nf90_get_var  (ncida,varid,real8_edge_a)
print*,'getvar ',rc
   call compress_1d_r8(nedgesa, nedgeso, real8_edge_a, real8_edge, edgeexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'dcEdge',NF90_DOUBLE, (/dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'dcEdge',varid)
print*,'inqvarid ',rc
   !   rc = nf90_enddef(ncido)
   real8_edge = real8_edge * sphere_radius
   rc = nf90mpi_put_var_all(ncido,varid,real8_edge)
print*,'putvar ',rc

   !double dvEdge(nEdges) ;
   rc = nf90_inq_varid(ncida,'dvEdge',varid)
   rc = nf90_get_var  (ncida,varid,real8_edge_a)
   call compress_1d_r8(nedgesa, nedgeso, real8_edge_a, real8_edge, edgeexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'dvEdge',NF90_DOUBLE, (/dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'dvEdge',varid)
   real8_edge = real8_edge * sphere_radius
   rc = nf90mpi_put_var_all(ncido,varid,real8_edge)

allocate(intwork2d(6,ncellso))
   !int edgesOnCell(nCells, maxEdges) ;
   rc = nf90_inq_varid(ncida,'edgesOnCell',varid2)
   !rc = nf90_redef(ncido)
   !rc = nf90_def_var(ncido,'edgesOnCell',NF90_INT, (/dimidme,dimidc/), varid)
   rc = nf90mpi_inq_varid(ncido,'edgesOnCell',varid)
   !rc = nf90_enddef(ncido)
   do k = 1,6
      do n = 1,ncellsa,10000
      ne = min(ncellsa,n+9999)
      nlen = min(10000,ncellsa+1-n)
      rc = nf90_get_var  (ncida,varid2,integer1d_a(n:ne),(/k,n/),(/1,nlen/))
      enddo
      integer1d = 1
      do n = 1,ncellso
         integer1d(n) = indexedgea2o(integer1d_a(indexcello2a(n)))
      intwork2d(k,n) = integer1d(n)
      enddo
   enddo
      rc = nf90mpi_put_var_all(ncido,varid,intwork2d(1:6,:))
deallocate(intwork2d)

allocate(intwork2d(12,nedgeso))
   !int edgesOnEdge(nEdges, maxEdges2) ;
   rc = nf90_inq_varid(ncida,'edgesOnEdge',varid2)
   !rc = nf90_redef(ncido)
   !rc = nf90_def_var(ncido,'edgesOnEdge',NF90_INT, (/dimidme2,dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'edgesOnEdge',varid)
   !rc = nf90_enddef(ncido)
   do k = 1,12
      do n = 1,nedgesa,10000
      ne = min(nedgesa,n+9999)
      nlen = min(10000,nedgesa+1-n)
      rc = nf90_get_var  (ncida,varid2,integer1d_edge_a(n:ne),(/k,n/),(/1,nlen/))
      enddo
      do n = 1,nedgeso
         integer1d_edge(n) = indexedgea2o(integer1d_edge_a(indexedgeo2a(n)))
      intwork2d(k,n) = integer1d_edge(n)
      enddo
   enddo
      rc = nf90mpi_put_var_all(ncido,varid,intwork2d(1:12,:))
deallocate(intwork2d)
   
allocate(intwork2d(3,nverticeso))
   !int edgesOnVertex(nVertices, vertexDegree) ;
   rc = nf90_inq_varid(ncida,'edgesOnVertex',varid2)
   !rc = nf90_redef(ncido)
   !rc = nf90_def_var(ncido,'edgesOnVertex',NF90_INT, (/dimidvd,dimidv/), varid)
   rc = nf90mpi_inq_varid(ncido,'edgesOnVertex',varid)
   !rc = nf90_enddef(ncido)
   do k = 1,3
      do n = 1,nverticesa,10000
      ne = min(nverticesa,n+9999)
      nlen = min(10000,nverticesa+1-n)
      rc = nf90_get_var  (ncida,varid2,integer1d_vert_a(n:ne),(/k,n/),(/1,nlen/))
      enddo
      do n = 1,nverticeso
         integer1d_vert(n) = indexedgea2o(integer1d_vert_a(indexverto2a(n)))
      intwork2d(k,n) = integer1d_vert(n)
      enddo
   enddo
      rc = nf90mpi_put_var_all(ncido,varid,intwork2d(1:3,:))
deallocate(intwork2d)


   !int indexToCellID(nCells) ;
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'indexToCellID',NF90_INT, (/dimidc/), varid)
   rc = nf90mpi_inq_varid(ncido,'indexToCellID',varid)
   !   rc = nf90_enddef(ncido)
      do n = 1,ncellso
         integer1d(n) = n
      enddo
      rc = nf90mpi_put_var_all(ncido,varid,integer1d)

   !int indexToEdgeID(nEdges) ;
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'indexToEdgeID',NF90_INT, (/dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'indexToEdgeID',varid)
   !   rc = nf90_enddef(ncido)
      do n = 1,nedgeso
         integer1d_edge(n) = n
      enddo
      rc = nf90mpi_put_var_all(ncido,varid,integer1d_edge)

   !int indexToVertexID(nVertices) ;
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'indexToVertexID',NF90_INT, (/dimidv/), varid)
   rc = nf90mpi_inq_varid(ncido,'indexToVertexID',varid)
   !   rc = nf90_enddef(ncido)
      do n = 1,nverticeso
         integer1d_vert(n) = n
      enddo
      rc = nf90mpi_put_var_all(ncido,varid,integer1d_vert)

      allocate(r8work2d(3,nverticeso))
   !double kiteAreasOnVertex(nVertices, vertexDegree) ;
   rc = nf90_inq_varid(ncida,'kiteAreasOnVertex',varid2)
   !rc = nf90_redef(ncido)
   !rc = nf90_def_var(ncido,'kiteAreasOnVertex',NF90_DOUBLE, (/dimidvd,dimidv/), varid)
   rc = nf90mpi_inq_varid(ncido,'kiteAreasOnVertex',varid)
   !rc = nf90_enddef(ncido)
   areatriangle = 0.0
   do k = 1,3
      rc = nf90_get_var  (ncida,varid2,real8_vert_a,(/k,1/),(/1,nverticesa/))
      do n = 1,nverticeso
         real8_vert(n) = real8_vert_a(indexverto2a(n))
         ! zero out the kites that belong to non-existent cellsoncells
         if(cellsonvertex(k,n) == 0) real8_vert(n) = 0.0
         areatriangle(n) = areatriangle(n) + real8_vert(n)
      enddo
      r8work2d(k,:) = real8_vert(:)
   enddo
   r8work2d = r8work2d * sphere_radius**2
      rc = nf90mpi_put_var_all(ncido,varid,r8work2d)
      deallocate(r8work2d)
   !double areaTriangle(nVertices) ;
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'areaTriangle',NF90_DOUBLE, (/dimidv/), varid)
   rc = nf90mpi_inq_varid(ncido,'areaTriangle',varid)
   !   rc = nf90_enddef(ncido)
   areatriangle = areatriangle * sphere_radius**2
   rc = nf90mpi_put_var_all(ncido,varid,areatriangle)
      
   !double latEdge(nEdges) ;
   rc = nf90_inq_varid(ncida,'latEdge',varid)
   rc = nf90_get_var  (ncida,varid,real8_edge_a)
   call compress_1d_r8(nedgesa, nedgeso, real8_edge_a, real8_edge, edgeexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'latEdge',NF90_DOUBLE, (/dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'latEdge',varid)
   !   rc = nf90_enddef(ncido)
   rc = nf90mpi_put_var_all(ncido,varid,real8_edge)

   !double latVertex(nVertices) ;
   rc = nf90_inq_varid(ncida,'latVertex',varid)
   rc = nf90_get_var  (ncida,varid,real8_vert_a)
   call compress_1d_r8(nverticesa, nverticeso, real8_vert_a, real8_vert, vertexexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'latVertex',NF90_DOUBLE, (/dimidv/), varid)
   rc = nf90mpi_inq_varid(ncido,'latVertex',varid)
   !   rc = nf90_enddef(ncido)
   rc = nf90mpi_put_var_all(ncido,varid,real8_vert)

   !double lonEdge(nEdges) ;
   rc = nf90_inq_varid(ncida,'lonEdge',varid)
   rc = nf90_get_var  (ncida,varid,real8_edge_a)
   call compress_1d_r8(nedgesa, nedgeso, real8_edge_a, real8_edge, edgeexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'lonEdge',NF90_DOUBLE, (/dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'lonEdge',varid)
   !   rc = nf90_enddef(ncido)
   rc = nf90mpi_put_var_all(ncido,varid,real8_edge)

   !double lonVertex(nVertices) ;
   rc = nf90_inq_varid(ncida,'lonVertex',varid)
   rc = nf90_get_var  (ncida,varid,real8_vert_a)
   call compress_1d_r8(nverticesa, nverticeso, real8_vert_a, real8_vert, vertexexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'lonVertex',NF90_DOUBLE, (/dimidv/), varid)
   rc = nf90mpi_inq_varid(ncido,'lonVertex',varid)
   !   rc = nf90_enddef(ncido)
   rc = nf90mpi_put_var_all(ncido,varid,real8_vert)

   !double meshDensity(nCells) ;
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'meshDensity',NF90_DOUBLE, (/dimidc/), varid)
   rc = nf90mpi_inq_varid(ncido,'meshDensity',varid)
   !   rc = nf90_enddef(ncido)
      real8_1d = 1.0
   rc = nf90mpi_put_var_all(ncido,varid,real8_1d)
   
   !int nEdgesOnCell(nCells) ;
   rc = nf90_inq_varid(ncida,'nEdgesOnCell',varid)
      do n = 1,ncellsa,10000
      ne = min(ncellsa,n+9999)
      nlen = min(10000,ncellsa+1-n)
      rc = nf90_get_var(ncida,varid,integer1d_a(n:ne),(/n/),(/nlen/))
      enddo
   call compress_1d_int(ncellsa, ncellso, integer1d_a, integer1d, cellexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'nEdgesOnCell',NF90_INT, (/dimidc/), varid)
   rc = nf90mpi_inq_varid(ncido,'nEdgesOnCell',varid)
      rc = nf90mpi_put_var_all(ncido,varid,integer1d)

   !int nEdgesOnEdge(nEdges) ;
   rc = nf90_inq_varid(ncida,'nEdgesOnEdge',varid)
      do n = 1,nedgesa,10000
      ne = min(nedgesa,n+9999)
      nlen = min(10000,nedgesa+1-n)
      rc = nf90_get_var(ncida,varid,integer1d_edge_a(n:ne),(/n/),(/nlen/))
      enddo
   call compress_1d_int(nedgesa, nedgeso, integer1d_edge_a, integer1d_edge, edgeexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'nEdgesOnEdge',NF90_INT, (/dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'nEdgesOnEdge',varid)
      rc = nf90mpi_put_var_all(ncido,varid,integer1d_edge)

      allocate(intwork2d(6,ncellso))
   !int verticesOnCell(nCells, maxEdges) ;
   rc = nf90_inq_varid(ncida,'verticesOnCell',varid2)
   !rc = nf90_redef(ncido)
   !rc = nf90_def_var(ncido,'verticesOnCell',NF90_INT, (/dimidme,dimidc/), varid)
   rc = nf90mpi_inq_varid(ncido,'verticesOnCell',varid)
   !rc = nf90_enddef(ncido)
   do k = 1,6
      do n = 1,ncellsa,10000
      ne = min(ncellsa,n+9999)
      nlen = min(10000,ncellsa+1-n)
      rc = nf90_get_var  (ncida,varid2,integer1d_a(n:ne),(/k,n/),(/1,nlen/))
      enddo
      do n = 1,ncellso
         integer1d(n) = indexverta2o(integer1d_a(indexcello2a(n)))
      enddo
      intwork2d(k,:) = integer1d
print*,'putvar k',k,rc
   enddo
   print*,'verticesOnCell',intwork2d(:,1)
      rc = nf90mpi_put_var_all  (ncido,varid,intwork2d)
      deallocate(intwork2d)

      allocate(intwork2d(2,nedgeso))
   !int verticesOnEdge(nEdges, TWO) ;
   rc = nf90_inq_varid(ncida,'verticesOnEdge',varid2)
   !rc = nf90_redef(ncido)
   !rc = nf90_def_var(ncido,'verticesOnEdge',NF90_INT, (/dimid2,dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'verticesOnEdge',varid)
   !rc = nf90_enddef(ncido)
   do k = 1,2
      do n = 1,nedgesa,10000
      ne = min(nedgesa,n+9999)
      nlen = min(10000,nedgesa+1-n)
      rc = nf90_get_var  (ncida,varid2,integer1d_edge_a(n:ne),(/k,n/),(/1,nlen/))
      enddo
      do n = 1,nedgeso
         integer1d_edge(n) = indexverta2o(integer1d_edge_a(indexedgeo2a(n)))
      enddo
      intwork2d(k,:) = integer1d_edge
   enddo
   print*,'verticesOnedge',intwork2d(:,1)
      rc = nf90mpi_put_var_all  (ncido,varid,intwork2d)
      deallocate(intwork2d)

      allocate(r8work2d(12,nedgeso))
   !double weightsOnEdge(nEdges, maxEdges2) ;
   rc = nf90_inq_varid(ncida,'edgesOnEdge',varid3)
   rc = nf90_inq_varid(ncida,'weightsOnEdge',varid2)
   !rc = nf90_redef(ncido)
   !rc = nf90_def_var(ncido,'weightsOnEdge',NF90_DOUBLE, (/dimidme2,dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'weightsOnEdge',varid)
   !rc = nf90_enddef(ncido)
   do k = 1,12
      do n = 1,nedgesa,10000
      ne = min(nedgesa,n+9999)
      nlen = min(10000,nedgesa+1-n)
      rc = nf90_get_var  (ncida,varid2,real8_edge_a(n:ne),(/k,n/),(/1,nlen/))
      rc = nf90_get_var  (ncida,varid3,integer1d_edge_a(n:ne),(/k,n/),(/1,nlen/))
      enddo
      real8_edge = 0.0
      do n = 1,nedgeso
         real8_edge(n) = real8_edge_a(indexedgeo2a(n))
         m = (integer1d_edge_a(indexedgeo2a(n)))
         if(m>0)then 
            !if(indexedgea2o(m)>0) then
            !   if(angleswap(indexedgea2o(m)))real8_edge(n) = -real8_edge(n)
            !endif
            if(angleswap(n))real8_edge(n) = -real8_edge(n)
         endif
      enddo
      r8work2d(k,:) = real8_edge
   enddo
      rc = nf90mpi_put_var_all  (ncido,varid,r8work2d)
      deallocate(r8work2d)
      
   !double xCell(nCells) ;
   rc = nf90_inq_varid(ncida,'xCell',varid)
   rc = nf90_get_var  (ncida,varid,real8_1d_a)
   call compress_1d_r8(ncellsa, ncellso, real8_1d_a, real8_1d, cellexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'xCell',NF90_DOUBLE, (/dimidc/), varid)
   rc = nf90mpi_inq_varid(ncido,'xCell',varid)
   !   rc = nf90_enddef(ncido)
   real8_1d = real8_1d * sphere_radius
   rc = nf90mpi_put_var_all(ncido,varid,real8_1d)

   !double xEdge(nEdges) ;
   rc = nf90_inq_varid(ncida,'xEdge',varid)
   rc = nf90_get_var  (ncida,varid,real8_edge_a)
   call compress_1d_r8(nedgesa, nedgeso, real8_edge_a, real8_edge, edgeexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'xEdge',NF90_DOUBLE, (/dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'xEdge',varid)
   !   rc = nf90_enddef(ncido)
   real8_edge = real8_edge * sphere_radius
   rc = nf90mpi_put_var_all(ncido,varid,real8_edge)

   !double xVertex(nVertices) ;
   rc = nf90_inq_varid(ncida,'xVertex',varid)
   rc = nf90_get_var  (ncida,varid,real8_vert_a)
   call compress_1d_r8(nverticesa, nverticeso, real8_vert_a, real8_vert, vertexexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'xVertex',NF90_DOUBLE, (/dimidv/), varid)
   rc = nf90mpi_inq_varid(ncido,'xVertex',varid)
   !   rc = nf90_enddef(ncido)
   real8_vert = real8_vert * sphere_radius
   rc = nf90mpi_put_var_all(ncido,varid,real8_vert)

   !double yCell(nCells) ;
   rc = nf90_inq_varid(ncida,'yCell',varid)
   rc = nf90_get_var  (ncida,varid,real8_1d_a)
   call compress_1d_r8(ncellsa, ncellso, real8_1d_a, real8_1d, cellexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'yCell',NF90_DOUBLE, (/dimidc/), varid)
   rc = nf90mpi_inq_varid(ncido,'yCell',varid)
   !   rc = nf90_enddef(ncido)
   real8_1d = real8_1d * sphere_radius
   rc = nf90mpi_put_var_all(ncido,varid,real8_1d)

   !double yEdge(nEdges) ;
   rc = nf90_inq_varid(ncida,'yEdge',varid)
   rc = nf90_get_var  (ncida,varid,real8_edge_a)
   call compress_1d_r8(nedgesa, nedgeso, real8_edge_a, real8_edge, edgeexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'yEdge',NF90_DOUBLE, (/dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'yEdge',varid)
   !   rc = nf90_enddef(ncido)
   real8_edge = real8_edge * sphere_radius
   rc = nf90mpi_put_var_all(ncido,varid,real8_edge)

   !double yVertex(nVertices) ;
   rc = nf90_inq_varid(ncida,'yVertex',varid)
   rc = nf90_get_var  (ncida,varid,real8_vert_a)
   call compress_1d_r8(nverticesa, nverticeso, real8_vert_a, real8_vert, vertexexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'yVertex',NF90_DOUBLE, (/dimidv/), varid)
   rc = nf90mpi_inq_varid(ncido,'yVertex',varid)
   !   rc = nf90_enddef(ncido)
   real8_vert = real8_vert * sphere_radius
   rc = nf90mpi_put_var_all(ncido,varid,real8_vert)

   !double zCell(nCells) ;
   rc = nf90_inq_varid(ncida,'zCell',varid)
   rc = nf90_get_var  (ncida,varid,real8_1d_a)
   call compress_1d_r8(ncellsa, ncellso, real8_1d_a, real8_1d, cellexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'zCell',NF90_DOUBLE, (/dimidc/), varid)
   rc = nf90mpi_inq_varid(ncido,'zCell',varid)
   !   rc = nf90_enddef(ncido)
   real8_1d = real8_1d * sphere_radius
   rc = nf90mpi_put_var_all(ncido,varid,real8_1d)

   !double zEdge(nEdges) ;
   rc = nf90_inq_varid(ncida,'zEdge',varid)
   rc = nf90_get_var  (ncida,varid,real8_edge_a)
   call compress_1d_r8(nedgesa, nedgeso, real8_edge_a, real8_edge, edgeexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'zEdge',NF90_DOUBLE, (/dimide/), varid)
   rc = nf90mpi_inq_varid(ncido,'zEdge',varid)
   !   rc = nf90_enddef(ncido)
   real8_edge = real8_edge * sphere_radius
   rc = nf90mpi_put_var_all(ncido,varid,real8_edge)

   !double zVertex(nVertices) ;
   rc = nf90_inq_varid(ncida,'zVertex',varid)
   rc = nf90_get_var  (ncida,varid,real8_vert_a)
   call compress_1d_r8(nverticesa, nverticeso, real8_vert_a, real8_vert, vertexexistsa)
   !rc = nf90_redef(ncido)
   !   rc = nf90_def_var(ncido,'zVertex',NF90_DOUBLE, (/dimidv/), varid)
   rc = nf90mpi_inq_varid(ncido,'zVertex',varid)
   !   rc = nf90_enddef(ncido)
   real8_vert = real8_vert * sphere_radius
   rc = nf90mpi_put_var_all(ncido,varid,real8_vert)
   rc = nf90mpi_sync(ncido)

!print*,'salinity'   
!   !   allocate(r8work2d(100,ncellso))
!   !double salinity(Time, nCells, nVertLevels) ;
!   rc = nf90mpi_inq_varid(ncido,'salinity',varid)
!   do k = 1,100
!      real8_1d = 0.0
!      do n = 1,ncellso
!         if(k==1) then
!            real8_1d(n) = get_ts(loncello(n),latcello(n),k,indxs)
!         else
!            if(bottomDepth(n) > refbottomdepth(k-1) ) &
!               real8_1d(n) = get_ts(loncello(n),latcello(n),k,indxs)
!         endif
!      enddo
!       rc = nf90mpi_put_var_all(ncido,varid,real8_1d, (/k,1,1/)*1_mpi_offset_kind, (/1,ncellso,1/)*1_mpi_offset_kind)
!   !   r8work2d(k,:) = real8_1d(:)
!   enddo
!   !   rc = nf90mpi_put_var_all(ncido,varid,r8work2d)
!   rc = nf90mpi_sync(ncido)

!print*,'temperature'   
!   !double temperature(Time, nCells, nVertLevels) ;
!   rc = nf90mpi_inq_varid(ncido,'temperature',varid)
!   do k = 1,100
!      real8_1d = 0.0
!      do n = 1,ncellso
!         if(k==1) then
!            real8_1d(n) = get_ts(loncello(n),latcello(n),k,indxt)
!         else
!            if(bottomDepth(n) > refbottomdepth(k-1) ) &
!               real8_1d(n) = get_ts(loncello(n),latcello(n),k,indxt)
!         endif
!      enddo
!       rc = nf90mpi_put_var_all(ncido,varid,real8_1d, (/k,1,1/)*1_mpi_offset_kind, (/1,ncellso,1/)*1_mpi_offset_kind)
!     !r8work2d(k,:) = real8_1d(:)
!   enddo
!   rc = nf90mpi_sync(ncido)
   !   rc = nf90mpi_put_var_all(ncido,varid,r8work2d)
   !   deallocate(r8work2d)

   !double seaSurfacePressure(Time, nCells) ;
   rc = nf90mpi_inq_varid(ncido,'seaSurfacePressure',varid)
      real8_1d = 0.0
      rc = nf90mpi_put_var_all(ncido,varid,real8_1d)
   
   !double ssh(Time, nCells) ;
   rc = nf90mpi_inq_varid(ncido,'ssh',varid)
      real8_1d = 0.0
      rc = nf90mpi_put_var_all(ncido,varid,real8_1d)

   !double refSSH(Time, nCells) ;
   rc = nf90mpi_inq_varid(ncido,'refSSH',varid)
      real8_1d = 0.0
      rc = nf90mpi_put_var_all(ncido,varid,real8_1d)

   !double oceanFracObserved(Time, nCells) ;
   rc = nf90mpi_inq_varid(ncido,'oceanFracObserved',varid)
      real8_1d = 0.0
      rc = nf90mpi_put_var_all(ncido,varid,real8_1d)
   rc = nf90mpi_sync(ncido)

print*,'normalVelocity'   
   !double normalVelocity(Time, nEdges, nVertLevels) ;
   rc = nf90mpi_inq_varid(ncido,'normalVelocity',varid)
   real8_edge = 0.0
   do nb =1,nedgeso,10000
      ne = min(nb+9999,nedgeso)
      !print*,'nb,ne',nb,ne,(ne+1-nb)
      rc = nf90mpi_put_var_all(ncido,varid,real8_edge(1:100*(ne+1-nb)),(/1,nb,1/)*1_mpi_offset_kind,(/100,ne+1-nb,1/)*1_mpi_offset_kind)
   enddo
   rc = nf90mpi_sync(ncido)
   
   !int modifySSHMask(Time, nCells) ;
   rc = nf90mpi_inq_varid(ncido,'modifySSHMask',varid)
      integer1d = 0
      rc = nf90mpi_put_var_all(ncido,varid,integer1d)

   !double deltaSSH(Time, nCells) ;
   rc = nf90mpi_inq_varid(ncido,'deltaSSH',varid)
      real8_1d = 0.0
   rc = nf90mpi_put_var_all(ncido,varid,real8_1d)
   
   !double density(Time, nCells, nVertLevels) ;
   rc = nf90mpi_inq_varid(ncido,'density',varid)
   real8_1d = 0.0
   do k = 1,100
      rc = nf90mpi_put_var_all(ncido,varid,real8_1d, (/k,1/)*1_mpi_offset_kind, (/1,nedgeso/)*1_mpi_offset_kind)
   enddo
   
   !double boundaryLayerDepth(Time, nCells) ;
   rc = nf90mpi_inq_varid(ncido,'boundaryLayerDepth',varid)
      real8_1d = 0.0
   rc = nf90mpi_put_var_all(ncido,varid,real8_1d)
   
print*,'fcell'   
   !double fCell(nCells) ;
   rc = nf90_inq_varid(ncida,'latCell',varid)
   rc = nf90_get_var  (ncida,varid,real8_1d_a)
   call compress_1d_r8(ncellsa, ncellso, real8_1d_a, real8_1d, cellexistsa)

   rc = nf90mpi_inq_varid(ncido,'fCell',varid)
   real8_1d = 0.0001458424 * sin(real8_1d)
   rc = nf90mpi_put_var_all(ncido,varid,real8_1d)
   
print*,'fVertex'   
      !double fVertex(nVertices) ;
   rc = nf90_inq_varid(ncida,'latVertex',varid)
   rc = nf90_get_var  (ncida,varid,real8_vert_a)
   call compress_1d_r8(nverticesa, nverticeso, real8_vert_a, real8_vert, vertexexistsa)

   rc = nf90mpi_inq_varid(ncido,'fVertex',varid)
   real8_vert = 0.0001458424 * sin(real8_vert)
   rc = nf90mpi_put_var_all(ncido,varid,real8_vert)
   
   
print*,'fEdge'   
   !double fEdge(nEdges) ;
   rc = nf90_inq_varid(ncida,'latEdge',varid)
   rc = nf90_get_var  (ncida,varid,real8_edge_a)
   call compress_1d_r8(nedgesa, nedgeso, real8_edge_a, real8_edge, edgeexistsa)

   rc = nf90mpi_inq_varid(ncido,'fEdge',varid)
   real8_edge = 0.0001458424 * sin(real8_edge)
   rc = nf90mpi_put_var_all(ncido,varid,real8_edge)
   rc = nf90mpi_sync(ncido)
   
print*,'vertCoordMovementWeights'   
   !double vertCoordMovementWeights(nVertLevels) ;
   rc = nf90mpi_inq_varid(ncido,'vertCoordMovementWeights',varid)
   real8_1d(1:100) = 0.0
   rc = nf90mpi_put_var_all(ncido,varid,real8_1d(1:100))
   
   !double refBottomDepth(nVertLevels) ;
   rc = nf90mpi_inq_varid(ncido,'refBottomDepth',varid)
   rc = nf90mpi_put_var_all(ncido,varid,refBottomDepth)
   rc = nf90mpi_sync(ncido)

print*,'restingThickness'   
   !double restingThickness(nCells, nVertLevels) ;
   rc = nf90mpi_inq_varid(ncido,'restingThickness',varid)
   do k = 1,100
      rc = nf90mpi_put_var_all(ncido,varid,restingThickness(k,:), (/k,1/)*1_mpi_offset_kind, (/1,ncellso/)*1_mpi_offset_kind)
   enddo
   rc = nf90mpi_sync(ncido)
      
print*,'restingThickness'   
   !double restingThickness(time, nCells, nVertLevels) ;
   rc = nf90mpi_inq_varid(ncido,'layerThickness',varid)
   do k = 1,100
      rc = nf90mpi_put_var_all(ncido,varid,restingThickness(k,:), (/k,1,1/)*1_mpi_offset_kind, (/1,ncellso,1/)*1_mpi_offset_kind)
   enddo
   rc = nf90mpi_sync(ncido)

 print*,'minLevelCell'   
  !double minLevelCell(nCells) ;
   rc = nf90mpi_inq_varid(ncido,'minLevelCell',varid)
   integer1d = 1
   rc = nf90mpi_put_var_all(ncido,varid,integer1d)
   rc = nf90mpi_sync(ncido)

   !double maxLevelCell(nCells) ;
   rc = nf90mpi_inq_varid(ncido,'maxLevelCell',varid)
   do n = 1,ncellso
      integer1d(n) = 0
      do k = 1,100
         if(restingthickness(k,n) > 0.0) integer1d(n) = k
      enddo
   enddo
   rc = nf90mpi_put_var_all(ncido,varid,integer1d)
   print*,'maxlevelcell',rc

      
!!!
! finalize

   rc = nf90_close(ncid)
   rc = nf90_close(ncida)
   rc = nf90mpi_close(ncido)

   call mpi_finalize(err)

contains
   
   subroutine compress_1d_r8(ncellsa, ncellso, xa, xo, existsa)
   integer, intent(in) :: ncellsa, ncellso
   real*8, intent(in), dimension(ncellsa) :: xa
   real*8, intent(out), dimension(ncellso) :: xo
   integer, intent(in), dimension(ncellsa) :: existsa
   
   integer no, na
   
   no = 0
   do na = 1,ncellsa
      if(existsa(na)>0) then
        no = no + 1
        xo(no) = xa(na)
      endif
   enddo
   end subroutine compress_1d_r8
   
   subroutine compress_1d_int(ncellsa, ncellso, xa, xo, existsa)
   integer, intent(in) :: ncellsa, ncellso
   integer, intent(in), dimension(ncellsa) :: xa
   integer, intent(out), dimension(ncellso) :: xo
   integer, intent(in), dimension(ncellsa) :: existsa
   
   integer no, na
   
   no = 0
   do na = 1,ncellsa
      if(existsa(na)>0) then
        no = no + 1
        xo(no) = xa(na)
      endif
   enddo
   end subroutine compress_1d_int
   
end program build_culled_file_from_bathy
