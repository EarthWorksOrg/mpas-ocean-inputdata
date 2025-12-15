program build_interpolated_velocities

! compile with mpi, hdf5, netcdf, and pnetcdf libraries
!
! this program works from restart of the coarser resolution and the interpolation weights file
!   and replaces the zero velocities in the finer resolution file with interpolated velocities.
! As the iterative solver can take a lot of time this program write to a work file and
!   can be restarted from it.
! setenv OMP_NUM_THREADS 16

! METHOD:
!   read coarse mesh velocities and compute coarse divergence and vorticity
!   interpolate divergence and vorticity to the fine mesh.
!   iterative jacobi solver compute velocity potential (chi) from divergence, and
!       stream function (psi) from vorticity on the fine mesh.
!   when converged, construct velocities from psi and chi and write to the
!       fine mesh restart

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
   character(len=128) :: workfile, terpfile, coarsefile, finefile

   ! netcdf-related variables
   !   'c' suffix is coarse grid
   !   'f' suffix is fine grid
   integer :: rc, ncid, ncidc, ncidf, ncidd
   integer :: dimid, varid, attid, dimidc, dimidv, dimide, &
              dimidt, dimidme, dimid2, dimidvd, dimidme2, dimidsl,   &
              dimid1, dimidcat, dimidily, dimidsly, varid2, dimidl, &
              dimidv1, dimidfg, dimid7, dimid4, dimidcc, dimidcv
   integer :: natts, maxiter, numiter, iter_glob
   real*8, dimension(101) :: resid_psi, resid_chi
   character(len=80) :: attname

   ! grid variables
   ! vertices
   integer :: ncellsc, ncellsf, maxlinks, nverticesc, nverticesf, nedgesc, nedgesf
   
   
   real*8 r8min, r8max, wrk1, wrk2, areatem, sphere_radius, rlim, rms_v
   integer :: ne, nlen
   
      
   ! loop indices
   integer :: n, na, m, id, k, iii, kk, indxo, l, i, icell, iedge, maxlevel
   
   integer  l_init
      
   ! iteration variables
   logical converged
   
   ! remapping arrays
   real*4, dimension(:,:), allocatable :: remap_wgts
   integer, dimension(:,:), allocatable :: remap_indices
   real*4, dimension(:,:), allocatable :: remap_wgtsv
   integer, dimension(:,:), allocatable :: remap_indicesv
   
   ! stream function and velocity potential, divergence and vorticity
   real*8, dimension(:), pointer :: chi_n, chi_f, psi_n, psi_f
   real*8, dimension(:), allocatable :: div_f, vort_f
   real*8, dimension(:,:), allocatable, target :: chi, psi
   integer :: indexn, indexf
   ! other solve variables
   real*8, dimension(:,:), allocatable :: lwghts_chi, lwghts_psi
   
   ! Initialization variables
   !coarse grid arrays
      integer, dimension(:), allocatable :: nedgesOnCellc, maxlevelcellc,       &
                minlevelvertextopc, minlevelcellc, maxLevelVertexBotc
      integer, dimension(:,:), allocatable :: edgesoncellc, cellsonedgec,       &
                cellsonvertexc, verticesoncellc, edgesOnVertexc, verticesonedgec
      real*8, dimension(:), allocatable :: AreaCellc, dvedgec, areatrianglec,   &
                invareacellc, invareatrianglec, dcEdgec
      real*8, dimension(:,:), allocatable :: edgeSignOnCellc, edgesignonvertexc

      real*8, dimension(:), allocatable :: vorticityc, velocityc, divergencec
   ! fine grid arrays
      integer, dimension(:), allocatable :: maxlevelcellf       
      integer, dimension(:,:), allocatable :: cellsonvertexf, edgesoncellf    
      real*8, dimension(:), allocatable :: AreaCellf, areatrianglef, dcedgef, dvedgef
      real*8, dimension(:,:), allocatable :: kiteareasf
   ! init work arrays
      integer, dimension(:), allocatable :: ncount, flag, assigned_cell
   
   ! solver mesh variables - and used in all phases of program
   integer, dimension(:), allocatable :: nedgesoncellf
   integer, dimension(:,:), allocatable :: cellsoncellf,  &
         verticesonedgef, edgesonvertexf
   
   ! mpi variables
   integer :: err, rank, nprocs, info
          integer(kind=MPI_OFFSET_KIND) start(3), count(3), dimtem
          integer(kind=MPI_OFFSET_KIND) malloc_size, sum_size
          
   real*8 :: d2r, pi
  
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
   maxlevel = 101  ! adjust this for debugging
!!!1
! preliminaries
   workfile = 'solver_work_file.nc'

   print*,'enter 1 to initialize solver, 2 to continue'
   read(5,*)l_init
             
   call MPI_Info_free(info, err)

   ! open the fine restart xfile
   print *,' enter the fine restart filename'
   read(5,'(a128)') finefile
   rc = nf90mpi_open(mpi_comm_world,trim(finefile), NF90_write,info,ncidf)
   print*,'open fine restart',rc
   if(rc .ne. 0)print*,trim(finefile)

   ! read the fine grid dimensions
   rc = nf90mpi_inq_dimid(ncidf,'nCells',dimid)
   rc = nf90mpi_inquire_dimension(ncidf,dimid,len=dimtem)
   ncellsf = dimtem
   rc = nf90mpi_inq_dimid(ncidf,'nVertices',dimid)
   rc = nf90mpi_inquire_dimension(ncidf,dimid, len=dimtem)
   nverticesf = dimtem
   rc = nf90mpi_inq_dimid(ncidf,'nEdges',dimid)
   rc = nf90mpi_inquire_dimension(ncidf,dimid, len=dimtem)
   nedgesf = dimtem
   print*,'fine grid dimensions',ncellsf,nverticesf,nedgesf

   ! allocate psi, chi, vort and div, and read grid metric variables for iteration
   allocate(psi(nverticesf,2))
   allocate(chi(ncellsf,2))
   allocate(vort_f(nverticesf))
   allocate(div_f(ncellsf))
   allocate(nEdgesOnCellf(ncellsf))
   rc = nf90mpi_inq_varid(ncidf,'nEdgesOnCell',varid)
   print*,'nedgesoncellf ',rc
   rc = nf90mpi_get_var_all(ncidf,varid,nedgesoncellf)
   print*,'nedgesoncellf ',rc, nedgesoncellf(1)
   allocate(cellsoncellf(6,ncellsf))
   rc = nf90mpi_inq_varid(ncidf,'cellsOnCell',varid)
   rc = nf90mpi_get_var_all(ncidf,varid,cellsoncellf)
   print*,'cellsoncellf',cellsoncellf(:,1651860)
   allocate(edgesonvertexf(3,nverticesf))
   rc = nf90mpi_inq_varid(ncidf,'edgesOnVertex',varid)
   rc = nf90mpi_get_var_all(ncidf,varid,edgesonvertexf)
   allocate(verticesonedgef(2,nedgesf))
   rc = nf90mpi_inq_varid(ncidf,'verticesOnEdge',varid)
   rc = nf90mpi_get_var_all(ncidf,varid,verticesonedgef)
   allocate(lwghts_chi(0:6,ncellsf))  ! zero index is the cell center, 1 to 6 is neighbors
   allocate(lwghts_psi(0:3,nverticesf))  ! zero index is the triangle center, 1 to 3 is neighbors
   
 !!!! begin init/restart block
     
   if(l_init == 1) then
      print *,' enter the coarse restart filename'
      read(5,'(a128)') coarsefile
      print *,' enter the interpolation weights filename'
      read(5,'(a128)') terpfile
      
      print*,'enter maximum permissible solution residual'
      read(5,*)rlim
      print*,'enter maximum permissible solution iterations'
      read(5,*)maxiter
      print*,'enter number of iterations between restarts'
      read(5,*)numiter
      print*,'enter number of levels to compute'
      read(5,*)maxlevel

      ! open the coarse restart file
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
   print*,'coarse grid dimensions',ncellsc,nverticesc,nedgesc
      
      iter_glob = 0
      resid_chi = 1.e9
      resid_psi = 1.e9
      ! create the coarse tem file
      rc = nf90_create('./coarse_tem_file.nc',nf90_clobber,ncidd)
      rc = nf90_def_dim(ncidd,'nCells',ncellsc,dimidcc)
      rc = nf90_def_dim(ncidd,'nVertices',nverticesc,dimidcv)
      rc = nf90_def_dim(ncidd,'nlev',101,dimid)
      rc = nf90_def_var(ncidd,'div',nf90_double,(/dimidcc,dimid/),varid=varid)
      rc = nf90_def_var(ncidd,'vort',nf90_double,(/dimidcv,dimid/),varid=varid)
      rc = nf90_enddef(ncidd)

      ! create the workfile
      rc = nf90_create(trim(workfile),ior(nf90_clobber,nf90_64bit_offset),ncid)
      rc = nf90_def_dim(ncid,'nCells',ncellsf,dimidc)
      rc = nf90_def_dim(ncid,'nVertices',nverticesf,dimidv)
      rc = nf90_def_dim(ncid,'nEdges',nedgesf,dimide)
      rc = nf90_def_dim(ncid,'nlev',101,dimid)
      rc = nf90_def_dim(ncid,'seven',7,dimid7)
      rc = nf90_def_dim(ncid,'four',4,dimid4)
      rc = nf90_def_var(ncid,'iter_glob',nf90_int,varid=varid)
      rc = nf90_def_var(ncid,'numiter',nf90_int,varid=varid)
      rc = nf90_def_var(ncid,'maxiter',nf90_int,varid=varid)
      rc = nf90_def_var(ncid,'rlim',nf90_double,varid=varid)
      rc = nf90_def_var(ncid,'resid_chi',nf90_double,(/dimid/),varid=varid)
      rc = nf90_def_var(ncid,'resid_psi',nf90_double,(/dimid/),varid=varid)
      rc = nf90_def_var(ncid,'lwghts_chi',nf90_double,(/dimid7,dimidc/),varid=varid)
      rc = nf90_def_var(ncid,'lwghts_psi',nf90_double,(/dimid4,dimidv/),varid=varid)
      rc = nf90_def_var(ncid,'chi',nf90_double,(/dimidc,dimid/),varid=varid)
      rc = nf90_def_var(ncid,'div',nf90_double,(/dimidc,dimid/),varid=varid)
      rc = nf90_def_var(ncid,'psi',nf90_double,(/dimidv,dimid/),varid=varid)
      rc = nf90_def_var(ncid,'vort',nf90_double,(/dimidv,dimid/),varid=varid)
      rc = nf90_enddef(ncid)
      
      rc = nf90_inq_varid(ncid,'iter_glob',varid)
      rc = nf90_put_var(ncid,varid,iter_glob)
      rc = nf90_inq_varid(ncid,'maxiter',varid)
      rc = nf90_put_var(ncid,varid,maxiter)
      rc = nf90_inq_varid(ncid,'numiter',varid)
      rc = nf90_put_var(ncid,varid,numiter)
      rc = nf90_inq_varid(ncid,'rlim',varid)
      rc = nf90_put_var(ncid,varid,rlim)
      rc = nf90_inq_varid(ncid,'resid_chi',varid)
      rc = nf90_put_var(ncid,varid,resid_chi)
      rc = nf90_inq_varid(ncid,'resid_psi',varid)
      rc = nf90_put_var(ncid,varid,resid_psi)

print*,'calling pre_init_solver'      
      call pre_init_solver()
print*,'exiting pre_init_solver'      
      
      ! initialize solvers
      allocate(velocityc(nedgesc))
      allocate(divergencec(ncellsc))
      allocate(vorticityc(nverticesc))
      !print *,'k loop'
      do k = 1,maxlevel
         if(k < 101) then
            rc = nf90mpi_inq_varid(ncidc,'normalVelocity',varid)
         !print*,'inqvarid velocityc',k,rc
            rc = nf90mpi_get_var_all(ncidc,varid,velocityc,        &
                  (/k,1,1/)*1_mpi_offset_kind,(/1,nedgesc,1/)*1_mpi_offset_kind)
         else
            rc = nf90mpi_inq_varid(ncidc,'normalBarotropicVelocity',varid)
         !print*,'inqvarid velocityc',k,rc
            rc = nf90mpi_get_var_all(ncidc,varid,velocityc)
         endif
         !print*,'read velocityc',k,rc
         
         rms_v = 0.0
         do i = 1,nedgesc
            rms_v = rms_v + velocityc(i)**2
         enddo
         print*,'rms_vel coarse ',k,rms_v/nedgesc
         call initialize_chi_solver(k)
         call initialize_psi_solver(k)
         call adjust_vort_div_means(k)
                  
         rc = nf90_inq_varid(ncid,'vort',varid)
         rc = nf90_put_var(ncid,varid,vort_f,(/1,k/),(/nverticesf,1/))
         rc = nf90_inq_varid(ncid,'div',varid)
         rc = nf90_put_var(ncid,varid,div_f,(/1,k/),(/ncellsf,1/))
         rc = nf90_inq_varid(ncid,'psi',varid)
         rc = nf90_put_var(ncid,varid,psi,(/1,k/),(/nverticesf,1/))
         rc = nf90_inq_varid(ncid,'chi',varid)
         rc = nf90_put_var(ncid,varid,chi,(/1,k/),(/ncellsf,1/))
      enddo
      deallocate(vorticityc,velocityc,divergencec)

      rc =  nf90_close(ncidd)
      
      call post_init_solver()

      rc = nf90_inq_varid(ncid,'lwghts_chi',varid)
      rc = nf90_put_var(ncid,varid,lwghts_chi)
      rc = nf90_inq_varid(ncid,'lwghts_psi',varid)
      rc = nf90_put_var(ncid,varid,lwghts_psi)
      rc = nf90_sync(ncid)

   else
      ! open the workfile and read some once-only variables
      rc = nf90_open(trim(workfile),nf90_write,ncid)
      rc = nf90_inq_varid(ncid,'lwghts_chi',varid)
      rc = nf90_get_var(ncid,varid,lwghts_chi)
      rc = nf90_inq_varid(ncid,'lwghts_psi',varid)
      rc = nf90_get_var(ncid,varid,lwghts_psi)
      rc = nf90_inq_varid(ncid,'resid_chi',varid)
      rc = nf90_get_var(ncid,varid,resid_chi)
      rc = nf90_inq_varid(ncid,'resid_psi',varid)
      rc = nf90_get_var(ncid,varid,resid_psi)
      rc = nf90_inq_varid(ncid,'iter_glob',varid)
      rc = nf90_get_var(ncid,varid,iter_glob)
      rc = nf90_inq_varid(ncid,'numiter',varid)
      rc = nf90_get_var(ncid,varid,numiter)
      rc = nf90_inq_varid(ncid,'maxiter',varid)
      rc = nf90_get_var(ncid,varid,maxiter)
      rc = nf90_inq_varid(ncid,'rlim',varid)
      rc = nf90_get_var(ncid,varid,rlim)
   endif

   print*,'done initialization/restart'   

!!!! end init/restart block

!!!! begin solver block
   
   converged = .false.
   do while (.not. converged)
      iter_glob = iter_glob + numiter
      converged = .true.
      ! 100 levels for normalVelocity, level 101 for normalBarotropicVelocity
      do k = 1,maxlevel
         ! iterative solver for chi
         if(resid_chi(k) > rlim .and. iter_glob <= maxiter) then
      
            indexn = 1
            indexf = 3-indexn
            chi_n => chi(:,indexn)
            chi_f => chi(:,indexf)
         
            rc = nf90_inq_varid(ncid,'div',varid)
            rc = nf90_get_var(ncid,varid,div_f,(/1,k/),(/ncellsf,1/))
            rc = nf90_inq_varid(ncid,'chi',varid)
            rc = nf90_get_var(ncid,varid,chi_n,(/1,k/),(/ncellsf,1/))
         
            call iterative_chi_solver(k)
            if(resid_chi(k) > rlim .and. iter_glob < maxiter) converged = .false.
         
            rc = nf90_put_var(ncid,varid,chi_n,(/1,k/),(/ncellsf,1/))
         
         endif

         ! iterative solver for psi
         if(resid_psi(k) > rlim .and. iter_glob <= maxiter) then
      
            indexn = 1
            indexf = 3-indexn
            psi_n => psi(:,indexn)
            psi_f => psi(:,indexf)
         
            rc = nf90_inq_varid(ncid,'vort',varid)
            rc = nf90_get_var(ncid,varid,vort_f,(/1,k/),(/nverticesf,1/))
            rc = nf90_inq_varid(ncid,'psi',varid)
            rc = nf90_get_var(ncid,varid,psi_n,(/1,k/),(/nverticesf,1/))
         
            call iterative_psi_solver(k)
            if(resid_psi(k) > rlim .and. iter_glob < maxiter)converged = .false.

            rc = nf90_put_var(ncid,varid,psi_n,(/1,k/),(/nverticesf,1/))
         
         endif
      enddo
         
      rc = nf90_inq_varid(ncid,'resid_chi',varid)
      rc = nf90_put_var(ncid,varid,resid_chi)
      rc = nf90_inq_varid(ncid,'resid_psi',varid)
      rc = nf90_put_var(ncid,varid,resid_psi)
      rc = nf90_inq_varid(ncid,'iter_glob',varid)
      rc = nf90_put_var(ncid,varid,iter_glob)

      rc = nf90_sync(ncid)

      print*,'iter, max resid_chi, max_resid_psi ', iter_glob,maxval(resid_chi(1:maxlevel)),maxval(resid_psi(1:maxlevel))  
   enddo
   
   print*,'done solution'
!!!! end solver block
stop   
!!!! begin final velocity block

   print*,'converged?',converged
   ! test for convergence, compute and write velocity
   if(converged) call compute_velocity_and_write()
   
   print*,'done velocity'
   
   ! close up the files we have written to
   rc = nf90_close(ncid)
   rc = nf90mpi_close(ncidf)
  
    call MPI_finalize(err)

   contains
      
!------------------------------------------

   subroutine iterative_chi_solver(k)
   integer, intent(in) :: k
   
   ! local variables
   integer iter, i, j
   real*8 del2, residtem
   
   !print*,'start chi solve, level ',k
   
   chi_f = chi_n
   
   !iterate this loop
   do iter = 1,numiter
      residtem = 0.0
      
!$omp  parallel do private(j,del2) reduction(max:residtem)
      do i = 1,ncellsf
         del2 = 0.
         do j = 1,nedgesOnCellf(i)
            if(cellsoncellf(j,i)>0) &
               del2 = del2 + lwghts_chi(j,i) * &
                  (chi_f(cellsoncellf(j,i)))
         enddo
         chi_n(i) = (del2-div_f(i))  
         residtem= max( residtem, abs( (del2-chi_f(i))-div_f(i) ) )
      enddo

!$omp end parallel do

      indexn = indexf
      indexf = 3-indexn
      chi_n => chi(:,indexn)
      chi_f => chi(:,indexf)
  enddo
  
  resid_chi(k) = residtem

   end subroutine iterative_chi_solver
   
!------------------------------------------

   subroutine iterative_psi_solver(k)
   integer, intent(in) :: k
   
   ! local variables
   integer iter, i, j
   real*8 del2, residtem
   
   !print*,'start psi solve, level ',k

   psi_f = psi_n
   
   !iterate this loop
   do iter = 1,numiter
      residtem = 0.0
      
!$omp parallel do private(j,del2) reduction(max:residtem)
      do i = 1,nverticesf
         del2 =  0.0
         do j = 1,3
            if(edgesOnVertexf(j,i)>0) then
               if(verticesOnEdgef(1,edgesOnVertexf(j,i)) == i) then
                  del2 = del2 +   &
                     lwghts_psi(j,i) * psi_f(verticesOnEdgef(2,edgesOnVertexf(j,i)))
               else
                  del2 = del2 +   &
                     lwghts_psi(j,i) * psi_f(verticesOnEdgef(1,edgesOnVertexf(j,i)))
               endif
            endif
         enddo
         residtem = max(residtem, abs( (del2-psi_f(i))-vort_f(i) ))
         psi_n(i) = del2-vort_f(i)
                          
      enddo
!$omp end parallel do
      indexn = indexf
      indexf = 3-indexn
      psi_n => psi(:,indexn)
      psi_f => psi(:,indexf)
   enddo

   resid_psi(k) = residtem

   end subroutine iterative_psi_solver
   
!------------------------------------------

   subroutine compute_velocity_and_write()

   ! compute fine grid velocities
   
   ! local variables
   real*8, dimension(:), allocatable :: dcedge, dvedge, velocity
   real*8, dimension(:,:), allocatable :: edgesignoncell, edgesignonvertex
   integer, dimension(:), allocatable :: maxlevelcell
   integer, dimension(:,:), allocatable :: cellsonedge, edgesoncell
   
   integer :: cell1, cell2, j, ivertex, vertex1, vertex2, iedge, iedgec, iedgev
                    
   
   allocate(cellsonedge(2,nedgesf))
   rc = nf90mpi_inq_varid(ncidf,'cellsOnEdge',varid)
   rc = nf90mpi_get_var_all(ncidf,varid,cellsonedge)
   allocate(edgesoncell(6,ncellsf))
   rc = nf90mpi_inq_varid(ncidf,'edgesOnCell',varid)
   rc = nf90mpi_get_var_all(ncidf,varid,edgesoncell)
   allocate(dcedge(nedgesf))
   rc = nf90mpi_inq_varid(ncidf,'dcEdge',varid)
   rc = nf90mpi_get_var_all(ncidf,varid,dcedge)
   allocate(dvedge(nedgesf))
   rc = nf90mpi_inq_varid(ncidf,'dvEdge',varid)
   rc = nf90mpi_get_var_all(ncidf,varid,dvedge)
   allocate(edgesignoncell(6,ncellsf))
   allocate(edgesignonvertex(3,nverticesf))
   allocate(velocity(nedgesf))
   allocate(maxlevelcell(ncellsf))
   rc = nf90mpi_inq_varid(ncidf,'maxLevelCell',varid)
   rc = nf90mpi_get_var_all(ncidf,varid,maxlevelcell)
   
!!$omp  parallel num_threads(1)

!!$omp  do private
      do iedge = 1, nedgesf
         dcedge(iedge) = 1./dcedge(iedge)
         dvedge(iedge) = 1./dvedge(iedge)
      enddo
      
!!$omp  do private(i,iedge)
      do iCell = 1, ncellsf
         edgeSignOnCell(:, icell)   = 0.0
         do i = 1, nedgesOnCellf(iCell)
            iEdge = edgesOnCell(i, iCell)

            if(iedge>0) then
            ! Vector points from cell 1 to cell 2
            if (iCell == cellsOnEdge(1, iEdge)) then
               edgeSignOnCell(i, iCell) = -1.0
            else
               edgeSignOnCell(i, iCell) =  1.0
            end if
            end if
         end do
      end do

!!$omp  do private(i,iedge)
      do iVertex = 1, nverticesf
         edgeSignOnVertex(:,ivertex) = 0.0
         do i = 1, 3
            iEdge = edgesOnVertexf(i, iVertex)

            if(iedge>0) then
            ! Vector points from vertex 1 to vertex 2
            if (iVertex == verticesOnEdgef(1,iEdge)) then
               edgeSignOnVertex(i,iVertex) = -1.0
            else
               edgeSignOnVertex(i,iVertex) =  1.0
            end if
            end if
         end do
      end do
!!$omp end parallel

   do k = 1,101
   
   print*,'write velocity k',k
   
   rc = nf90_inq_varid(ncid,'psi',varid)
   rc = nf90_get_var(ncid,varid,psi_f,(/1,k/),(/nedgesf,1/))
   rc = nf90_inq_varid(ncid,'chi',varid)
   rc = nf90_get_var(ncid,varid,chi_f,(/1,k/),(/ncellsf,1/))
   

!!$omp  parallel num_threads(1)
!!$omp  do private(dctem,dvtem,cell1,cell2,vertex1,vertex2,j,k,edgesign_temc,edgesign_temv)
      do i = 1,nedgesf
         velocity(i) = 0.0
         cell1 = cellsonedge(1,i)
         cell2 = cellsonedge(2,i)
         vertex1 = verticesonedgef(1,i)
         vertex2 = verticesonedgef(2,i)
         !TODO check signs
         if(minval(cellsonedge(:,i)) > 0) then
            do j = 1,nedgesOnCellf(cell1)
               if(i==edgesoncell(j,cell1))iedgec=j
            enddo
            do j = 1,3
               if(i==edgesOnVertexf(j,vertex1))iedgev=j
            enddo
            !divergent component
               velocity(i) = velocity(i) + &
                  (chi_f(cell1) - chi_f(cell2)) * dcedge(i) * edgeSignOnCell(iedgec,cell1)
         
            !irrotational component
               velocity(i) = velocity(i) - &
                  (psi_f(vertex1) - psi_f(vertex2)) * dvedge(i) * edgeSignOnVertex(iedgev,vertex1)
         
         endif
      
         !zero velocities at edges of basins
         if(cellsonedge(2,i)==0) then
            velocity(i) = 0.0
         else
            if(k>maxlevelcell(cellsonedge(1,i)) .and. k<=maxlevelcell(cellsonedge(2,i)) .or. &
               k>maxlevelcell(cellsonedge(2,i)) .and. k<=maxlevelcell(cellsonedge(1,i)) ) velocity(i) = 0.0
         endif
      
      enddo

!!$omp end parallel

         rms_v = 0.0
         do i = 1,nedgesf
            rms_v = rms_v + velocity(i)**2
         enddo
         print*,'rms_vel fine ',k,rms_v/nedgesf


      if(k < 101) then
         rc = nf90mpi_inq_varid(ncidf,'normalVelocity',varid)
         rc = nf90mpi_put_var_all(ncidf,varid,velocity,        &
                  (/k,1,1/)*1_mpi_offset_kind,(/1,nedgesf,1/)*1_mpi_offset_kind)
      else
         rc = nf90mpi_inq_varid(ncidf,'normalBarotropicVelocity',varid)
         rc = nf90mpi_put_var_all(ncidf,varid,velocity)
      endif

   enddo

   deallocate(maxlevelcell)
   deallocate(cellsonedge)
   deallocate(edgesoncell)
   deallocate(dcedge)
   deallocate(dvedge)
   deallocate(edgeSignOnCell, edgeSignOnVertex)
   
   end subroutine compute_velocity_and_write
   
!------------------------------------------

   subroutine pre_init_solver( )!{{{

      integer nctem, ivertex, icell, iedge, i, j, cell1, cell2, &
              iedgec, iedgev, vertex1, vertex2
      
      ! coarse grid mesh variables independent of level
      allocate(nedgesOnCellc(ncellsc))
      allocate(maxlevelcellc(ncellsc))
      allocate(edgesoncellc(6,ncellsc))
      allocate(edgesonvertexc(3,nverticesc))
      allocate(areacellc(ncellsc))
      allocate(invareacellc(ncellsc))
      allocate(areatrianglec(nverticesc))
      allocate(invareatrianglec(nverticesc))
      allocate(dvedgec(nedgesc))
      allocate(dcedgec(nedgesc))
      allocate(edgeSignOnCellc(6,ncellsc))
      allocate(edgeSignOnvertexc(3,nverticesc))
      allocate(cellsonedgec(2,nedgesc))
      allocate(minlevelvertextopc(nverticesc))
      allocate(maxlevelvertexbotc(nverticesc))
      allocate(minlevelcellc(ncellsc))
      allocate(cellsonvertexc(3,nverticesc))
      allocate(verticesoncellc(6,ncellsc))
      allocate(verticesonedgec(2,nedgesc))
      rc = nf90mpi_inq_varid(ncidc,'nEdgesOnCell',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,nedgesOnCellc)
      rc = nf90mpi_inq_varid(ncidc,'maxLevelCell',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,maxlevelcellc)
      rc = nf90mpi_inq_varid(ncidc,'edgesOnCell',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,edgesoncellc)
      rc = nf90mpi_inq_varid(ncidc,'dvEdge',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,dvedgec)
      rc = nf90mpi_inq_varid(ncidc,'dcEdge',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,dcedgec)
      rc = nf90mpi_inq_varid(ncidc,'cellsOnEdge',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,cellsonedgec)
      rc = nf90mpi_inq_varid(ncidc,'areaCell',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,areacellc)
      rc = nf90mpi_inq_varid(ncidc,'areaTriangle',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,areatrianglec)
      rc = nf90mpi_inq_varid(ncidc,'cellsOnEdge',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,cellsOnEdgec)
      rc = nf90mpi_inq_varid(ncidc,'cellsOnVertex',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,cellsOnVertexc)
      rc = nf90mpi_inq_varid(ncidc,'verticesOnCell',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,verticesOnCellc)
      rc = nf90mpi_inq_varid(ncidc,'verticesOnEdge',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,verticesOnedgec)
      rc = nf90mpi_inq_varid(ncidc,'edgesOnVertex',varid)
      rc = nf90mpi_get_var_all  (ncidc,varid,edgesonvertexc)
      invareacellc = 1.0/areacellc
      invareaTrianglec = 1.0 / areaTrianglec
 
      ! fine grid mesh variables
      allocate(maxlevelcellf(ncellsf))
      allocate(areacellf(ncellsf))
      allocate(areatrianglef(nverticesf))
      allocate(kiteareasf(3,nverticesf))
      allocate(cellsonvertexf(3,nverticesf))
      allocate(edgesoncellf(6,ncellsf))
      allocate(ncount(nverticesf))
      allocate(assigned_cell(ncellsf))
      allocate(flag(ncellsf))
      allocate(dvedgef(nedgesf))
      allocate(dcedgef(nedgesf))
      rc = nf90mpi_inq_varid(ncidf,'areaCell',varid)
      rc = nf90mpi_get_var_all  (ncidf,varid,areacellf)
      rc = nf90mpi_inq_varid(ncidf,'areaTriangle',varid)
      rc = nf90mpi_get_var_all  (ncidf,varid,areatrianglef)
      rc = nf90mpi_inq_varid(ncidf,'kiteAreas',varid)
      rc = nf90mpi_get_var_all  (ncidf,varid,kiteAreasf)
      rc = nf90mpi_inq_varid(ncidf,'cellsOnVertex',varid)
      rc = nf90mpi_get_var_all  (ncidf,varid,cellsOnVertexf)
      rc = nf90mpi_inq_varid(ncidf,'edgesOnCell',varid)
      rc = nf90mpi_get_var_all  (ncidf,varid,edgesoncellf)
      rc = nf90mpi_inq_varid(ncidf,'maxLevelCell',varid)
      rc = nf90mpi_get_var_all  (ncidf,varid,maxLevelCellf)
      rc = nf90mpi_inq_varid(ncidf,'dvEdge',varid)
      rc = nf90mpi_get_var_all  (ncidf,varid,dvedgef)
      rc = nf90mpi_inq_varid(ncidf,'dcEdge',varid)
      rc = nf90mpi_get_var_all  (ncidf,varid,dcedgef)

! compute
!!$omp  parallel 

      
!!$omp  do private(i, iedge)
      do iCell = 1, nCellsc
         do i = 1, nedgesOnCellc(iCell)
            iEdge = edgesOnCellc(i, iCell)
            edgeSignOnCellc (i,icell)  = 0.0

            ! Vector points from cell 1 to cell 2
               if (iCell == cellsOnEdgec(1, iEdge)) then
                  edgeSignOnCellc(i, iCell) = -1.0
               else
                  edgeSignOnCellc(i, iCell) =  1.0
               end if
         end do
      end do


!!$omp  do private(i, iedge)
      do iVertex = 1, nVerticesc
         do i = 1, 3
            edgeSignOnVertexc(i,ivertex) = 0.0
            iEdge = edgesOnVertexc(i, iVertex)

            if(iedge > 0) then
            ! Vector points from vertex 1 to vertex 2
            if (iVertex == verticesOnEdgec(1,iEdge)) then
               edgeSignOnVertexc(i,iVertex) = -1.0
            else
               edgeSignOnVertexc(i,iVertex) =  1.0
            end if
            endif
            
         end do
      end do


      ! minLevelVertexTop is the minimum (shallowest) of surrounding cells
      ! if the water column is dry, set minLevelCell outside of bounds
      ! to prevent dry cells from determining minimum
!!$omp  do 
      do iCell = 1, nCellsc
         if (maxLevelCellc(iCell) == 0) minLevelCellc(iCell) = 101 
      end do
!!$omp  do private(i)
      do iVertex = 1,nVerticesC
         minLevelVertexTopc(iVertex) = 100
         do i = 1, 3
            if(cellsOnVertexc(i,ivertex)>0) then
               minLevelVertexTopc(iVertex) = &
                  min( minLevelVertexTopc(iVertex), &
                       minLevelCellc(cellsOnVertexc(i,iVertex)))
            endif
         end do
      end do


      ! maxLevelVertexBot is the maximum (deepest) of surrounding cells
!!$omp  do private(i)
      do iVertex = 1,nVerticesc
         maxLevelVertexBotc(iVertex) = 0
         do i = 1, 3
            if(cellsOnVertexc(i,ivertex)>0) then
               maxLevelVertexBotc(iVertex) = &
                  max( maxLevelVertexBotc(iVertex), &
                       maxLevelCellc(cellsOnVertexc(i,iVertex)))
            endif
         end do
      end do

      
!!$omp  do private(j,iedge, lwghtsinv)
   do i = 1,ncellsf
      lwghts_chi(0:6,i) = 0.0
      do j = 1,nedgesOnCellf(i)
         if(cellsoncellf(j,i)>0) then
            iedge = edgesoncellf(j,i) 
            lwghts_chi(j,i) = dvedgef(iedge) / dcedgef(iedge)
            lwghts_chi(0,i) = lwghts_chi(0,i) + lwghts_chi(j,i)
         endif
      enddo
      lwghts_chi(0,i) = 1./lwghts_chi(0,i)
      lwghts_chi(1:6,i) = lwghts_chi(1:6,i) * lwghts_chi(0,i)
   enddo


!!$omp  do private(j,iedge)
   do i = 1,nverticesf
      lwghts_psi(0:3,i) = 0.0
      do j = 1,3
         if(edgesOnVertexf(j,i)>0) then
            iedge = edgesOnVertexf(j,i)
            lwghts_psi(j,i) = dcedgef(iedge) / dvedgef(iedge)
            lwghts_psi(0,i) = lwghts_psi(0,i) + lwghts_psi(j,i)
         endif
      enddo
      lwghts_psi(0,i) = 1./lwghts_psi(0,i)
      lwghts_psi(1:3,i) = lwghts_psi(1:3,i) * lwghts_psi(0,i)
   enddo
   
!!$omp end parallel 

      ALLOCATE(remap_indices(3,ncellsf))  !
      ALLOCATE(remap_wgts(3,ncellsf))     !
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
   rc = nf90_close(nctem)


   end subroutine pre_init_solver
   
!------------------------------------------

   subroutine post_init_solver( )!{{{

      ! deallocate mesh variables
      deallocate(nedgesOnCellc)
      deallocate(maxlevelcellc)
      deallocate(edgesoncellc)
      deallocate(areacellc)
      deallocate(dvedgec)
      deallocate(dcedgec)
      deallocate(edgeSignOnCellc)
      deallocate(cellsonedgec)
      deallocate(invareacellc)
      deallocate(areatrianglec)
      deallocate(invareatrianglec)
      deallocate(edgeSignOnvertexc)
      deallocate(minlevelvertextopc)
      deallocate(maxlevelvertexbotc)
      deallocate(minlevelcellc)
      deallocate(cellsonvertexc)
      deallocate(verticesoncellc)
      deallocate(verticesonedgec)
      deallocate(edgesOnVertexc)

      deallocate(maxlevelcellf)
      deallocate(areacellf)
      deallocate(areatrianglef)
      deallocate(kiteareasf)
      deallocate(cellsonvertexf)
      deallocate(edgesoncellf)
      deallocate(ncount)
      deallocate(assigned_cell)
      deallocate(flag)
      deallocate(dvedgef)
      deallocate(dcedgef)

      deallocate(remap_indices)  !
      deallocate(remap_wgts)     !
      deallocate(remap_indicesv)  !
      deallocate(remap_wgtsv)     !
      
   end subroutine post_init_solver
   
!------------------------------------------

   subroutine ocn_relativeVorticity_circulation( k )!{{{
      integer, intent(in) :: k

      ! This is computed on the coarse grid to be interpolated

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iVertex, iEdge, i
      real*8 ::  r_tmp, rms_v
      

      rms_v = 0.0
      !!$omp parallel do schedule(runtime) private(i, iEdge, r_tmp)
      do iVertex = 1, nVerticesc
         Vorticityc(iVertex) = 0.0
         do i = 1, 3
            iEdge = edgesOnVertexc(i, iVertex)
            if(iedge>0) then
               if ( k >= minlevelvertextopc(ivertex) .and. k <= maxLevelVertexBotc(iVertex) ) then
                  r_tmp = dcEdgec(iEdge) * Velocityc(iEdge)
                  Vorticityc(iVertex) = Vorticityc(iVertex) +              &
                       edgeSignOnVertexc(i, iVertex) * r_tmp * invareaTrianglec(iVertex)
               endif
            endif
         end do
         rms_v = rms_v + vorticityc(ivertex)**2
      end do
      !!$omp end parallel do
      print*,'rms_v coarse ',k,rms_v/nverticesc
      

!--------------------------------------------------------------------

   end subroutine ocn_relativeVorticity_circulation!}}}

!--------------------------------------------------------------------

   subroutine ocn_divergence( k )
      integer, intent(in) :: k

      integer :: icell, iedge, i
      real*8 :: r_tmp, rms_d
      
      rms_d = 0.0
      !!$omp parallel do schedule(runtime) &
      !!$omp private(i, iEdge, r_tmp )
      do iCell = 1, nCellsc
         divergencec( iCell) = 0.0
         do i = 1, nedgesOnCellc(iCell)
            iEdge = edgesOnCellc(i, iCell)
            if(iedge > 0) then
               if (k <= maxlevelcellc(icell) .or. k==101) then
                  r_tmp = dvEdgec(iEdge) * Velocityc(iEdge) * invareacellc(icell)
                  divergencec(iCell) = divergencec(iCell) -                 &
                           edgeSignOnCellc(i, iCell) * r_tmp
               end if
            endif
         end do
         rms_d = rms_d + divergencec(icell)**2
      end do
      !!$omp end parallel do
      print*,'rms_d coarse ',k,rms_d/ncellsc

   end subroutine ocn_divergence

!--------------------------------------------------------------------
   
   subroutine initialize_psi_solver(k)
   integer, intent(in) :: k
   
   integer i, l, m
      real*8 vort_bar, rms_v
      
   ! compute coarse grid vorticity 
   if(k < 101) then
      call ocn_relativeVorticity_circulation( k )
   else
      vorticityc(:) = 0.0
   endif

   
   !print*,'extrema vorticity_c ',k,minval(vorticityc),maxval(vorticityc)   
   
   ! interpolate vort from coarse grid vertices to fine grid vertices
   rms_v = 0.0
      !!$omp parallel do private(l,m,vort_bar)
      DO i = 1,nverticesf
         vort_f(i) = 0.0
         DO l = 1,3
            IF(remap_indicesv(l,i) > 0)     &
               vort_f(i) = vort_f(i) +     &
                             vorticityc(remap_indicesv(l,i)) * remap_wgtsv(l,i)
         ENDDO
         IF(remap_indicesv(4,i) > 0) then
            vort_bar = 0.0
            do m = 1,nedgesOnCellc(remap_indicesv(4,i))
               vort_bar = vort_bar + vorticityc(verticesoncellc(m,remap_indicesv(4,i)))
            enddo
            vort_bar = vort_bar / nedgesOnCellc(remap_indicesv(4,i))
            vort_f(i) = vort_f(i) +     &
                           vort_bar * remap_wgtsv(4,i)
         endif
         rms_v = rms_v + vort_f(i)**2
      ENDDO
      !!$omp end parallel do
      
   !print*,'extrema vorticity_f ',k,minval(vort_f),maxval(vort_f)
      print*,'rms_v fine ',k,rms_v/nverticesf

         rc = nf90_inq_varid(ncidd,'vort',varid)
         rc = nf90_put_var(ncidd,varid,vorticityc,(/1,k/),(/nverticesc,1/))
   
   end subroutine initialize_psi_solver

!------------------------------------------
   
   subroutine initialize_chi_solver(k)
   integer, intent(in) :: k
   
   integer i, j, m, cell1, cell2, vertex1, vertex2,   &
           icell, iedge, iedgec, iedgev, iter, ivertex
   real*8 rms_d
   
   call ocn_divergence(k)
   
   !print*,'extrema divergence_c ',minval(divergencec),maxval(divergencec)
   rms_d = 0.0
   
   ! interpolate div from coarse grid cells to fine grid cells

      !!$omp parallel do private(l)
      DO i = 1,ncellsf
         div_f(i) = 0.0
         DO l = 1,3
            IF(remap_indices(l,i) > 0)     &
               div_f(i) = div_f(i) +     &
                             divergencec(remap_indices(l,i)) * remap_wgts(l,i)
         ENDDO
         rms_d = rms_d + div_f(i)**2
      ENDDO
      !!$omp end parallel do
   !print*,'extrema divergence_f ',minval(div_f),maxval(div_f)
      print*,'rms_d fine ',k,rms_d/ncellsf
      
         rc = nf90_inq_varid(ncidd,'div',varid)
         rc = nf90_put_var(ncidd,varid,divergencec,(/1,k/),(/ncellsc,1/))

   end subroutine initialize_chi_solver

!------------------------------------------

   subroutine adjust_vort_div_means(k)
   
   integer, intent(in) :: k
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  get basin info for depth

! At each level we want
!    number of mask cells, number of real cells
!       list of real cell indices
!    number of basins
!       for each basin we want
!           number of cells in the basin
!           a list of member cell indices in the basin
!
! Sanity checks:
!   number of real cells equals sum across basins of number of basin members
!   each real cell occurs in exactly one and only one basin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer i, nbasins, istart, icount, iter, count1, nmasked, nunmasked, nn, n
   logical lcount1
   real*8 divmean, vortmean, areamean, temmean1, temmean2
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$omp parallel

         !assign zero to the bathymetry masked cells
!!$omp  do 
         do i = 1, ncellsf
            assigned_cell(i) = -1   
            if( k > maxlevelcellf(i)) assigned_cell(i) = 0
         enddo
         
         nbasins = 0
         do while ( sum(assigned_cell(:)) .ne. sum(abs(assigned_cell(:))) ) 
            nbasins = nbasins + 1
!print*,'nbasins',nbasins
            
            ! Start a new basin search with the first unassigned cell
            do i = 1,ncellsf
               if (assigned_cell(i) < 0) then
                  istart = i
                  exit
               endif
            enddo
!!$omp  do 
            do i = 1,ncellsf
               flag(i) = 0
               if (assigned_cell(i) .ge. 0) flag(i) = -1
            enddo
            flag(istart) = 1
!!print*,'flag initialized'
            iter = 0
            count1 = 1
            
            DO WHILE (count1>0)
               iter = iter + 1
!!print*,'iter,count1',iter,count1
               lcount1 = count1==1 .and. iter>1
!!$omp  do private(m) reduction(+:count1) 
               DO i = 1,ncellsf
                  IF(flag(i)==1) THEN
                     flag(i) = 2
!if(lcount1)print*,'flag1to2',i
                     assigned_cell(i) = nbasins
                     count1 = count1 - 1
                     DO m = 1,nedgesOnCellf(i)
                        IF(cellsoncellf(m,i)>0 ) then
                        if ( flag(cellsoncellf(m,i))==0)THEN
                        !!$OMP critical
                           flag(cellsoncellf(m,i)) = 1
!if(lcount1)print*,'flagto1',m,i
                           count1 = count1 + 1
                        !!$OMP end critical
                        ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
!            if(nbasins==1)then
!            ! assigned_cell counts
!               acount = 0
!               do i = 1,ncellsf
!                  if(assigned_cell(i)==-1)acount = acount + 1
!               enddo
!               !print*,'assigned_cell -1 count ',acount
!               acount = 0
!!$omp  do reduction(+:acount) 
!               do i = 1,ncellsf
!                  if(assigned_cell(i)==0)acount = acount + 1
!               enddo
!               !print*,'assigned_cell 0 count ',acount
!               acount = 0
!!$omp  do reduction(+:acount) 
!               do i = 1,ncellsf
!                  if(assigned_cell(i)==1)acount = acount + 1
!               enddo
!               !print*,'assigned_cell 1 count ',acount
!            ! flag counts
!               acount = 0
!!$omp  do reduction(+:acount) 
!               do i = 1,ncellsf
!                  if(flag(i)==-1)acount = acount + 1
!               enddo
!               !print*,'flag -1 count ',acount
!               acount = 0
!!$omp  do reduction(+:acount) 
!               do i = 1,ncellsf
!                  if(flag(i)==0)acount = acount + 1
!               enddo
!              !print*,'flag 0 count ',acount
!               acount = 0
!!$omp  do reduction(+:acount) 
!               do i = 1,ncellsf
!                  if(flag(i)==1)acount = acount + 1
!               enddo
!              !print*,'flag 1 count ',acount
!               acount = 0
!!$omp  do reduction(+:acount) 
!               do i = 1,ncellsf
!                  if(flag(i)==2)acount = acount + 1
!               enddo
!               !print*,'flag 2 count ',acount
!            endif
            if(iter > 5000) stop   
            ENDDO
            !PRINT*,'iter = ',iter         
            ! starting a index istart iterate until all cells in its basin are found
         enddo

!print*,'Level',k,'summary'  
   nmasked = 0
   nunmasked = 0
   do i = 1,ncellsf
      if(assigned_cell(i) == 0) nmasked = nmasked + 1
      if(assigned_cell(i) > 0) nunmasked = nunmasked + 1
   enddo
!print*,'   masked, real, total cells ',  nmasked, nunmasked, nmasked + nunmasked      
!print*,'   number of basins ', nbasins     
      do nn = 1,nbasins
   nunmasked = 0.0
   do i = 1,ncellsf
      if(assigned_cell(i) == nn) nunmasked = nunmasked + 1
   enddo
!print*,'  Basin ',nn,' count = ',nunmasked

            divmean = 0.0
            vortmean = 0.0
            areamean = 0.0
!!$omp  do reduction(+:areamean,divmean) 
            do n = 1,ncellsf
               if(k.le.maxlevelcellf(n)) then
               if(assigned_cell(n)==nn)then
                  areamean = areamean + areacellf(n)
                  divmean = divmean + areacellf(n) * div_f(n)
               endif
               endif
            enddo
            divmean = divmean / areamean
            temmean1 = 0.0
            temmean2 = 0.0
!!$omp  do reduction(+:temmean1,temmean2)
            do n = 1,ncellsf
               if(k.le.maxlevelcellf(n) .and. assigned_cell(n)==nn)then
                  div_f(n) = div_f(n) - divmean 
                  temmean1 = temmean1 + div_f(n) * areacellf(n)
                  temmean2 = temmean2 + areacellf(n)
               endif
            enddo
            
!            print*,'k,basin, divmean, maxabs',k,nn,temmean1/temmean2,maxval(abs(div_f(k,:)))
!!$omp  do private(m) reduction(+:vortmean) 
            do n = 1,nverticesf
               ncount(n) = 0
               do m = 1,3
                  if(cellsonvertexf(m,n)>0) then
                  if(assigned_cell(cellsonvertexf(m,n))==nn) then
                  if(k .le. maxlevelcellf(cellsonvertexf(m,n)) )   &
                     ncount(n) = ncount(n) + 1
                  endif
                  endif
               enddo
               if (ncount(n) > 1) then
                  do m = 1,3
                     if(cellsonvertexf(m,n)>0) then
                     if(assigned_cell(cellsonvertexf(m,n))==nn) then
                     if(k .le. maxlevelcellf(cellsonvertexf(m,n)))   &
                        vortmean = vortmean + vort_f(n) * kiteareasf(m,n)
                     endif
                     endif
                  enddo
               endif
            enddo
            vortmean = vortmean / areamean
!!$omp  do 
            do n = 1,nverticesf
               if (ncount(n) > 1) &
                  vort_f(n) = vort_f(n) - vortmean
            enddo
         enddo ! do nn, loop over basins
!!$omp end parallel

!!$omp  paralllel
!!$omp  do private
   do i = 1,ncellsf
      div_f(i) = div_f(i) * areacellf(i) *  lwghts_chi(0,i)
   enddo


!!$omp parallel do 
   do i = 1,nverticesf
      vort_f(i) = vort_f(i) * lwghts_psi(0,i)  * areatrianglef(i)
   enddo

!!$omp end parallel
!print*,'div vort mean adjusted'

   end subroutine adjust_vort_div_means

!------------------------------------------

end program build_interpolated_velocities
