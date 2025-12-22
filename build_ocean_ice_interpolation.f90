program build_ocean_ice_interpolation

! compile with netcdf and, optionally openmp

use netcdf

implicit none

!!!
! Declarations

   ! master ocean data files
   character(len=128) :: finefile, coarsefile, outfile
   !resolution specific variables
   integer caseres

   ! netcdf-related variables
   !   'c' suffix is coarse grid
   !   'f'  suffix is finer grid
   integer :: ncidf, ncidc, rc, dimid, varid
   integer :: dimid1, dimid2, dimid3, dimid4, dimid5

   ! grid variables
   
   real*8 pi, d2r
   
   ! remapping arrays
   real*4, dimension(:,:), allocatable :: remap_wgts, remap_wgts_a
   integer, dimension(:,:), allocatable :: remap_indices, remap_indices_a
   real*4, dimension(:,:), allocatable :: remap_wgtsv
   integer, dimension(:,:), allocatable :: remap_indicesv
   
   integer :: err, nctem
  
      INTEGER :: ncellsc, nverticesc, nedgesc
      INTEGER :: ncellsf, nverticesf, nedgesf
      

!!!
! Executable code begins
   pi = 4._8 * atan(1._8)
   d2r = pi/180.
!!!
! preliminaries
   print*,'enter the input ocean restart file name'
   read(5,'(a128)')coarsefile
   print*,'enter the ocean mesh file name'
   read(5,'(a128)')finefile
   print*,'enter the interpolation weights file name to be created'
   read(5,'(a128)')outfile

   ! open the necessary netcdf files
   rc = nf90_open(trim(coarsefile),nf90_nowrite,ncidc)
   if(rc .ne. 0) stop 'trouble opening the input ocean restart - stopping'

   rc = nf90_open(trim(finefile),nf90_nowrite,ncidf)
   if(rc .ne. 0) stop 'trouble opening the ocean mesh file - stopping'

      ! read the coarse grid dimensions
      rc = nf90_inq_dimid(ncidc,'nCells',dimid)
      rc = nf90_inquire_dimension(ncidc,dimid,len=ncellsc)
      rc = nf90_inq_dimid(ncidc,'nVertices',dimid)
      rc = nf90_inquire_dimension(ncidc,dimid, len=nverticesc)
      rc = nf90_inq_dimid(ncidc,'nEdges',dimid)
      rc = nf90_inquire_dimension(ncidc,dimid, len=nedgesc)

print*,'coarse grid',ncellsc,nverticesc,nedgesc
      ! read the fine grid dimensions
      rc = nf90_inq_dimid(ncidf,'nCells',dimid)
      rc = nf90_inquire_dimension(ncidf,dimid,len=ncellsf)
      rc = nf90_inq_dimid(ncidf,'nVertices',dimid)
      rc = nf90_inquire_dimension(ncidf,dimid, len=nverticesf)
      rc = nf90_inq_dimid(ncidf,'nEdges',dimid)
      rc = nf90_inquire_dimension(ncidf,dimid, len=nedgesf)
print*,'fine grid',ncellsf,nverticesf,nedgesf

print*,'start initialize_interpolate_variable'
      ALLOCATE(remap_indices(3,ncellsf))  !
      ALLOCATE(remap_wgts(3,ncellsf))     !
      ALLOCATE(remap_indices_a(3,ncellsf))  !
      ALLOCATE(remap_wgts_a(3,ncellsf))     !
      ALLOCATE(remap_indicesv(4,nverticesf))  !
      ALLOCATE(remap_wgtsv(4,nverticesf))     !

      rc = nf90_create(trim(outfile),ior(nf90_clobber,nf90_netcdf4),nctem)
   if(rc .ne. 0) stop 'trouble creating the interpolation weights file - stopping'
      rc = nf90_def_dim(nctem,'nCells',ncellsf,dimid1)
      print*,'defdim1 ',rc
      rc = nf90_def_dim(nctem,'nVertices',nverticesf,dimid2)
      print*,'defdim2 ',rc
      rc = nf90_def_dim(nctem,'THREE',3,dimid3)
      print*,'defdim3 ',rc
      rc = nf90_def_dim(nctem,'FOUR',4,dimid4)
      print*,'defdim4 ',rc
      rc = nf90_def_dim(nctem,'nLevels',100,dimid5)
      print*,'defdim5 ',rc
      rc = nf90_def_var(nctem,'remap_indices',NF90_INT,(/dimid3,dimid1,dimid5/),varid)
      print*,'defvar1 ',rc
      rc = nf90_def_var(nctem,'remap_wgts',NF90_float,(/dimid3,dimid1,dimid5/),varid)
      print*,'defvar2 ',rc
      rc = nf90_def_var(nctem,'remap_indicesv',NF90_INT,(/dimid4,dimid2/),varid)
      print*,'defvar3 ',rc
      rc = nf90_def_var(nctem,'remap_wgtsv',NF90_float,(/dimid4,dimid2/),varid)
      print*,'defvar4 ',rc
      rc = nf90_enddef(nctem)
      print*,'enddef ',rc

      CALL initialize_interpolate_variable( ncellsf, nverticesf )
      
      rc = nf90_inq_varid(nctem,'remap_indicesv',varid)
      print*,'inqvar3 ',rc
      rc = nf90_put_var(nctem,varid,remap_indicesv)
      print*,'putvar3 ',rc
      rc = nf90_inq_varid(nctem,'remap_wgtsv',varid)
      print*,'inqvar4 ',rc
      rc = nf90_put_var(nctem,varid,remap_wgtsv)
      print*,'putvar4 ',rc
      rc = nf90_close(nctem)
      print*,'close ',rc
print*,'done initialize_interpolate_variable'     
      
   contains
      
!------------------------------------------

      SUBROUTINE initialize_interpolate_variable( ncellsf, nverticesf )
      ! This routine will initialize a horizontal interpolation.
      ! Method:
      ! 1) For each fine grid cell location identify the nearest
      !    coarse grid vertex.
      ! 2) For the identified fine grid cell / coarse grid vertex pair
      !    generate interpolation weights for the three coarse grid cells
      !    that neighbor the vertex.
      ! 3) For ocean depths, this method is modified to only search for coarse
      !    grid vertices that are above the ocean floor.
      
      ! Argument list
      INTEGER, INTENT(IN) :: ncellsf, nverticesf
      
      ! local variables
      ! fine grid
      REAL*8, DIMENSION(ncellsf) :: xcell_f, ycell_f, zcell_f ! cell locations in cartesian coordinates
      REAL*8, DIMENSION(nverticesf) :: xvert_f, yvert_f, zvert_f ! vertex locations in cartesian coordinates
      ! coarse grid
      REAL*8, DIMENSION(:), ALLOCATABLE :: xcell_c, ycell_c, zcell_c ! cell locations in cartesian coordinates
      REAL*8, DIMENSION(:), ALLOCATABLE :: xvert_c, yvert_c, zvert_c ! vertex locations in cartesian coordinates
      INTEGER, DIMENSION(:), ALLOCATABLE :: maxlevelcellc, nedgesoncell, maxlevelcellf
      INTEGER, DIMENSION(:), ALLOCATABLE :: maxlevelvertexc
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: cellsonvertexc, verticesonedge, verticesoncell
      REAL*8 :: radius, dotmax, dotloc, d1, d2, d3,  &
                x1,y1,x2,y2,x3,y3, xp, yp, denom, wsum, wsum1
      
      INTEGER :: rc, varid, dimid, wsumloc
      INTEGER :: n, m, k, mmax, zeroloc(3), l, l1, l2, l3, lmax, zeroloc0, zeroloc1, zeroloc2, zeroloc3
      logical l_useabove, l_match
      REAL*8, DIMENSION(:), ALLOCATABLE :: wsummax ! cell locations in cartesian coordinates

      print*,'enter initialize_interpolate_variable'
      ! read the fine grid variables
      rc = nf90_inq_varid(ncidf,'xCell',varid)
      rc = nf90_get_var(ncidf,varid,xcell_f)
      rc = nf90_inq_varid(ncidf,'yCell',varid)
      rc = nf90_get_var(ncidf,varid,ycell_f)
      rc = nf90_inq_varid(ncidf,'zCell',varid)
      rc = nf90_get_var(ncidf,varid,zcell_f)
      rc = nf90_inq_varid(ncidf,'xVertex',varid)
      rc = nf90_get_var(ncidf,varid,xvert_f)
      rc = nf90_inq_varid(ncidf,'yVertex',varid)
      rc = nf90_get_var(ncidf,varid,yvert_f)
      rc = nf90_inq_varid(ncidf,'zVertex',varid)
      rc = nf90_get_var(ncidf,varid,zvert_f)
      
      ! read the coarse grid variables
      rc = nf90_get_att(ncidf, nf90_global, 'sphere_radius', radius)
      print*,'fine radius',radius
      radius = 1./radius
      
      !map to unit sphere
      xcell_f = xcell_f * radius
      ycell_f = ycell_f * radius
      zcell_f = zcell_f * radius
      xvert_f = xvert_f * radius
      yvert_f = yvert_f * radius
      zvert_f = zvert_f * radius
      
      ALLOCATE(maxlevelcellc(ncellsc))
      ALLOCATE(maxlevelcellf(ncellsf))
      ALLOCATE(wsummax(ncellsf))
      ALLOCATE(maxlevelvertexc(nverticesc))
      ALLOCATE(cellsonvertexc(3,nverticesc))
      ALLOCATE(verticesonedge(2,nverticesc))
      ALLOCATE(verticesoncell(6,ncellsc))
      ALLOCATE(nedgesoncell(ncellsc))
      ALLOCATE(xcell_c(ncellsc))
      ALLOCATE(ycell_c(ncellsc))
      ALLOCATE(zcell_c(ncellsc))
      ALLOCATE(xvert_c(nverticesc))
      ALLOCATE(yvert_c(nverticesc))
      ALLOCATE(zvert_c(nverticesc))

      rc = nf90_inq_varid(ncidc,'xCell',varid)
      rc = nf90_get_var(ncidc,varid,xcell_c)
      rc = nf90_inq_varid(ncidc,'yCell',varid)
      rc = nf90_get_var(ncidc,varid,ycell_c)
      rc = nf90_inq_varid(ncidc,'zCell',varid)
      rc = nf90_get_var(ncidc,varid,zcell_c)
      rc = nf90_inq_varid(ncidc,'xVertex',varid)
      rc = nf90_get_var(ncidc,varid,xvert_c)
      rc = nf90_inq_varid(ncidc,'yVertex',varid)
      rc = nf90_get_var(ncidc,varid,yvert_c)
      rc = nf90_inq_varid(ncidc,'zVertex',varid)
      rc = nf90_get_var(ncidc,varid,zvert_c)
      rc = nf90_inq_varid(ncidc,'maxLevelCell',varid)
      rc = nf90_get_var(ncidc,varid,maxlevelcellc)
      rc = nf90_inq_varid(ncidf,'maxLevelCell',varid)
      rc = nf90_get_var(ncidf,varid,maxlevelcellf)
      rc = nf90_inq_varid(ncidc,'cellsOnVertex',varid)
      rc = nf90_get_var(ncidc,varid,cellsonvertexc)
      rc = nf90_inq_varid(ncidc,'verticesOnCell',varid)
      rc = nf90_get_var(ncidc,varid,verticesOnCell)
      rc = nf90_inq_varid(ncidc,'nEdgesOnCell',varid)
      rc = nf90_get_var(ncidc,varid,nEdgesOnCell)
      rc = nf90_get_att(ncidc, nf90_global, 'sphere_radius', radius)
      print*,'coarse radius',radius
      radius = 1./radius

      !map to unit sphere
      xcell_c = xcell_c * radius
      ycell_c = ycell_c * radius
      zcell_c = zcell_c * radius
      xvert_c = xvert_c * radius
      yvert_c = yvert_c * radius
      zvert_c = zvert_c * radius

print*,'coarse grid',ncellsc,nverticesc,nedgesc
print*,'fine grid',ncellsf,nverticesf,nedgesf
      
      remap_indicesv = 0
      remap_wgtsv = 0.0
      ! for each vertex, find the max of the maxlevelcell of the neighbor cellsoncell
      DO m = 1,nverticesc
         maxlevelvertexc(m) = 0
         DO l = 1,3
            if(cellsonvertexc(l,m)>0)  &
            maxlevelvertexc(m) = max(maxlevelvertexc(m), maxlevelcellc(cellsonvertexc(l,m)))
         ENDDO
      ENDDO
      
!print*,'maxlevelcell',maxlevelcell(1563)
      DO k = 1,100
         remap_indices = 0
         remap_wgts = 0.0
         wsummax = 0.0
         IF(k>1 .and. k .le. minval(maxlevelcellf)) then
            remap_indices(:,:) = remap_indices_a(:,:)
            remap_wgts(:,:) = remap_wgts_a(:,:)
         ELSE
print*,'layer ', k
         ! Find the coarse grid vertex nearest the fine grid cell
         zeroloc0=0
          zeroloc1=0
           zeroloc2=0
            zeroloc3=0
!$omp parallel do private(dotmax,mmax,m,dotloc,zeroloc,l,l1,l2,l3,wsum,d1,d2,d3,l_useabove,l_match,denom,x1,y1,x2,y2,x3,y3,xp,yp) reduction(+ : zeroloc0,zeroloc1,zeroloc2,zeroloc3)
         DO n = 1,ncellsf
         
           l_useabove = .false.
           if(k > 1) then
              l_match = .true.
              do l = 1,3
                 if(remap_indices_a(l,n)>0) then
                    if(k>maxlevelcellc(remap_indices_a(l,n))) then 
                       l_match = .false.
                       exit
                    endif
                 endif
              enddo
              if(l_match) l_useabove = .true.
           endif
         
           if(l_useabove ) then
            
            remap_indices(:,n) = remap_indices_a(:,n)
            remap_wgts(:,n) = remap_wgts_a(:,n)
 
           else 
            
            dotmax = -1
            mmax = 0
            IF(k .le. maxlevelcellf(n)) THEN
               DO m = 1,nverticesc
                  IF(k .LE. maxlevelvertexc(m)) THEN
                     dotloc = xcell_f(n)*xvert_c(m) + &
                              ycell_f(n)*yvert_c(m) + &
                              zcell_f(n)*zvert_c(m)
                     IF(dotloc > dotmax) THEN
                        dotmax = dotloc
                        mmax = m
                     ENDIF
                  ENDIF
               ENDDO
!               if(n==224127) then
!                  print*,'dotmax,mmax ',dotmax,mmax
!                  print*,'cellsonvertexc',cellsonvertexc(:,mmax)
!                  print*,'xyzc ',xvert_c(mmax),yvert_c(mmax),zvert_c(mmax)
!                  print*,'xyzf ',xcell_f(n),ycell_f(n),zcell_f(n)
!               endif
               zeroloc = 0
               DO l = 1,3
                  IF(cellsonvertexc(l,mmax)>0)THEN
                     IF(k .le. maxlevelcellc(cellsonvertexc(l,mmax)))THEN
                        remap_indices(l,n) = cellsonvertexc(l,mmax)
                     ELSE
                        zeroloc(l) = 1
                     ENDIF
                  ELSE
                     zeroloc(l) = 1
                  ENDIF
               ENDDO
               !if(n==7310)print*,'7310,zeroloc',zeroloc
               !if(n==7310)print*,'7310,remap_indices(l,n,k) ',remap_indices(:,n,k) 
				!print*,'mmax',mmax
				!print*,'dotmax',dotmax
				!print*,'maxlevelvertexc(m)',maxlevelvertexc(mmax)
				!print*,'cellsonvertexc(m)',cellsonvertexc(:,mmax)
    			!print*,'zeroloc', zeroloc
               SELECT CASE( SUM(zeroloc))
                  CASE(0) ! all cell neighbors of vertex contribute to interpolation
                     ! compute a triangle in 2D for barycentric coordinates
                     zeroloc0 = zeroloc0 + 1
                     l1 = remap_indices(1,n)
                     l2 = remap_indices(2,n)
                     l3 = remap_indices(3,n)
                     CALL compute_triangle(            &
                             (/xcell_c(l1),ycell_c(l1),zcell_c(l1)/) ,  &
                             (/xcell_c(l2),ycell_c(l2),zcell_c(l2)/) ,  &
                             (/xcell_c(l3),ycell_c(l3),zcell_c(l3)/) ,  &
                             (/xcell_f(n),ycell_f(n),zcell_f(n)/) ,     &
                             x1, y1, x2, y2, x3, y3, xp, yp )
                     ! now that I have a triangle in a 2D plane I can use the 
                     !   formulae in https://codeplea.com/triangular-interpolation
                     !   for the weights
                     denom = 1./ ( (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3) )
                     
                     remap_wgts(1,n) = denom * ( (y2-y3)*(xp-x3) + (x3-x2)*(yp-y3) )
                     remap_wgts(2,n) = denom * ( (y3-y1)*(xp-x3) + (x1-x3)*(yp-y3) )
                     remap_wgts(3,n) = 1. - (remap_wgts(1,n)+remap_wgts(2,n))
                     
                !if(n==7310)print*,'7310,remap_wgts1 ',remap_wgts(:,n,k) 
     				!print*,'p1',x1,y1
      				!print*,'p2',x2,y2
      				!print*,'p3',x3,y3
      				!print*,'pp',xp,yp
      				!print*,'denom',denom
      				!print*,'remap_wgts',remap_wgts(:,1,1)
                     wsum = sum(abs(remap_wgts(:,n)))
                     
                     wsummax(n) = wsum
                     !!$OMP critical
                     !if(wsum > wsummax) then
                     !   wsummax = wsum
                     !   wsumloc = n
                     !endif
                     !!$OMP end critical
                     
                     !if( n == 224127) then
                     !   print*,'remap_wgts 224127 ',remap_wgts(:,n,k) 
                     !   print*,'remap_indices 224127 ',remap_indices(:,n,k) 
                     !endif
                     
                     ! correct if sum of abs(weights) > 1 (this means one weight is negative)
                     IF(wsum > 1.0000001) THEN
                     !print*,'case0',n, wsum
                        DO m = 1,3
                           IF(remap_wgts(m,n) < 0.) THEN
                              remap_wgts(m,n) = 0.0
                              wsum = sum(remap_wgts(:,n))
                              remap_wgts(:,n) = remap_wgts(:,n) / wsum
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF
                !if(n==7310)print*,'7310,remap_wgts2 ',remap_wgts(:,n,k) 
                  CASE(1) ! two cell neighbors of vertex contribute to interpolation
                     zeroloc1 = zeroloc1 + 1
                     !print*,'case1',n
                     IF( zeroloc(1) == 1) THEN
                        l1 = 2
                        l2 = 3
                     ENDIF
                     IF( zeroloc(2) == 1) THEN
                        l1 = 1
                        l2 = 3
                     ENDIF
                     IF( zeroloc(3) == 1) THEN
                        l1 = 1
                        l2 = 2
                     ENDIF
                     m = remap_indices(l1,n)
                     d1 = xcell_f(n)*xcell_c(m) + &
                          ycell_f(n)*ycell_c(m) + &
                          zcell_f(n)*zcell_c(m)
                     d1 = acos(max(-1.0,min(1.0,d1)))
                     m = remap_indices(l2,n)
                     d2 = xcell_f(n)*xcell_c(m) + &
                          ycell_f(n)*ycell_c(m) + &
                          zcell_f(n)*zcell_c(m)
                     d2 = acos(max(-1.0,min(1.0,d2)))
                     d3 = 1. / (d1+d2)
                     remap_wgts(l1,n) = d2 * d3
                     remap_wgts(l2,n) = d1 * d3
                  CASE(2) ! one cell neighbors of vertex contribute to interpolation
                     !print*,'case2',n
                     zeroloc2 = zeroloc2 + 1
                     DO l = 1,3
                        IF( zeroloc(l) == 0) remap_wgts(l,n) = 1.0
                     ENDDO
               END SELECT
            ENDIF
            
          endif
          
         ENDDO  ! n loop
!$OMP End parallel do
         print*,'wsummax',maxval(wsummax),maxloc(wsummax)
         print*,'zeroloc',zeroloc0,zeroloc1,zeroloc2,zeroloc3
      ENDIF
      
      ! one more check for negative weights
      if(minval(remap_wgts(:,:)) < 0.0)then
         print*,'k, minwgt ',k,minval(remap_wgts(:,:))
!$omp parallel do private(wsum,l)
         do n = 1,ncellsf
            if(minval(remap_wgts(:,n))<0.0) then
               wsum = 0.0
               do l = 1,3
                  remap_wgts(l,n) = max(0.0,remap_wgts(l,n))
                  wsum = wsum + remap_wgts(l,n)
               enddo
               remap_wgts(:,n) = remap_wgts(:,n) / wsum
            endif
                !if(n==7310)print*,'7310,remap_wgts3 ',remap_wgts(:,n,k) 
         enddo
!$OMP End parallel do
      endif
      
      rc = nf90_inq_varid(nctem,'remap_indices',varid)
      print*,'inqvar1 ',rc
      rc = nf90_put_var(nctem,varid,remap_indices,(/1,1,k/),(/3,ncellsf,1/))
      print*,'putvar1 ',rc
      rc = nf90_inq_varid(nctem,'remap_wgts',varid)
      print*,'inqvar2 ',rc
      rc = nf90_put_var(nctem,varid,remap_wgts,(/1,1,k/),(/3,ncellsf,1/))
      print*,'putvar2 ',rc

      rc = nf90_sync(nctem)
      
      remap_indices_a = remap_indices
      remap_wgts_a = remap_wgts
   enddo ! k loop

!      print*,'remap_indices ',minval(remap_indices(:,:,k)),maxval(remap_indices(:,:,k))
!      print*,'remap_wgts ',minval(remap_wgts(:,:,k)),maxval(remap_wgts(:,:,k))
!      do n = 1,ncellsf
!         write(83,*)n,remap_indices(:,n,1)
!         write(83,*)n,remap_wgts(:,n,1)
!      enddo
!stop 'debug'
!return
!!!!!!! 
!k = 223  
      ! now setup the interpolation from coarse to fine-grid vertices
!$omp parallel do private(dotmax,mmax,dotloc,lmax,l,l1,l2,l3,x1,y1,x2,y2,x3,y3,xp,yp,wsum1,denom,m)
      do n = 1,nverticesf
         !    find the nearest coarse grid cell center
         dotmax = -1
         mmax = 0
         DO m = 1,ncellsc
            dotloc = xcell_c(m)*xvert_f(n) + &
                     ycell_c(m)*yvert_f(n) + &
                     zcell_c(m)*zvert_f(n)
            IF(dotloc > dotmax) THEN
               dotmax = dotloc
               mmax = m
            ENDIF
         ENDDO
         remap_indicesv(4,n) = mmax
         
         dotmax = -1
         lmax = 0
         do l = 1,nEdgesOnCell(mmax)
         !    find the nearest vertex belonging to this cell
            dotloc = xvert_c(verticesOnCell(l,mmax))*xvert_f(n) + &
                     yvert_c(verticesOnCell(l,mmax))*yvert_f(n) + &
                     zvert_c(verticesOnCell(l,mmax))*zvert_f(n)
            IF(dotloc > dotmax) THEN
               dotmax = dotloc
               lmax = l
            ENDIF
         enddo
         remap_indicesv(1,n) = verticesOnCell(lmax,mmax)
         if(lmax==nEdgesOnCell(mmax)) then
            remap_indicesv(2,n) = verticesOnCell(1,mmax)
         else
            remap_indicesv(2,n) = verticesOnCell(lmax+1,mmax)
         endif
         if(lmax==1) then
            remap_indicesv(3,n) = verticesOnCell(nEdgesOnCell(mmax),mmax)
         else
            remap_indicesv(3,n) = verticesOnCell(lmax-1,mmax)
         endif
         
         ! find if the fine grid vertex lies within triangle formed by 
         !   indices 1,2,3, or 2,3,4
         ! We will do this by converting to barycentric coordinates.
                     ! compute a triangle in 2D for barycentric coordinates for
                     !   indices 1,2,3
                     l1 = remap_indicesv(1,n)
                     l2 = remap_indicesv(2,n)
                     l3 = remap_indicesv(3,n)
!print*,'triangle 123 ',n
!if(n==k) then
!   print*,'l1,l2,l3 ',l1,l2,l3
!   print*,'xyz_c(l1)',xvert_c(l1),yvert_c(l1),zvert_c(l1)
!   print*,'xyz_c(l2)',xvert_c(l2),yvert_c(l2),zvert_c(l2)
!   print*,'xyz_c(l3)',xvert_c(l3),yvert_c(l3),zvert_c(l3)
!   print*,'xyz_f(n)',xvert_f(n),yvert_f(n),zvert_f(n)
!endif
                     CALL compute_triangle(            &
                             (/xvert_c(l1),yvert_c(l1),zvert_c(l1)/) ,  &
                             (/xvert_c(l2),yvert_c(l2),zvert_c(l2)/) ,  &
                             (/xvert_c(l3),yvert_c(l3),zvert_c(l3)/) ,  &
                             (/xvert_f(n),yvert_f(n),zvert_f(n)/) ,     &
                             x1, y1, x2, y2, x3, y3, xp, yp )
                     ! now that I have a triangle in a 2D plane I can use the 
                     !   formulae in https://codeplea.com/triangular-interpolation
                     !   for the weights
                     denom = 1./ ( (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3) )
                     
                     remap_wgtsv(1,n) = denom * ( (y2-y3)*(xp-x3) + (x3-x2)*(yp-y3) )
                     remap_wgtsv(2,n) = denom * ( (y3-y1)*(xp-x3) + (x1-x3)*(yp-y3) )
                     remap_wgtsv(3,n) = 1. - (remap_wgtsv(1,n)+remap_wgtsv(2,n))
                     remap_wgtsv(4,n) = 0.0
                     remap_indicesv(4,n) = 0
                     
                     wsum1 = sum(abs(remap_wgtsv(:,n)))
!if(n==k) then
!   print*,'x1,y1',x1,y1
!   print*,'x2,y2',x2,y2
!   print*,'x3,y3',x3,y3
!   print*,'xp,yp',xp,yp
!   print*,'remap_indicesv',remap_indicesv(:,n)
!   print*,'remap_wgtsv',remap_wgtsv(:,n)
!endif
                     ! check if sum of abs(weights) > 1 (this means one weight is negative)
                     IF(wsum1 > 1.000001) THEN
                        ! the fine grid vertex is in triangle 2,3,4
                        !   convert this one to barycentric coordinates
                        remap_indicesv(4,n) = mmax
                        remap_wgtsv(1,n) = 0.0
                        l1 = remap_indicesv(2,n)
                        l2 = remap_indicesv(3,n)
                        l3 = remap_indicesv(4,n)
!if(n==k) then
!   print*,'l1,l2,l3 ',l1,l2,l3
!   print*,'xyz_c(l1)',xvert_c(l1),yvert_c(l1),zvert_c(l1)
!   print*,'xyz_c(l2)',xvert_c(l2),yvert_c(l2),zvert_c(l2)
!   print*,'xyz_c(l3)',xcell_c(l3),ycell_c(l3),zcell_c(l3)
!   print*,'xyz_f(n)',xvert_f(n),yvert_f(n),zvert_f(n)
!endif
                        CALL compute_triangle(            &
                                (/xvert_c(l1),yvert_c(l1),zvert_c(l1)/) ,  &
                                (/xvert_c(l2),yvert_c(l2),zvert_c(l2)/) ,  &
                                (/xcell_c(l3),ycell_c(l3),zcell_c(l3)/) ,  &
                                (/xvert_f(n),yvert_f(n),zvert_f(n)/) ,     &
                                x1, y1, x2, y2, x3, y3, xp, yp )
                        ! now that I have a triangle in a 2D plane I can use the 
                        !   formulae in https://codeplea.com/triangular-interpolation
                        !   for the weights
                        denom = 1./ ( (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3) )
                     
                        remap_wgtsv(2,n) = denom * ( (y2-y3)*(xp-x3) + (x3-x2)*(yp-y3) )
                        remap_wgtsv(3,n) = denom * ( (y3-y1)*(xp-x3) + (x1-x3)*(yp-y3) )
                        remap_wgtsv(4,n) = 1. - (remap_wgtsv(2,n)+remap_wgtsv(3,n))

                        wsum = sum(abs(remap_wgtsv(:,n)))
                        ! check if sum of abs(weights) > 1 (this means some assumption earlier failed)
                        IF(wsum > 1.000001) THEN
                           ! copy value from nearest coarse grid vertex
                           dotmax = -1.
                           lmax = 0
                           remap_wgtsv(1,n) = 1.0
                           remap_wgtsv(2:4,n) = 0.0
                           remap_indicesv(4,n) = 0
                           do m = 1,nverticesc
                              dotloc = xvert_c(m)*xvert_f(n) + &
                                       yvert_c(m)*yvert_f(n) + &
                                       zvert_c(m)*zvert_f(n)
                              IF(dotloc > dotmax) THEN
                                 dotmax = dotloc
                                 lmax = m
                              ENDIF
                           enddo
                           remap_indicesv(1,n) = lmax
!                           print*,' nearest neighbor ',n, dotmax
                        endif
                     ENDIF
!if(n==k) then
!   print*,'x1,y1',x1,y1
!   print*,'x2,y2',x2,y2
!   print*,'x3,y3',x3,y3
!   print*,'xp,yp',xp,yp
!   print*,'remap_indicesv',remap_indicesv(:,n)
!   print*,'remap_wgtsv',remap_wgtsv(:,n)
!   print*,'wsum,wsum1',wsum,wsum1
!stop
!endif
      enddo
!$OMP End parallel do
      if(minval(remap_indicesv) < 0.0) then
         print*,'minwgtv ',minval(remap_indicesv)
!$OMP parallel do Private(l, wsum)
         do n = 1,ncellsf
            if(minval(remap_wgtsv(:,n))<0.0) then
               wsum = 0.0
               do l = 1,3
                  remap_wgtsv(l,n) = max(0.0,remap_wgtsv(l,n))
                  wsum = wsum + remap_wgtsv(l,n)
               enddo
               remap_wgtsv(:,n) = remap_wgtsv(:,n) / wsum
            endif
         enddo

!$OMP End parallel do
      endif
      
      END SUBROUTINE initialize_interpolate_variable

!------------------------------------------

      SUBROUTINE compute_triangle(             &
                             v1, v2, v3, vp ,  &
                             x1, y1, x2, y2, x3, y3, xp, yp )
      ! This routine transforms the triangle formed by the cells neighboring
      !  a vertex to 2D space so that we can then compute the barycentric 
      !  coordinates to get interpolation weights.
      
      ! Argument list
      REAL*8, INTENT(IN), DIMENSION(3) :: v1, v2, v3 ! cell locations on coarse grid, 
                                                     !   now vertices of our triangle
      REAL*8, INTENT(IN), DIMENSION(3) :: vp ! cell location of fine grid cell, assumed
                                            ! within the above triangle 
      REAL*8, INTENT(OUT) :: x1, y1, x2, y2, x3, y3  ! new 2D locations of v1, v2, v3 above
      REAL*8, INTENT(OUT) :: xp, yp  ! new 2D location of vp above
      
      ! Local variables
      REAL*8 :: sdist12, sdist23, sdist31, sdist1p, sdist2p ! arc distances
      REAL*8 :: dist12, dist23, dist31, dist1p, dist2p      ! distances in the plane
      REAL*8 :: b1, b2, h, g1, g2, q
      
      ! We will assume vertices proceed counter-clockwise
      !spherical distance is arcsin(dotproduct(v_a,v_b))
      !  this is also the subtended angle on the unit sphere
      sdist12 = acos(max(-1.0,min(1.0,dot_product(v1,v2))))
      sdist23 = acos(max(-1.0,min(1.0,dot_product(v2,v3))))
      sdist31 = acos(max(-1.0,min(1.0,dot_product(v3,v1))))
      sdist1p = acos(max(-1.0,min(1.0,dot_product(v1,vp))))
      sdist2p = acos(max(-1.0,min(1.0,dot_product(v2,vp))))
      
      dist12 = 2. * sin(sdist12*0.5)
      dist23 = 2. * sin(sdist23*0.5)
      dist31 = 2. * sin(sdist31*0.5)
      dist1p = 2. * sin(sdist1p*0.5)
      dist2p = 2. * sin(sdist2p*0.5)
                                        
      x1 = 0.0; y1 = 0.0 ! we place the first vertex here   
      
      b1 = (dist31**2 + dist12**2 - dist23**2 ) / (2.*dist12) 
      h = sqrt(max(0.,dist31**2 - b1**2))
               
      x2 = dist12; y2 = 0.
      x3 = b1; y3 = h        
       
      g1 = (dist1p**2 + dist12**2 - dist2p**2 ) / (2.*dist12) 
      q = sqrt(max(0.,dist1p**2 - g1**2))
               
      xp = g1; yp = q

      END SUBROUTINE compute_triangle

!------------------------------------------

      SUBROUTINE ll2xyz(lon,lat,x,y,z)
      ! This routine takes input lon and lat (in degrees) and returns
      !   the cartesian location on the unit sphere.

      REAL*8, INTENT(IN) :: lon, lat
      REAL*8, INTENT(OUT) :: x, y, z

      REAL*8 :: coslat, lonloc, latloc

      lonloc = lon * d2r
      latloc = lat * d2r
      coslat = cos(latloc)

      x = cos(lonloc) * coslat
      y = sin(lonloc) * coslat
      z = sin(latloc)
 
      END SUBROUTINE ll2xyz
      

!------------------------------------------
   
end program build_ocean_ice_interpolation
