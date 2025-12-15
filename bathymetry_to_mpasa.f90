PROGRAM bathymetry_to_mpasa
! This program will take the 15 second bathymetry data and average it to a quasi-uniform (QU) mpas-atmosphere
!   grid. A negative result is the bathymetry at an ocean point, zero is a land point.
! The program will prompt the user to select a target QU grid

      USE netcdf

   IMPLICIT NONE

      ! kind variables
      INTEGER, PARAMETER :: real8 = SELECTED_REAL_KIND(12) ! 8 byte real
      INTEGER, PARAMETER :: real4 = SELECTED_REAL_KIND( 6) ! 4 byte real

      CHARACTER(len=128) :: meshfile, outfile

      ! variables related to the bathymetry input file
      REAL(real8), DIMENSION(86401) :: bath_tem
      INTEGER*2, DIMENSION(:,:), ALLOCATABLE :: bath_in
      REAL(real8), DIMENSION(86401) :: lon_b
      REAL(real8), DIMENSION(43201) :: lat_b

      ! variables related to the mpasa grid input file
      REAL(real8), DIMENSION(:), ALLOCATABLE :: lon_a, lat_a, x_a, y_a, z_a, bath_a
      REAL(real8), DIMENSION(:), ALLOCATABLE :: lonv_a, latv_a
      REAL(real8), DIMENSION(:), ALLOCATABLE :: x_n, y_n, z_n
      INTEGER, DIMENSION(:), ALLOCATABLE :: nedgesoncell, flag
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: cellsoncell, verticesoncell
      INTEGER :: ncells, nvertices, maxedges

      ! variables for netcdf use
      INTEGER :: ncid, varid, rc, dimid, ncid2, dimidi

      INTEGER :: caseres
      INTEGER :: n, i, j, m, nedges, imax, imin, jmax, jmin, np, count1, &
                 ntot, ncount, imin1, imax1, item, iterh, iter1, iterg,  &
                 n1, n2, nland1, nland2, fillh, fill1, itero
      REAL(real8) :: mag, lon, dist
      REAL(real8), PARAMETER :: r2d = 180. / (4.*ATAN(1.)), d2r = 1./r2d
      REAL(real8) :: lonmin, lonmax, latmin, latmax

       integer :: testindx
       
! Begin executable code

      PRINT *,'   Enter the global mesh (mpasa) filename '
      READ(*,'(a128)') meshfile
      ! read from the input mpasa grid file
      rc = nf90_open(trim(meshfile), nf90_nowrite,ncid)
      if(rc .ne. 0) stop ' troubles opening the global mesh (mpasa) file - stopping'

      PRINT *,'   Enter the new mpasa bathymetry filename '
      READ(*,'(a128)') outfile
      rc = nf90_create(trim(outfile),nf90_clobber,ncid2)

      rc = nf90_inq_dimid(ncid,'nCells',dimid)
      rc = nf90_inquire_dimension(ncid,dimid,len=ncells)
      ALLOCATE(bath_a(ncells))
      ALLOCATE(lon_a(ncells))
      ALLOCATE(lat_a(ncells))
      ALLOCATE(  x_a(ncells))
      ALLOCATE(  y_a(ncells))
      ALLOCATE(  z_a(ncells))
      rc = nf90_inq_varid(ncid,'lonCell',varid)
      rc = nf90_get_var(ncid,varid, lon_a)
      rc = nf90_inq_varid(ncid,'latCell',varid)
      rc = nf90_get_var(ncid,varid, lat_a)
      rc = nf90_inq_varid(ncid,'xCell',varid)
      rc = nf90_get_var(ncid,varid,   x_a)
      rc = nf90_inq_varid(ncid,'yCell',varid)
      rc = nf90_get_var(ncid,varid,   y_a)
      rc = nf90_inq_varid(ncid,'zCell',varid)
      rc = nf90_get_var(ncid,varid,   z_a)
      
      rc = nf90_inq_dimid(ncid,'nVertices',dimid)
      rc = nf90_inquire_dimension(ncid,dimid,len=nVertices)
      ALLOCATE(lonv_a(nVertices))
      ALLOCATE(latv_a(nVertices))
      rc = nf90_inq_varid(ncid,'lonVertex',varid)
      rc = nf90_get_var(ncid,varid, lonv_a)
      rc = nf90_inq_varid(ncid,'latVertex',varid)
      rc = nf90_get_var(ncid,varid, latv_a)

      rc = nf90_inq_dimid(ncid,'maxEdges',dimid)
      rc = nf90_inquire_dimension(ncid,dimid,len=maxedges)
      ALLOCATE(cellsoncell(maxedges,ncells))
      ALLOCATE(verticesoncell(maxedges,ncells))
      ALLOCATE(nedgesoncell(ncells))
      ALLOCATE(x_n(maxedges))
      ALLOCATE(y_n(maxedges))
      ALLOCATE(z_n(maxedges))
      rc = nf90_inq_varid(ncid,'cellsOnCell',varid)
      DO n = 1,ncells,100000
      print*,'cellsOnCell block ',n,' of ', (ncells-1)/100000+1
         n1 = MIN(n+99999,ncells)
         n2 = 1+n1 - n
         rc = nf90_get_var(ncid,varid, cellsoncell(:,n:n1),(/1,n/),(/maxedges,n2/))
      ENDDO
      rc = nf90_inq_varid(ncid,'nEdgesOnCell',varid)
      DO n = 1,ncells,100000
      print*,'nEdgesOnCell block ',n,' of ', (ncells-1)/100000+1
         n1 = MIN(n+99999,ncells)
         n2 = 1+n1 - n
         rc = nf90_get_var(ncid,varid, nEdgesOnCell(n:n1),(/n/),(/n2/))
      ENDDO
      rc = nf90_inq_varid(ncid,'verticesOnCell',varid)
      DO n = 1,ncells,100000
      print*,'verticesOnCell block ',n,' of ', (ncells-1)/100000+1
         n1 = MIN(n+99999,ncells)
         n2 = 1+n1 - n
         rc = nf90_get_var(ncid,varid, verticesoncell(:,n:n1),(/1,n/),(/maxedges,n2/))
      ENDDO

      rc = nf90_close(ncid)

      ! convert lat/lon to degress, compute cartesian coordinates, normalize cartesian coordinates
      DO n = 1,ncells
         lon_a(n) = lon_a(n) * r2d
         lat_a(n) = lat_a(n) * r2d
      ENDDO
      DO n = 1,nvertices
         lonv_a(n) = lonv_a(n) * r2d
         latv_a(n) = latv_a(n) * r2d
      ENDDO
      DO n = 1,ncells
         mag = 1./SQRT(x_a(n)**2 + y_a(n)**2 + z_a(n)**2)
         x_a(n) = x_a(n) * mag
         y_a(n) = y_a(n) * mag
         z_a(n) = z_a(n) * mag
      ENDDO

      testindx = 0
            print*,'testindx0',testindx
      do n = 1,ncells
         if( abs(lon_a(n)-70.682)<0.05 .and. abs(lat_a(n)-72.998)<0.05 )then
            testindx = n
            print*,'testindx',testindx
         endif
      enddo

      ! read from the bathymetry input file
      rc = nf90_open('./inputdata/SRTM15_plus_earth_relief_15s_original.nc',  &
                     nf90_nowrite,ncid)
      if(rc .ne. 0) stop ' troubles opening the 15s input bathymetry file - stopping'

      rc = nf90_inq_varid(ncid,'lon',varid)
      rc = nf90_get_var(ncid,varid,lon_b)
      rc = nf90_inq_varid(ncid,'lat',varid)
      rc = nf90_get_var(ncid,varid,lat_b)
      rc = nf90_inq_varid(ncid,'z',varid)
      ALLOCATE(bath_in(86401,43201))
      DO j = 1,43201
         if(mod(j,500)==1)print*,'reading lat ',j
         rc = nf90_get_var(ncid,varid,bath_tem,(/1,j/),(/86401,1/))
         bath_in(:,j) = bath_tem
      ENDDO
      rc = nf90_close(ncid)

      ! Match the range of longitudes between the input grids
      DO n = 1,ncells
         IF(lon_a(n) >=180.) lon_a(n) = lon_a(n) - 360.
      ENDDO
      DO n = 1,nvertices
         IF(lonv_a(n) >=180.) lonv_a(n) = lonv_a(n) - 360.
      ENDDO

      np = ncells + 1
      ! compute the bathymetry for each mpas-a grid cell
!$OMP Parallel DO Private(nedges, m, x_n, y_n, z_n, imin, imax, jmin, jmax, latmin, latmax, lonmin, lonmax, lon, imax1, imin1, ncount, ntot)
      DO n = 1,ncells
         
         nedges = nedgesoncell(n)
         latmin = lat_a(n)
         latmax = lat_a(n)
         lonmin = lon_a(n)
         lonmax = lon_a(n)
         DO m = 1,nedges
            x_n(m) = x_a(cellsoncell(m,n))
            y_n(m) = y_a(cellsoncell(m,n))
            z_n(m) = z_a(cellsoncell(m,n))
            latmin = MIN(latmin,latv_a(verticesoncell(m,n)))
            latmax = MAX(latmax,latv_a(verticesoncell(m,n)))
            lonmin = MIN(lonmin,lonv_a(verticesoncell(m,n)))
            lonmax = MAX(lonmax,lonv_a(verticesoncell(m,n)))
         ENDDO
         imin = MAX(INT((lonmin+180)*240) + 1,1)
         imax = MIN(INT((lonmax+180)*240) + 1,86400)
         jmin = MAX(INT((latmin+90)*240) + 1,1)
         jmax = MIN(INT((latmax+90)*240) + 1,43201)
         
         ntot = 0
         ncount = 0
         bath_a(n) = 0
         IF(jmin < 3 ) THEN
            if(n==3)print*,'n=3 south pole adjustment'
            !south pole adjustment
            bath_a(n) = 0.0
            print*,'south pole ',n
         ELSE IF (jmax > 43199) THEN
            !north pole adjustment
            bath_a(n) = sum(real(bath_in(1:86400,jmin:jmax),real8)) / (86400*(1+jmax-jmin))
            np = n
            print*,'north pole ',np, sum(bath_in(1:86400,jmin:jmax)), (86400*(1+jmax-jmin))
            print*,'north pole ',bath_a(n),jmin,jmax
         ELSE IF(lonmax-lonmin > 180.) THEN
            if(n==3)print*,'n=3 straddle date line'
            !straddling date line adjustment
            imin = 1
            imax = 86400

            imax1 = 1
            imin1 = 86400
            DO m = 1,nedges
               lon = ATAN2 (x_n(m),y_n(m)) * r2d
               IF (lon >= 180 ) lon = lon - 360.
               item = (lon+180.)*240 + 1
               IF (lon < 0) THEN
                  imax1 = MAX(imax1, item)
               ELSE
                  imin1 = MIN(imin1, item)
               ENDIF
            ENDDO
            
            CALL compute_bathymetry(nedges, imin, imax1, jmin, jmax, &
                                    x_a(n),y_a(n),z_a(n),       &
                                    x_n(1:nedges),y_n(1:nedges),z_n(1:nedges), &
                                    lon_b, lat_b,                &
                                    !x_b(imin:imax1,jmin:jmax),   &
                                    !y_b(imin:imax1,jmin:jmax),   &
                                    !z_b(imin:imax1,jmin:jmax),   &
                                    bath_in(imin:imax1,jmin:jmax), bath_a(n), &
                                    ntot, ncount )
            CALL compute_bathymetry(nedges, imin1, imax, jmin, jmax, &
                                    x_a(n),y_a(n),z_a(n),       &
                                    x_n(1:nedges),y_n(1:nedges),z_n(1:nedges), &
                                    lon_b, lat_b,                &
                                    !x_b(imin1:imax,jmin:jmax),   &
                                    !y_b(imin1:imax,jmin:jmax),   &
                                    !z_b(imin1:imax,jmin:jmax),   &
                                    bath_in(imin1:imax,jmin:jmax), bath_a(n), &
                                    ntot, ncount )
         ELSE
            if(n==3)print*,'n=3 default'
            CALL compute_bathymetry(nedges, imin, imax, jmin, jmax, &
                                    x_a(n),y_a(n),z_a(n),       &
                                    x_n(1:nedges),y_n(1:nedges),z_n(1:nedges), &
                                    lon_b, lat_b,                &
                                    !x_b(imin:imax,jmin:jmax),   &
                                    !y_b(imin:imax,jmin:jmax),   &
                                    !z_b(imin:imax,jmin:jmax),   &
                                    bath_in(imin:imax,jmin:jmax), bath_a(n), &
                                    ntot, ncount )
         ENDIF
         if(n .ne. np) then
         IF(ncount > ntot/2) THEN
            if(n==3)print*,'n=3 ncount large'
            bath_a(n) = bath_a(n) / REAL(ncount,real8)
         ELSE
            if(n==3)print*,'n=3 ncount small'
            bath_a(n) = 0.0
         ENDIF
         endif
         IF(MOD(n,1000)==1) PRINT*,'n,bathy ',n,bath_a(n)
         if(n==3)print*,'n=3',np,bath_a(n)
      ENDDO
      print*,'northpole 0',np, bath_a(np)

! initial bathymetry_to_mpasa
      itero = 1
      ! define and write out the mpas-a bathymetry_to_mpasa
      rc = nf90_def_dim(ncid2,'nCells',ncells,dimid)
      rc = nf90_def_dim(ncid2,'iter',nf90_unlimited,dimidi)
      rc = nf90_def_var(ncid2,'bathymetry',nf90_double,(/dimid/),varid)
         rc = nf90_put_att(ncid2,varid,'units','m')
         rc = nf90_put_att(ncid2,varid,'long_name','ocean bathymetry')
      rc = nf90_enddef(ncid2)
      !rc = nf90_put_var(ncid2,varid,bath_a,(/1,itero/),(/ncells,1/))
      
      !print*,'init bathy(2341120)',bath_a(2341120)
      !print*,'init bathyn(2341120)',bath_a(cellsoncell(1:6,2341120))

      
! This block to be iterated as necessary
      iterg = 0
      fill1 = 1
      fillh = 1
      ALLOCATE(flag(ncells))   
      

print*,'northpole1',bath_a(np)
   DO WHILE (fillh + fill1 > 0)
      iterg = iterg + 1
      
      ! keep only the main contiguous ocean (mco) points 

      ! flag is to keep track of things
      ! -1 is land
      !  0 is provisionally ocean
      !  1 is an ocean neighbor of a proven cell and whose neighbors will
      !    be identified next iteration
      !  2 is an ocean cell that is part of the mco
      !  We start by assigning a priori the above identified North Pole flag(np) = 1
      !  In an iterative loop we:
      !    Loop over all cells. In this loop
      !      If we find a flag=1, change it to 2
      !      Reassign all neighbors with flag=0 to 1
      !    End cell loop
      !    If the count of flag=1 cells is zero exit iterative loop
      !  Set bathymetry of cells with flag = 0 to 0.0
      !  This procedure should eliminate all inland depressions or seas 
      !    that are not contiguous with the main ocean (eg Death Valley or Dead Sea)
      
      nland1 = 0
      DO n = 1,ncells
         IF(bath_a(n) > -0.5) THEN
            flag(n) = -1  ! already land
            nland1 = nland1 + 1
         ELSE
            flag(n) = 0   ! provisionally ocean
         ENDIF
      ENDDO
      count1 = 1
      flag(np) = 1
      iterh = 0
      DO WHILE (count1>0) 
         iterh = iterh + 1
         DO n = 1,ncells
            IF(flag(n)==1) THEN
               flag(n) = 2
               count1 = count1 - 1
               DO m = 1,nedgesoncell(n)
                  IF(flag(cellsoncell(m,n))==0)THEN
                     flag(cellsoncell(m,n)) = 1
                     count1 = count1 + 1
                  ENDIF
             ENDDO
            ENDIF
         ENDDO
      ENDDO
      nland2 = 0
      DO n = 1,ncells
         IF(flag(n)==0) THEN 
            bath_a(n) = 0.0
            flag(n)=-1
         ENDIF
         IF(flag(n) < 0) nland2 = nland2 + 1
      ENDDO
      PRINT*,'hole fill iter = ',iterh
      fillh = nland2 - nland1
      PRINT*,'hole number filled = ',fillh

      itero = itero + 1
      
      print*,'northpole2',bath_a(np)
      
      ! Now fill cells that have only one edge in contact with other cells
      fill1 = 0
      DO n = 1,ncells
         IF(flag(n)==2) THEN
            count1 = 0
            DO m = 1,nedgesoncell(n)
               IF(flag(cellsoncell(m,n))==2)THEN
                  count1 = count1 + 1
               ENDIF
            ENDDO
            IF(count1 == 1) THEN
               fill1 = fill1 + 1
               flag(n) = -1
               bath_a(n) = 0.0
            ENDIF
         ENDIF
      ENDDO
      
      PRINT*,'inlet number filled = ',fill1
      itero = itero + 1
      
      ! We have to iterate for we may have created cells no longer part of the mco.

   ENDDO
     
   print*,'northpole3',bath_a(np)
      PRINT*,'global iter = ',iterg

!!15km - fill/dredge ice-problematic spots
!     print*,'15km adjustment'
!     print *,'lon ',minval(lon_a),maxval(lon_a)
!     print *,'lat ',minval(lat_a),maxval(lat_a)
!      do n = 1,ncells
!         if ( lon_a(n) > 76.7 .and. lon_a(n) < 77.9 .and.  &
!              lat_a(n) > 68.1 .and. lat_a(n) < 69.1 ) then             
!            print*,'15km fill 1 ',n
!             bath_a(n) = 0.0
!         endif
!         if ( lon_a(n) > 77.4 .and. lon_a(n) < 78.3 .and.  &
!              lat_a(n) > 71.0 .and. lat_a(n) < 71.3 ) then
!            print*,'15km fill 2 ',n
!             bath_a(n) = 0.0
!         endif
!         if ( lon_a(n) > 70.4 .and. lon_a(n) < 70.5 .and.  &
!              lat_a(n) > 72.9 .and. lat_a(n) < 73.0) then
!            print*,'15km fill 4 ',n
!             bath_a(n) = 0.0
!         endif
!         if ( lon_a(n) > 70.0 .and. lon_a(n) < 70.1 .and.  &
!              lat_a(n) > 72.9 .and. lat_a(n) < 73.0) then
!            print*,'15km fill 5 ',n
!             bath_a(n) = 0.0
!         endif
!         if ( lon_a(n) > -91.7 .and. lon_a(n) < -80.1 .and.  &
!              lat_a(n) > 80.0 .and. lat_a(n) < 81.5 ) then
!            print*,'15km fill 3 ',n
!             bath_a(n) = 0.0
!         endif
!         if ( lon_a(n) > -80.0 .and. lon_a(n) < -79.9 .and.  &
!              lat_a(n) > 80.65 .and. lat_a(n) < 80.75) then
!            print*,'15km fill 6 ',n
!             bath_a(n) = 0.0
!         endif
!         if ( lon_a(n) > -93.3 .and. lon_a(n) < -93.2 .and.  &
!              lat_a(n) > 81.1 .and. lat_a(n) < 81.2 ) then
!            print*,'15km fill 7 ',n
!             bath_a(n) = 0.0
!         endif
!         dist = (lon_a(n)+93.23)**2
!         dist = dist + (lat_a(n)-81.15)**2
!         if( sqrt(dist) < 0.1 ) then
!           print*,'filled desired spot',n
!           bath_a(n) = -20.
!         endif
!      enddo

! write the output
      rc = nf90_put_var(ncid2,varid,bath_a)
      print*,'put rc ',rc
     
! End iterated block.
         
      rc = nf90_close(ncid2)

   CONTAINS

!------------------------------------------

      SUBROUTINE ll2xyz(lon,lat,x,y,z)
      ! This routine takes input lon and lat (in degrees) and returns
      !   the cartesian location on the unit sphere.

      REAL(real8), INTENT(IN) :: lon, lat
      REAL(real8), INTENT(OUT) :: x, y, z

      REAL(real8) :: coslat, lonloc, latloc

      lonloc = lon * d2r
      latloc = lat * d2r
      coslat = cos(latloc)

      x = cos(lonloc) * coslat
      y = sin(lonloc) * coslat
      z = sin(latloc)
 
      END SUBROUTINE ll2xyz

!------------------------------------------

      SUBROUTINE compute_bathymetry(nedges, imin, imax, jmin, jmax,  &
                                    x_a,y_a,z_a,   &
                                    x_n,y_n,z_n,   &
                                    lon_b, lat_b,  &
!                                    x_b,y_b,z_b,   &
                                    bath_in, bath_a, ntot, ncount )

      ! Compute the bathymetry at an mpas-a cell
      !  Approximations:
      !   the mpas-a cell bathymetry is an unweighted average of all input points that lie
      !    closest to its center.

      ! argument list
      INTEGER, INTENT(IN) :: nedges   ! number of edges (and vertices and neighboring cells)
      INTEGER, INTENT(IN) :: imin, imax, jmin, jmax  ! the bounds of the input bathymetry arrays
      REAL(real8), INTENT(IN) :: x_a, y_a, z_a  ! cartesian locations of mpas-s cell
      REAL(real8), INTENT(IN), DIMENSION(nedges) :: x_n, y_n, z_n  ! cartesian locations of mpas-s cell neighbors
      REAL(real8), INTENT(IN), DIMENSION(86401) :: lon_b
      REAL(real8), INTENT(IN), DIMENSION(43201) :: lat_b
      INTEGER*2, INTENT(IN), DIMENSION(imin:imax,jmin:jmax) :: bath_in
      REAL(real8), INTENT(INOUT) :: bath_a
      INTEGER, INTENT(INOUT) :: ncount, ntot

      ! local variables
      REAL(real8) :: dot_a, dot_n
      REAL(real8) :: x_b, y_b, z_b  ! cartesian locations of input bathymetry     
      INTEGER i, j, m
      LOGICAL dotloc

         DO j = jmin,jmax
            DO i = imin,imax
               CALL ll2xyz(lon_b(i),lat_b(j),x_b,y_b,z_b)
               dot_a = x_a*x_b + y_a*y_b + z_a*z_b
               DO m = 1,nedges
                  dot_n = x_n(m)*x_b + y_n(m)*y_b + z_n(m)*z_b
                  dotloc = dot_n >= dot_a
                  IF(dotloc) EXIT
               ENDDO
               IF(dotloc) CYCLE
               ntot = ntot + 1 !total input points in cell
               IF(bath_in(i,j) > -1.) CYCLE
               ncount = ncount + 1 ! total ocean points in cell
               bath_a = bath_a + bath_in(i,j)
            ENDDO
         ENDDO
      
      END SUBROUTINE compute_bathymetry
!------------------------------------------

END PROGRAM bathymetry_to_mpasa
