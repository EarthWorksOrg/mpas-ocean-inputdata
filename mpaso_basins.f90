program mpaso_basins

use netcdf

implicit none

! This program will construct a basins file, using the template file to 
!    supply information about the groups, and the mpaso grid file for grid info.

! compile with netcdf library


   character*128 infile,     & ! mpaso ocean mesh file
                 outfile       ! basin file to write
                 
   ! netcdf variables
   integer ncidin, ncidout, rc, dimid, varid
   ! specific to output file
   integer dimidc, dimide, dimidv, dimidt, dimidr, varidout
   
   ! dimension sizes
   integer ncells, nvertices, nedges, maxedges, maxvertices,   &
           nregions, nstring, ntransects, ntem, edge1, vertex1
   
   
   ! mesh locations
   real*8, dimension(:), allocatable ::  lonc, latc, lontem, xc, yc, zc
   
   ! transect nominal start and final points
   real*8 :: lonw, lone, lats, latn, latp, lonp
   
   ! computed transect/basin files
   integer, dimension(:,:), allocatable :: regionCellMasks, cellsoncell
   integer, dimension(:), allocatable :: nedgesoncell
   
   ! misc. variables
   real*8 r2d, cmag
   integer n
   
!!! BEGIN EXECUTABLE CODE
   
   print*,'enter the ocean mesh file name'
   read(5,*)infile
   print*,'enter the ocean basin and transects file name to be created'
   read(5,*)outfile
   !infile = '/Volumes/disk11/mpas-ocean/poutfiles/ocean.QU.015km.250922.nc'
   !outfile = '15km/QU015_mocBasinsAndTransects20250922.nc'
   
   rc = nf90_open(trim(infile),nf90_nowrite,ncidin)
   print*,'open mesh file ',rc, trim(infile)
   rc = nf90_open(trim(outfile),nf90_write,ncidout)
   print*,'open transect file ',rc, trim(outfile)

   ! dimensions from outfile
   rc = nf90_inq_dimid(ncidout,'nRegions',dimidr)
   rc = nf90_inquire_dimension(ncidout,dimidr,len=nregions)
   
   if(nregions .ne.5) stop 'we expected 5 regions'
      
   ! dimensions from infile
   rc = nf90_inq_dimid(ncidin,'nCells',dimid)
   rc = nf90_inquire_dimension(ncidin,dimid,len=ncells)
   

   rc = nf90_redef(ncidout)
   rc = nf90_def_dim(ncidout,'nCells',ncells,dimidc)

   ! read mesh variables from infile   
   allocate(lonc(ncells))
   allocate(lontem(ncells))
   allocate(latc(ncells))
   allocate(xc(ncells))
   allocate(yc(ncells))
   allocate(zc(ncells))
   allocate(nedgesoncell(ncells))
   allocate(cellsoncell(6,ncells))
   
   rc = nf90_inq_varid(ncidin,'nEdgesOnCell',varid)
   rc = nf90_get_var(ncidin,varid,nedgesOnCell)
   rc = nf90_inq_varid(ncidin,'cellsOnCell',varid)
   rc = nf90_get_var(ncidin,varid,cellsoncell)

   rc = nf90_inq_varid(ncidin,'lonCell',varid)
   rc = nf90_get_var(ncidin,varid,lonc)
   rc = nf90_inq_varid(ncidin,'latCell',varid)
   rc = nf90_get_var(ncidin,varid,latc)

   rc = nf90_inq_varid(ncidin,'xCell',varid)
   rc = nf90_get_var(ncidin,varid,xc)
   rc = nf90_inq_varid(ncidin,'yCell',varid)
   rc = nf90_get_var(ncidin,varid,yc)
   rc = nf90_inq_varid(ncidin,'zCell',varid)
   rc = nf90_get_var(ncidin,varid,zc)
   
   ! convert radians to degrees
   r2d = 180./ (4.*atan(1.))
   lonc = lonc * r2d
   latc = latc * r2d
   
   where (lonc < 0.0) lonc = lonc + 360

   ! normalize cartesian locations to the unit sphere
   do n = 1,ncells
      cmag = sqrt(xc(n)**2+yc(n)**2+zc(n)**2)
      xc(n) = xc(n) / cmag
      yc(n) = yc(n) / cmag
      zc(n) = zc(n) / cmag
   enddo
   
   
   allocate(regionCellMasks(nregions,ncells))
   regionCellMasks = 0
   
   ! define more transect variables
   rc = nf90_def_var(ncidout,'regionCellMasks',nf90_int,(/dimidr,dimidc/),varid)
   
   rc = nf90_enddef(ncidout)
   
   ! here is the bulk of the work, region by region
   ! specify longitude bounds so that lone > lonw
   ! for each region, first we do a rough mask, then fine tune
   
    ! Atlantic_MOC",
   n = 1
   lontem = lonc
   where(lonc > 180.)lontem = lontem - 360.
   
   ! fill blocks
   lonw = -60.
   lone = 19.
   lats = -34.
   latn = 10.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -60.
   lone = -6.
   lats = 5.
   latn = 60.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -60.
   lone = 0.
   lats = 41.
   latn = 60.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -100.
   lone = 32.
   lats = 48.
   latn = 67.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -130.
   lone = 25.
   lats = 66.
   latn = 71.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -129.
   lone = 24.
   lats = 71.
   latn = 72.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -129.
   lone = 22.
   lats = 72.
   latn = 73.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -129.
   lone = 20.
   lats = 73.
   latn = 74.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -129.
   lone = 18.
   lats = 74.
   latn = 75.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -127.
   lone = 16.
   lats = 75.
   latn = 76.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -127.
   lone = 16.
   lats = 76.
   latn = 77.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -125.
   lone = 16.
   lats = 77.
   latn = 78.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -123.
   lone = 16.
   lats = 78.
   latn = 79.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -110.
   lone = -50.
   lats = 79.
   latn = 80.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -105.
   lone = -50.
   lats = 80.
   latn = 81.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -100.
   lone = -50.
   lats = 81.
   latn = 82.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -90.
   lone = -50.
   lats = 82.
   latn = 83.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lonw = -80.
   lone = -50.
   lats = 83.
   latn = 84.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   ! fill basins
      ! west atlantic
   lonp = -70.
   latp = 30.
   call fill_basin(lonp,latp,regionCellMasks(n,:))
      ! gulf of st lawrence
   lonp = -60.
   latp = 47.
   call fill_basin(lonp,latp,regionCellMasks(n,:))
   
   ! AtlanticMed_MOC",
   n = 2
   lats = -34.
   lonw = -60.
   lone = 30.
   latn = 5.
   lontem = lonc
   where(lonc > 180.)lontem = lontem - 360.
   regionCellMasks(2,:) = regionCellMasks(1,:)
   !call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   ! fill the mediterranean
   lonp = 10
   latp = 40.
   call fill_basin(lonp,latp,regionCellMasks(n,:))
   
   ! IndoPacific_MOC",
   n = 3
   lats = -34.
   latn = 9.
   lonw = 26.
   lone = 300.
   lontem = lonc
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lats = 9.
   latn = 24.
   lonw = 43.
   lone = 260.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lats = 52.
   latn = 66.
   lonw = 160.
   lone = 225.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lats = 24.
   latn = 30.
   lonw = 56.
   lone = 90.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   ! fill the n pacific
   lonp = 180.
   latp = 50.
   call fill_basin(lonp,latp,regionCellMasks(n,:))
   ! fill the campeche
   lonp = 265
   latp = 11.
   call fill_basin(lonp,latp,regionCellMasks(n,:))
   ! fill the gulf of californiat
   lonp = 249.5
   latp = 26.5
   call fill_basin(lonp,latp,regionCellMasks(n,:))
   ! unfill the south of australia
   lonp = 130.
   latp = -33.
   call unfill_basin(lonp,latp,regionCellMasks(n,:))
   
   ! Pacific_MOC",
   n = 4
   lats = -6.
   latn = 9.
   lonw = 104.75
   lone = 280.
   lontem = lonc
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lats = 52.
   latn = 66.
   lonw = 160.
   lone = 225.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lats = 9.
   latn = 24.
   lonw = 98.
   lone = 260.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lats = -3.5
   latn = 3.
   lonw = 103.5
   lone = 108.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   ! fill the n pacific
   lonp = 180.
   latp = 50.
   call fill_basin(lonp,latp,regionCellMasks(n,:))
   ! fill the campeche
   lonp = 265
   latp = 11.
   call fill_basin(lonp,latp,regionCellMasks(n,:))
   ! fill the gulf of californiat
   lonp = 249.5
   latp = 26.5
   call fill_basin(lonp,latp,regionCellMasks(n,:))
   ! fill the south of panama
   lonp = 281.
   latp = 5.
   call fill_basin(lonp,latp,regionCellMasks(n,:))
   ! fill along thailand
   lonp = 103.
   latp = 7.
   call fill_basin(lonp,latp,regionCellMasks(n,:))
   
   ! Indian_MOC" ;
   n = 5
   lontem = lonc
   lats = -6.
   latn = -3.
   lonw = 38.
   lone = 104.75
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lats = -3.
   latn = 2.
   lonw = 38.
   lone = 103.5
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lats = 2.
   latn = 14.
   lonw = 43.25
   lone = 50.5
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   lats = 21.
   latn = 28.
   lonw = 56.
   lone = 60.
   call mask_basin(lonw,lone,lats,latn,regionCellMasks(n,:))
   ! fill indian ocean
   lonp = 63.
   latp = 10.
   call fill_basin(lonp,latp,regionCellMasks(n,:))
   
        
   rc = nf90_inq_varid(ncidout,'regionCellMasks',varid)
   rc = nf90_put_var(ncidout,varid,regioncellmasks)
   

   rc = nf90_close(ncidout)

contains

   subroutine mask_basin(lonw,lone,lats,latn,regionCellMasks)
   real*8, intent(in) :: lonw, lone, lats, latn
   integer, intent(out), dimension(ncells) :: regionCellMasks
   
   integer n
   
   do n = 1,ncells
      if(lontem(n) > lonw .and. lontem(n) < lone .and. latc(n) >= lats .and. latc(n) < latn) regionCellMasks(n) = 1
   enddo
   
   end subroutine mask_basin

   subroutine fill_basin(lonp,latp,regionCellMasks)
   real*8, intent(in) :: lonp, latp
   integer, intent(inout), dimension(ncells) :: regionCellMasks
   
   integer n, np, m, k, count
   real*8 xp, yp, zp, dotloc, dotmax
   integer, dimension(:), allocatable :: cstatus
   
   ! we will find the cell nearest lonp and latp, then mask it
   ! then we will iterate and mask unmasked neighbors until there
   ! are no more unmasked neighbors
   
   allocate(cstatus(ncells))

   !convert starting point to cartesian unit sphere
   xp = cos((lonp)/r2d) * cos(latp/r2d)
   yp = sin((lonp)/r2d) * cos(latp/r2d)
   zp =                   sin(latp/r2d)

   dotmax = 0.0
   np = 0
   do n = 1,ncells
      dotloc = xp*xc(n) + yp*yc(n) + zp*zc(n)
      if(dotloc > dotmax) then
         dotmax = dotloc
         np = n
      endif
   enddo
   
   cstatus = regionCellMasks*2
   ! initialize cstatus
   ! cstatus = 0 - unmasked
   ! cstatus = 1 - reset these cells in the next iteration and check their neighbors
   ! cstatus = 2 - mask is set and off bounds for searching
   
   cstatus(np) = 1
   
   count = 1  ! this is the number of cells with cstatus = 1
   do while( count > 0)
      do n = 1,ncells
         if(cstatus(n) ==1) then
            ! we will set cstatus = 0 neighbors to 1 and set this cstatus to 2
            cstatus(n) = 2
            regionCellMasks(n) = 1
            count = count - 1
            do m = 1,nedgesoncell(n)
               k = cellsoncell(m,n)
               if(k > 0) then
                  if(cstatus(k) == 0) then
                     cstatus(k) = 1
                     count = count + 1
                  endif
               endif
            enddo
         endif
      enddo
      !print*,'count ',count
   enddo
   
   deallocate(cstatus)
   end subroutine fill_basin


   subroutine unfill_basin(lonp,latp,regionCellMasks)
   real*8, intent(in) :: lonp, latp
   integer, intent(inout), dimension(ncells) :: regionCellMasks
   
   integer n, np, m, k, count
   real*8 xp, yp, zp, dotloc, dotmax
   integer, dimension(:), allocatable :: cstatus
   
   ! we will find the cell nearest lonp and latp, then unmask it
   ! then we will iterate and unmask unmasked neighbors until there
   ! are no more masked neighbors
   
   allocate(cstatus(ncells))

   !convert starting point to cartesian unit sphere
   xp = cos((lonp)/r2d) * cos(latp/r2d)
   yp = sin((lonp)/r2d) * cos(latp/r2d)
   zp =                   sin(latp/r2d)

   dotmax = 0.0
   np = 0
   do n = 1,ncells
      dotloc = xp*xc(n) + yp*yc(n) + zp*zc(n)
      if(dotloc > dotmax) then
         dotmax = dotloc
         np = n
      endif
   enddo
   
   cstatus = (1-regionCellMasks)*2
   ! initialize cstatus
   ! cstatus = 0 - masked
   ! cstatus = 1 - reset these cells in the next iteration and check their neighbors
   ! cstatus = 2 - mask is set and off bounds for searching
   
   cstatus(np) = 1
   
   count = 1  ! this is the number of cells with cstatus = 1
   do while( count > 0)
      do n = 1,ncells
         if(cstatus(n) ==1) then
            ! we will set cstatus = 0 neighbors to 1 and set this cstatus to 2
            cstatus(n) = 2
            regionCellMasks(n) = 0
            count = count - 1
            do m = 1,nedgesoncell(n)
               k = cellsoncell(m,n)
               if(k > 0) then
                  if(cstatus(k) == 0) then
                     cstatus(k) = 1
                     count = count + 1
                  endif
               endif
            enddo
         endif
      enddo
      !print*,'count ',count
   enddo
   
   deallocate(cstatus)
   end subroutine unfill_basin


end program mpaso_basins
