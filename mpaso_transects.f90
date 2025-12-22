program mpaso_transects

use netcdf

implicit none

! This program will construct a transects file, using the template file to 
!    supply information about the groups, and the mpaso grid file for grid info.

! compile with netcdf library

! number and transects and regions assumed to be identical

! METHOD:
!   Get region and transect info from transect using ncks
!     ncks -x -v transectEdgeMasks,transectEdgeMaskSigns,transectVertexMasks,regionCellMasks,transectVertexGlobalIDs,transectEdgeGlobalIDs \
!        demo/QU60_mocBasinsAndTransects20210623.nc <output_file>
!   Define dimensions -
!     nedges, ncells, and nvertices from mesh (infile)
!     maxedgesintransect, maxverticesintransect after all transects computed
!   Variables -
!     compute transectedgemasks, transectedgemasksigns, regioncellMasks,
!             transectvertexglobalids, transectedgeglobalids
!   For each transect:
!     Find the starting point 
!     Proceed south along coast to transect latitude 
!     Proceed east along transect latitude until land or final point encountered
!     As we proceed we save the transectedgeglobalids and the transectvertexglobalids
!     Compute the mask variables 

   character*128 infile,     & ! mpaso ocean mesh file
                 outfile       ! transect file to write
                 
   integer, parameter :: idlim = 100000
   ! netcdf variables
   integer ncidin, ncidout, rc, dimid, varid
   ! specific to output file
   integer dimidc, dimide, dimidv, dimidt, dimidr, varidout
   
   ! dimension sizes
   integer ncells, nvertices, nedges, maxedges, maxvertices,   &
           nregions, nstring, ntransects, ntem, edge1, vertex1
   integer, dimension(:), allocatable :: vcount, ecount
   
   ! mesh indices
   integer, dimension(:), allocatable :: vertexid, edgeid
   
   ! mesh neighbor variables
   integer, dimension(:,:), allocatable :: edgesOnVertex, verticesOnEdge,   &
            cellsOnCell, cellsOnEdge, cellsOnVertex, edgesoncell
   integer, dimension(:), allocatable :: nedgesOnCell
   
   ! mesh locations
   real*8, dimension(:), allocatable :: lone, late, lonv, latv, lonc, latc, angle
   real*8, dimension(:), allocatable :: xv, yv, zv, xc, yc, zc
   
   ! transect nominal start and final points
   real*8, dimension(:), allocatable :: lon1, lonee, lat1, lonw, translat
   
   ! computed transect/basin files
   integer, dimension(:), allocatable :: transectsInGroup, transectEdgeMasks, &
            regionCellMasks, transectVertexMasks,     &
            temvertex, temedge, transectEdgeMaskSigns
   integer, dimension(:,:), allocatable :: transectEdgeGlobalIDs, transectVertexGlobalIDs
   
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
   rc = nf90_inq_dimid(ncidout,'nTransects',dimidt)
   rc = nf90_inquire_dimension(ncidout,dimidr,len=ntransects)
      
   ! dimensions from infile
   rc = nf90_inq_dimid(ncidin,'nCells',dimid)
   rc = nf90_inquire_dimension(ncidin,dimid,len=ncells)
   rc = nf90_inq_dimid(ncidin,'nEdges',dimid)
   rc = nf90_inquire_dimension(ncidin,dimid,len=nEdges)
   rc = nf90_inq_dimid(ncidin,'nVertices',dimid)
   rc = nf90_inquire_dimension(ncidin,dimid,len=nVertices)
   

   rc = nf90_redef(ncidout)
   rc = nf90_def_dim(ncidout,'nVertices',nvertices,dimidv)
   rc = nf90_def_dim(ncidout,'nEdges',nedges,dimide)

   ! read mesh variables from infile   
   allocate(vertexid(nvertices))
   allocate(edgeid(nedges))
   allocate(lone(nedges))
   allocate(late(nedges))
   allocate(angle(nedges))
   allocate(lonc(ncells))
   allocate(latc(ncells))
   allocate(lonv(nvertices))
   allocate(latv(nvertices))
   allocate(xv(nvertices))
   allocate(yv(nvertices))
   allocate(zv(nvertices))
   allocate(xc(ncells))
   allocate(yc(ncells))
   allocate(zc(ncells))
   
   rc = nf90_inq_varid(ncidin,'indexToEdgeID',varid)
   rc = nf90_get_var(ncidin,varid,edgeid)
   rc = nf90_inq_varid(ncidin,'indexToVertexID',varid)
   rc = nf90_get_var(ncidin,varid,vertexid)
   rc = nf90_inq_varid(ncidin,'lonEdge',varid)
   rc = nf90_get_var(ncidin,varid,lone)
   rc = nf90_inq_varid(ncidin,'latEdge',varid)
   rc = nf90_get_var(ncidin,varid,late)
   rc = nf90_inq_varid(ncidin,'angleEdge',varid)
   rc = nf90_get_var(ncidin,varid,angle)
   rc = nf90_inq_varid(ncidin,'lonVertex',varid)
   rc = nf90_get_var(ncidin,varid,lonv)
   rc = nf90_inq_varid(ncidin,'latVertex',varid)
   rc = nf90_get_var(ncidin,varid,latv)
   rc = nf90_inq_varid(ncidin,'lonCell',varid)
   rc = nf90_get_var(ncidin,varid,lonc)
   rc = nf90_inq_varid(ncidin,'latCell',varid)
   rc = nf90_get_var(ncidin,varid,latc)
   rc = nf90_inq_varid(ncidin,'xVertex',varid)
   rc = nf90_get_var(ncidin,varid,xv)
   rc = nf90_inq_varid(ncidin,'yVertex',varid)
   rc = nf90_get_var(ncidin,varid,yv)
   rc = nf90_inq_varid(ncidin,'zVertex',varid)
   rc = nf90_get_var(ncidin,varid,zv)
   rc = nf90_inq_varid(ncidin,'xCell',varid)
   rc = nf90_get_var(ncidin,varid,xc)
   rc = nf90_inq_varid(ncidin,'yCell',varid)
   rc = nf90_get_var(ncidin,varid,yc)
   rc = nf90_inq_varid(ncidin,'zCell',varid)
   rc = nf90_get_var(ncidin,varid,zc)
   
   ! convert radians to degrees
   r2d = 180./ (4.*atan(1.))
   lone = lone * r2d
   late = late * r2d
   lonv = lonv * r2d
   latv = latv * r2d
   lonc = lonc * r2d
   latc = latc * r2d
   
   where (lonc < 0.0) lonc = lonc + 360
   where (lonv < 0.0) lonv = lonv + 360
   where (lone < 0.0) lone = lone + 360
   
   ! normalize cartesian locations to the unit sphere
   do n = 1,nvertices
      cmag = sqrt(xv(n)**2+yv(n)**2+zv(n)**2)
      xv(n) = xv(n) / cmag
      yv(n) = yv(n) / cmag
      zv(n) = zv(n) / cmag
   enddo
   do n = 1,ncells
      cmag = sqrt(xc(n)**2+yc(n)**2+zc(n)**2)
      xc(n) = xc(n) / cmag
      yc(n) = yc(n) / cmag
      zc(n) = zc(n) / cmag
   enddo
   
   ! get mesh variables
   allocate(edgesOnCell(6,ncells))
   rc = nf90_inq_varid(ncidin,'edgesOnCell',varid)
   rc = nf90_get_var(ncidin,varid,edgesOnCell)
   allocate(cellsOnCell(6,ncells))
   rc = nf90_inq_varid(ncidin,'cellsOnCell',varid)
   rc = nf90_get_var(ncidin,varid,cellsonCell)
   allocate(nEdgesOnCell(ncells))
   rc = nf90_inq_varid(ncidin,'nEdgesOnCell',varid)
   rc = nf90_get_var(ncidin,varid,nEdgesOnCell)
   allocate(cellsOnEdge(2,nedges))
   rc = nf90_inq_varid(ncidin,'cellsOnEdge',varid)
   rc = nf90_get_var(ncidin,varid,cellsOnEdge)
   allocate(cellsOnVertex(3,nvertices))
   rc = nf90_inq_varid(ncidin,'cellsOnVertex',varid)
   rc = nf90_get_var(ncidin,varid,cellsOnVertex)
   allocate(edgesOnVertex(3,nvertices))
   rc = nf90_inq_varid(ncidin,'edgesOnVertex',varid)
   rc = nf90_get_var(ncidin,varid,edgesOnVertex)
   allocate(verticesOnEdge(2,nedges))
   rc = nf90_inq_varid(ncidin,'verticesOnEdge',varid)
   rc = nf90_get_var(ncidin,varid,verticesOnEdge)
   
   ! define more transect variables
   rc = nf90_def_var(ncidout,'transectEdgeMasks',nf90_int,(/dimidt,dimide/),varid)
   rc = nf90_def_var(ncidout,'transectEdgeMaskSigns',nf90_int,(/dimidt,dimide/),varid)
   rc = nf90_def_var(ncidout,'transectVertexMasks',nf90_int,(/dimidt,dimidv/),varid)
   
   rc = nf90_enddef(ncidout)

   !define the transects start and end points

   allocate(lon1(ntransects))
   allocate(lonw(ntransects))
   allocate(lat1(ntransects))
   allocate(lonee(ntransects))
   allocate(translat(ntransects))
   
   ! set starting points, longitude limits, and the nominal transect latitude
   !   Keep starting point sufficiently away from land and from 0 degrees longitude
   ! Atlantic_MOC",
   lon1(1) = 340.
   lat1(1) = -33.
   translat(1) = -34.
   lonee(1) = 18.5
   lonw(1) = 307.
   ! AtlanticMed_MOC",
   lon1(2) = 340.
   lat1(2) = -33.
   translat(2) = -34.
   lonee(2) = 18.5
   lonw(2) = 307.
   ! IndoPacific_MOC",
   lon1(3) = 90.
   lat1(3) = -33.
   translat(3) = -34.0
   lonee(3) = 288.5
   lonw(3) = 25.
   ! Pacific_MOC",
   lon1(4) = 200.
   lat1(4) = -5.
   translat(4) = -6.
   lonee(4) = 279.5
   lonw(4) = 104.75
   ! Indian_MOC" ;
   lon1(5) = 60.
   lat1(5) = -5.
   translat(5) = -6.
   lonee(5) = 104.75
   lonw(5) = 43.
   
   
   ! the globalids will be overallocated. the actual values will be concatenated
   !   within the declared vector.
   allocate(transectEdgeMasks(nedges))
   allocate(transectEdgeMaskSigns(nedges))
   allocate(transectVertexMasks(nvertices))
   allocate(regionCellMasks(ncells))
   allocate(transectVertexGlobalIDs(idlim,ntransects))
   allocate(transectEdgeGlobalIDs(idlim,ntransects))
   allocate(ecount(ntransects))
   allocate(vcount(ntransects))
   !
   maxedges = 0
   maxvertices = 0
   vertex1 = 1
   edge1 = 1
   do n = 1,5
   
      print *,'transect ',n
      rc = nf90_inq_varid(ncidout,'regionCellMasks',varid)
      rc = nf90_get_var(ncidout,varid,regionCellMasks,(/n,1/),(/1,ncells/))
      
      call compute_transect_region(lon1(n),lat1(n),lonee(n),lonw(n),   &
                                   translat(n),                        &
                                   ecount(n), vcount(n),               &
                                   transectEdgeMasks,                  &
                                   transectEdgeMaskSigns,              &
                                   transectVertexMasks,                &
                                   regionCellMasks,                    &
                                   transectVertexGlobalIDs(:,n),  &
                                   transectEdgeGlobalIDs(:,n)  )
      
      
      maxedges = max(maxedges, ecount(n) )
      edge1 = edge1 + ecount(n)
      maxvertices = max(maxvertices, vcount(n) )
      vertex1 = vertex1 + vcount(n)
      rc = nf90_inq_varid(ncidout,'transectEdgeMasks',varid)
      rc = nf90_put_var(ncidout,varid,transectEdgeMasks,(/n,1/),(/1,nedges/))
      rc = nf90_inq_varid(ncidout,'transectEdgeMaskSigns',varid)
      rc = nf90_put_var(ncidout,varid,transectEdgeMaskSigns,(/n,1/),(/1,nedges/))
      rc = nf90_inq_varid(ncidout,'transectVertexMasks',varid)
      rc = nf90_put_var(ncidout,varid,transectVertexMasks,(/n,1/),(/1,nvertices/))
   enddo
   
   rc = nf90_redef(ncidout)
   
   print *,'ecount',ecount
   print *,'vcount',vcount
   ! write the transect file global ids
   rc = nf90_def_dim(ncidout,'maxVerticesInTransect',maxvertices,dimid)
   rc = nf90_def_var(ncidout,'transectVertexGlobalIDs',nf90_int,(/dimid,dimidt/),varid)
   rc = nf90_def_dim(ncidout,'maxEdgesInTransect',maxedges,dimid)
   rc = nf90_def_var(ncidout,'transectEdgeGlobalIDs',nf90_int,(/dimid,dimidt/),varid)
      
   rc = nf90_enddef(ncidout)
   allocate(temvertex(maxvertices))
   allocate(temedge(maxedges))
   vertex1 = 1
   edge1 = 1
   do n = 1,ntransects
      temvertex = 0
      temvertex(1:vcount(n)) = transectVertexGlobalIDs(1:vcount(n),n)
      rc = nf90_inq_varid(ncidout,'transectVertexGlobalIDs',varid)
      rc = nf90_put_var(ncidout,varid,temvertex,(/1,n/),(/maxvertices,1/))

      temedge = 0
      temedge(1:ecount(n)) = transectEdgeGlobalIDs(1:ecount(n),n)
      rc = nf90_inq_varid(ncidout,'transectEdgeGlobalIDs',varid)
      rc = nf90_put_var(ncidout,varid,temedge,(/1,n/),(/maxedges,1/))

      edge1 = edge1 + ecount(n)
      vertex1 = vertex1 + vcount(n)
   enddo

   rc = nf90_close(ncidout)

contains

   subroutine compute_transect_region(lon1, lat1, lone, lonw,          &
                                   translat,                           &
                                   ecount, vcount,                     &
                                   transectEdgeMasks,                  &
                                   transectEdgeMaskSigns,              &
                                   transectVertexMasks,                &
                                   regionCellMasks,                    &
                                   transectVertexGlobalIDs,            &
                                   transectEdgeGlobalIDs  )
                                   
   ! This is the main routine for defining a transect and a region
   
   real*8, intent(in) :: lon1, lat1   ! starting point for transect
   real*8, intent(in) :: lone, lonw   ! eastern and western limits to transect
   real*8, intent(in) :: translat     ! nominal latitude of transect
   integer, intent(out) :: ecount, vcount ! count of number of edges (vertices) in the transect
   integer, intent(out), dimension(nedges) :: transectEdgeMasks
   integer, intent(out), dimension(nedges) :: transectEdgeMaskSigns
   integer, intent(out), dimension(nvertices) :: transectVertexMasks
   integer, intent(in), dimension(ncells) :: regionCellMasks
   integer, intent(out), dimension(:) :: transectVertexGlobalIDs
   integer, intent(out), dimension(:) :: transectEdgeGlobalIDs
   
   real*8 dotmax, dotloc, x1, y1, z1, dellat, dellatmin
   integer n, n1, currcell, nextcell, lastcell, m, k , &
           edge1, nextedge, vertexeast, vertexwest,  &
           currvert, lastvert, nextvert, curredge
   logical headsouth, headeast, headwest, restrictlat
   
   integer, dimension(idlim,2) :: vertids, edgeids
   integer, dimension(2) :: vcounttem, ecounttem
   
   ! initialization
   ecount = 0
   vcount = 0
   ecounttem = 0
   vcounttem = 1
   transectEdgeMasks = 0
   transectEdgeMaskSigns = 0
   transectVertexMasks = 0

   ! the starting point was chosen in open ocean where one can intersect the transect
   !   by going south  
   
   ! find the nearest cell with mask = 1
   x1 = cos((lon1)/r2d) * cos(lat1/r2d)
   y1 = sin((lon1)/r2d) * cos(lat1/r2d)
   z1 =                   sin(lat1/r2d)
   dotmax = -1.
   do n = 1,ncells
      dotloc = x1*xc(n) + y1*yc(n) + z1*zc(n)
      if(dotloc > dotmax) then
         dotmax = dotloc
         n1 = n
      endif
   enddo
   
   ! find the southernmost cell neighbor
   headsouth = .true.
   currcell = n1
   lastcell = n1
   do while (headsouth)
   
      ! if neighbor cell mask = 0 then the edge between these cells is our starting point
      !  we will compute a segment heading west, and a segment heading east, then
      !  combine them
      dellatmin = 999.
      nextcell = 0
      do m = 1,nedgesOnCell(currcell)
         k = cellsOnCell(m,currcell)
         if(k.ne.0 .and. k.ne.lastcell) then
            dellat = latc(k) - latc(currcell)
            if( dellat < dellatmin ) then
               dellatmin = dellat
               nextcell = k
               edge1 = edgesOnCell(m,currcell)
            endif
         endif
      enddo
      
      if (nextcell == 0) stop 'heading south failure'
      
      if ( regionCellMasks(nextcell) == 1) then
         lastcell = currcell
         currcell = nextcell
      else
         headsouth = .false.
      endif
   
   enddo
   
   ! test that we are neat translat
   if( abs(translat-latc(currcell)) > 0.5) stop  'did not reach translat'
   
   ! we identify the eastern and western vertices
   if(lonv(verticesOnEdge(1,edge1))>lonv(verticesOnEdge(2,edge1))) then
      vertexeast = verticesOnEdge(1,edge1)
      vertexwest = verticesOnEdge(2,edge1)
   else
      vertexeast = verticesOnEdge(2,edge1)
      vertexwest = verticesOnEdge(1,edge1)
   endif
   
   vertids(1,1) = vertexeast
   vertids(1,2) = vertexwest
   
   ! proceed east from the eastern vertex and compute the eastern segment
   !  along the regionCellMasks 0/1 boundary
   
   ! we may head north. If we are more than three degrees from the eastern
   !  longitude limit we may depart from translat as far as necessary.
   !  Else we will end the segment two degrees from translat.
   headeast = .true.
   currvert = vertexeast
   lastvert = vertexwest
   restrictlat = .false.
   curredge = edge1
   do while( headeast )
      do m = 1,3
         k = edgesOnVertex(m,currvert)
         if( k .ne. curredge .and. k .ne. 0) then
            if(cellsonedge(1,k) == 0) then
               if( regionCellMasks(cellsonedge(2,k)) == 1)     &
                     nextedge = k
            elseif(cellsonedge(2,k) == 0) then
               if( regionCellMasks(cellsonedge(1,k)) == 1)     &
                     nextedge = k
            else 
               if(regionCellMasks(cellsonedge(1,k))+regionCellMasks(cellsonedge(2,k)) == 1)   &
                     nextedge = k
            endif
         endif
      enddo
      do m = 1,2
         k = verticesOnEdge(m,nextedge)
         if(k .ne. currvert) then
            nextvert = k
         endif
      enddo
      curredge = nextedge
      lastvert = currvert
      currvert = nextvert
      ecounttem(1) = ecounttem(1) + 1
      vcounttem(1) = vcounttem(1) + 1
      if( vcounttem(1) > idlim) stop 'count exceeds limit east'
      vertids(vcounttem(1),1) = currvert
      edgeids(ecounttem(1),1) = nextedge
      if(abs(lone-lonv(currvert)) < 3.)restrictlat = .true.
      if(restrictlat) then
         if(abs(translat-latv(currvert)) > 2.) headeast = .false.
      endif
   enddo
   
   ! Now we start from the original western vertex and proceed west
   headwest = .true.
   currvert = vertexwest
   lastvert = vertexeast
   restrictlat = .false.
   curredge = edge1
   do while( headwest )
      do m = 1,3
         k = edgesOnVertex(m,currvert)
         if( k .ne. curredge .and. k .ne. 0) then
            if(cellsonedge(1,k) == 0) then
               if( regionCellMasks(cellsonedge(2,k)) == 1)     &
                     nextedge = k
            elseif(cellsonedge(2,k) == 0) then
               if( regionCellMasks(cellsonedge(1,k)) == 1)     &
                     nextedge = k
            else 
               if(regionCellMasks(cellsonedge(1,k))+regionCellMasks(cellsonedge(2,k)) == 1)   &
                     nextedge = k
            endif
         endif
      enddo
      do m = 1,2
         k = verticesOnEdge(m,nextedge)
         if(k .ne. currvert) then
            nextvert = k
         endif
      enddo
      curredge = nextedge
      lastvert = currvert
      currvert = nextvert
      ecounttem(2) = ecounttem(2) + 1
      vcounttem(2) = vcounttem(2) + 1
      if( vcounttem(2) > idlim) stop 'count exceeds limit west'
      vertids(vcounttem(2),2) = currvert
      edgeids(ecounttem(2),2) = nextedge
      if(abs(lonw-lonv(currvert)) < 3.)restrictlat = .true.
      if(restrictlat) then
         if(abs(translat-latv(currvert)) > 2.) headwest = .false.
      endif
   enddo
   
   print*,'ecounttem',ecounttem
   print*,'vcounttem',vcounttem
   ! stitch the two segments together
   vcount = 1
   transectVertexGlobalIDs(vcount) = vertids(vcounttem(2),2)
   do vcount = 2,vcounttem(2)
      transectVertexGlobalIDs(vcount) = vertids(vcounttem(2)+1-vcount,2)
      transectedgeglobalIDs(vcount-1) = edgeids(vcounttem(2)+1-vcount,2)
   enddo
   transectedgeglobalIDs(vcounttem(2)) = edge1

   vcount = vcounttem(2) + 1
   transectVertexGlobalIDs(vcount) = vertexeast
   do n = 2,vcounttem(1)
      transectedgeglobalIDs(vcount) = edgeids(n-1,1)
      vcount = vcount + 1
      transectVertexGlobalIDs(vcount) = vertids(n,1)
   enddo  
   
   ecount = vcount - 1
   do n = 1,vcount
       transectVertexMasks(transectVertexGlobalIDs(n)) = 1
   enddo
   
   do n = 1,ecount
       k = transectEdgeGlobalIDs(n)
       transectEdgeMasks(k) = 1
       transectEdgeMaskSigns(k) = 1
       if ( angle(k) < 0. ) transectEdgeMaskSigns(k) = -1
       if (minval(cellsOnEdge(:,k))==0)transectEdgeMaskSigns(k) = -1

   enddo
   
   end subroutine compute_transect_region

end program mpaso_transects
