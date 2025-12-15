program mpas2esmf
! This program takes the grid information from an MPAS ocean domain mesh file
!   and converts it to an ESMF compliant grid file that NUOPC can read.
!

use netcdf

implicit none

! variable declarations

   ! netcdf variables
   integer ncidin, ncidout, rc
   integer dimid, dimid_nodecount, dimid_elementcount, &
           dimid_maxnodepelement, dimid_coorddim
   integer varid
   character(len=128) :: filein, fileout

   ! dimensions
   integer nodeCount    ! nVertices on input
   integer elementCount ! nCells on input
   integer, parameter :: maxNodePElement = 6 ! QU
   !integer, parameter :: maxNodePElement = 7 ! EC
   integer, parameter :: coordDim = 2

   ! input arrays
   real*8, allocatable, dimension(:) :: latCell, lonCell,      &
                                        latVertex, lonVertex,  &
                                        areaCell
   integer, allocatable, dimension(:,:) :: verticesOnCell
   integer, allocatable, dimension(:) :: nEdgesOnCell
   ! input attribute
   real*8 radius

   ! output arrays
   real*8, allocatable, dimension(:,:) :: nodeCoords, centerCoords
   real*8, allocatable, dimension(:) :: elementArea
   integer, allocatable, dimension(:) :: elementMask
   integer, allocatable, dimension(:,:) :: elementConn
   integer*1, allocatable, dimension(:) :: numElementConn

   ! work variables
   integer i, n, ncon
   real*8 r2d, areafac

   !test unique nodes
   logical :: nodeunique
   integer m, ne, nlen, k

! begin executable code

   ! assign filenames, open the input and create the output files
   print*,'enter the ocean domain mesh file name'
   read(5,'(a128)')filein
   print*,'enter the ocean domain ESMF mesh file name'
   read(5,'(a128)')fileout

   rc = nf90_open(trim(filein),nf90_nowrite,ncidin)
   if(rc .ne. 0) stop 'troubles with the input ocean domain mesh file - stopping'
   rc = nf90_create(trim(fileout),ior(nf90_clobber,nf90_netcdf4),ncidout)
   if(rc .ne. 0) stop 'troubles creating the ocean domain ESMF mesh file - stopping'

   ! get dimensions and allocate arrays
   rc = nf90_inq_dimid(ncidin,'nCells', dimid)
      rc = nf90_inquire_dimension(ncidin, dimid, len=elementCount)
   rc = nf90_inq_dimid(ncidin,'nVertices', dimid)
      rc = nf90_inquire_dimension(ncidin, dimid, len=nodeCount)

   allocate(latCell(elementCount))
   allocate(lonCell(elementCount))
   allocate(areaCell(elementCount))
   allocate(nEdgesOnCell(elementCount))
   allocate(latVertex(nodeCount))
   allocate(lonVertex(nodeCount))
   allocate(verticesOnCell(maxNodePElement, elementCount))

   allocate(nodeCoords(coordDim, nodeCount))
   allocate(centerCoords(coordDim, elementCount))
   allocate(elementArea(elementCount))
   allocate(elementMask(elementCount))
   allocate(numElementConn(elementCount))
   allocate(elementConn(maxNodePElement, elementCount))

   ! define the output file
   ! dimensions
   rc = nf90_def_dim(ncidout, 'nodeCount', nodeCount, dimid_nodeCount)
   print*,'def_dim ',rc
   rc = nf90_def_dim(ncidout, 'elementCount', elementCount, dimid_elementCount)
   print*,'def_dim ',rc
   rc = nf90_def_dim(ncidout, 'maxNodePElement', maxNodePElement, dimid_maxNodePElement)
   print*,'def_dim ',rc
   rc = nf90_def_dim(ncidout, 'coordDim', coordDim, dimid_coordDim)
   print*,'def_dim ',rc

   ! variables
   rc = nf90_def_var(ncidout, 'nodeCoords', nf90_double, (/dimid_coordDim,dimid_nodeCount/), varid )
   print*,'def_var ',rc
      rc = nf90_put_att(ncidout, varid, 'units', 'degrees')
   rc = nf90_def_var(ncidout, 'elementConn', nf90_int, (/dimid_maxNodePElement,dimid_elementCount/), varid )
      rc = nf90_put_att(ncidout, varid, 'long_name', 'Node Indices that define the element connectivity')
      rc = nf90_put_att(ncidout, varid, '_FillValue', -1)
   rc = nf90_def_var(ncidout, 'numElementConn', nf90_byte, (/dimid_elementCount/), varid )
   print*,'def_var ',rc
      rc = nf90_put_att(ncidout, varid, 'long_name', 'Number of nodes per element')
   rc = nf90_def_var(ncidout, 'centerCoords', nf90_double, (/dimid_coordDim,dimid_elementCount/), varid )
   print*,'def_var ',rc
      rc = nf90_put_att(ncidout, varid, 'units', 'degrees')
   rc = nf90_def_var(ncidout, 'elementArea', nf90_double, (/dimid_elementCount/), varid )
   print*,'def_var ',rc
      rc = nf90_put_att(ncidout, varid, 'units', 'radians^2')
      rc = nf90_put_att(ncidout, varid, 'long_name', 'area weights')
   rc = nf90_def_var(ncidout, 'elementMask', nf90_int, (/dimid_elementCount/), varid )
   print*,'def_var ',rc

   ! global attributes
   rc = nf90_put_att(ncidout, nf90_global, 'gridType', 'unstructured')
   print*,'put_att ',rc
   rc = nf90_put_att(ncidout, nf90_global, 'version', '1.0')
   print*,'put_att ',rc
   rc = nf90_put_att(ncidout, nf90_global, 'inputFile', trim(filein) )
   print*,'put_att ',rc
   rc = nf90_put_att(ncidout, nf90_global, 'timeGenerated', '08/04/21')
   print*,'put_att ',rc
   rc = nf90_put_att(ncidout, nf90_global, 'history', 'created by Dazlich mpas2esmf program')
   print*,'put_att ',rc

   rc = nf90_enddef(ncidout)

   ! read
   rc = nf90_inq_varid(ncidin,'latCell',varid)
      rc = nf90_get_var(ncidin, varid, latCell)
print*,'getvar latcell ',rc
   rc = nf90_inq_varid(ncidin,'lonCell',varid)
      rc = nf90_get_var(ncidin, varid, lonCell)
print*,'getvar loncell ',rc
   rc = nf90_inq_varid(ncidin,'latVertex',varid)
      rc = nf90_get_var(ncidin, varid, latVertex)
print*,'getvar latvertex ',rc
   rc = nf90_inq_varid(ncidin,'lonVertex',varid)
      rc = nf90_get_var(ncidin, varid, lonVertex)
print*,'getvar lonvertex ',rc
   rc = nf90_inq_varid(ncidin,'verticesOnCell',varid)
      do n = 1,elementCount,10000
      ne = min(elementCount,n+9999)
      nlen = min(10000,elementCount+1-n)
      rc = nf90_get_var(ncidin, varid, verticesOnCell(:,n:ne),(/1,n/),(/maxNodePElement,nlen/))
print*,'getvar verticesoncell ',rc
      enddo
   rc = nf90_inq_varid(ncidin,'areaCell',varid)
      rc = nf90_get_var(ncidin, varid, areaCell)
print*,'getvar areacell ',rc
   rc = nf90_inq_varid(ncidin,'nEdgesOnCell',varid)
      rc = nf90_get_var(ncidin, varid, nEdgesOnCell)
print*,'getvar nedgesoncell ',rc
   rc = nf90_get_att(ncidin, nf90_global, 'sphere_radius',radius)

   ! copy, convert units and compute numElementConn
   r2d = 180._8 / (4._8 * atan(1._8))
   areafac = 1. / (radius**2)
print*,'radius ',radius
   do i = 1,nodeCount
      nodeCoords(1,i) = lonVertex(i) * r2d
      nodeCoords(2,i) = latVertex(i) * r2d
   enddo

   do i = 1,elementCount
      centerCoords(1,i) = lonCell(i) * r2d
      centerCoords(2,i) = latCell(i) * r2d
      elementArea(i) = areaCell(i) * areafac
      elementMask(i) = 1
      ncon = 0
      
      !! find n pole
      !if(centerCoords(2,i) > 89.9) then
      !  print*,'n. pole at cell ',i
      !  do m = 1,6
      !    print*,'node ',verticesoncell(m,i),nodecoords(:,verticesoncell(m,i))
      !  enddo
      !endif
      !nodeunique = .true.
      do n = 1,nEdgesOnCell(i)
         !if( verticesOnCell(n,i) <= nodeCount ) then
            ncon = ncon + 1
            elementConn(ncon,i) = verticesOnCell(n,i)
         !endif
      !   !uniqueness test
      !   do m = 1,n-1
      !      if(verticesOnCell(n,i) == verticesOnCell(m,i)) nodeunique = .false.
      !   enddo
      enddo
      !if(.not.nodeunique) print*,'nodes not unique for cell ',i,'nodes ',verticesOnCell(:,i)
      do n = nEdgesOnCell(i)+1,maxNodePElement
         elementConn(n,i) = -1
      enddo
      numElementConn(i) = ncon
   enddo

   ! write
   rc = nf90_inq_varid(ncidout,'nodeCoords',varid)
      rc = nf90_put_var(ncidout,varid, nodeCoords)
   rc = nf90_inq_varid(ncidout,'elementConn',varid)
      rc = nf90_put_var(ncidout,varid, elementConn)
   rc = nf90_inq_varid(ncidout,'numElementConn',varid)
      rc = nf90_put_var(ncidout,varid, numElementConn)
   rc = nf90_inq_varid(ncidout,'centerCoords',varid)
      rc = nf90_put_var(ncidout,varid, centerCoords)
   rc = nf90_inq_varid(ncidout,'elementArea',varid)
      rc = nf90_put_var(ncidout,varid, elementArea)
   rc = nf90_inq_varid(ncidout,'elementMask',varid)
      rc = nf90_put_var(ncidout,varid, elementMask)
 
   ! close the files
   rc = nf90_close(ncidin)
   rc = nf90_close(ncidout)

   !check that the global area of the ocean domain is reasonable
   print*,'4*pi, sum(elementArea)',16.*atan(1.),sum(elementArea)

end program mpas2esmf
