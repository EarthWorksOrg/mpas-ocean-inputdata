program create_graph
! This program takes the grid information from an MPAS ocean mesh file    
!   and converts it to an ocean domain graph file that metis can perform decompositions on.
!

use netcdf

implicit none

! variable declarations

   ! netcdf variables
   integer ncid, rc
   integer dimid
   integer varid
   character(len=128) :: filein, fileout

   ! dimensions
   integer nedges
   integer ncells

   ! input arrays
   integer, allocatable, dimension(:) :: nedgesoncell
   integer, allocatable, dimension(:,:) :: cellsOnCell
   integer :: cells(6)

   ! work variables
   integer i, n, ne, nlen

! begin executable code

   ! assign filenames, open the input and create the output files
   print*,'enter the ocean domain mesh file name'
   read(5,'(a128)')filein
   print*,'enter the ocean domain graph file name'
   read(5,'(a128)')fileout


   rc = nf90_open(trim(filein),nf90_nowrite,ncid)
   if(rc .ne. 0) stop 'trouble with the input ocean mesh file - stopping'

   ! get dimensions and allocate arrays
   rc = nf90_inq_dimid(ncid,'nCells', dimid)
      rc = nf90_inquire_dimension(ncid, dimid, len=nCells)
   rc = nf90_inq_dimid(ncid,'nEdges', dimid)
      rc = nf90_inquire_dimension(ncid, dimid, len=nEdges)

   allocate(nEdgesOnCell(ncells))
   allocate(cellsOnCell(6,ncells))

   ! read
   rc = nf90_inq_varid(ncid,'nEdgesOnCell',varid)
     do n=1,ncells,10000
      ne = min(ncells,n+9999)
      nlen = min(10000,ncells+1-n)
      rc = nf90_get_var(ncid, varid, nedgesoncell(n:ne),(/n/),(/nlen/))
     enddo
   rc = nf90_inq_varid(ncid,'cellsOnCell',varid)
     do n=1,ncells,10000
      ne = min(ncells,n+9999)
      nlen = min(10000,ncells+1-n)
      rc = nf90_get_var(ncid, varid, cellsOnCell(:,n:ne),(/1,n/),(/6,nlen/))
     enddo

   open(unit=7,file=trim(fileout),form='formatted')
   write(7,*) ncells, nedges

   ! copy, convert units and compute numElementConn
   do i = 1,ncells
      nedges = 0
      do n = 1,nedgesoncell(i)
         if(cellsoncell(n,i) > 0) then
            nedges = nedges + 1
            cells(nedges) = cellsoncell(n,i)
         endif
      enddo
      write(7,*)cells(1:nedges)
!print*,cellsoncell(:,i)
   enddo

 
   ! close the files
   rc = nf90_close(ncid)


end program create_graph
