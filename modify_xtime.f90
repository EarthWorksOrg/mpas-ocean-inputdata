program modify_ic
   use netcdf
implicit none
   character(len=128)::ocnfile,icefile
   integer :: ncells, nedges, rc, ncidin, ncidout, varid, dimid, n, idcells, idtime
   character(len=64) :: xtime

!   ocnfile = 'omip_015.lev2/run/omip_015.lev2.mpaso.rst.0001-01-01_00000.nc'  
   print*,'enter the new seaice restart file name'       
   read(5,'(a128)')icefile
   !icefile = 'omipv2.5_015.mpassi.rst.0001-01-01_00000.nc'

!   rc = nf90_open(trim(ocnfile),nf90_write,ncidin)
!print*,'ocn open'
!   rc = nf90_inq_varid(ncidin,'xtime',varid)
!print*,'ocn inqvarid ',rc
!   rc = nf90_get_var(ncidin,varid,xtime)
!print*,'ocn getvar ',rc
!   xtime(1:10)='0001-01-01'
!   rc = nf90_put_var(ncidin,varid,xtime)
!print*,'ocn putvar ',rc
!   rc = nf90_close(ncidin)
!print*,'ocn close'

   rc = nf90_open(trim(icefile),nf90_write,ncidin)
print*,'ice open'
   rc = nf90_inq_varid(ncidin,'xtime',varid)
print*,'ice inqvarid ',rc
   rc = nf90_get_var(ncidin,varid,xtime)
print*,'ice getvar ',rc
   xtime(1:10)='0001-01-01'
   rc = nf90_put_var(ncidin,varid,xtime)
print*,'ice putvar ',rc
   rc = nf90_close(ncidin)
print*,'ice close'

end program modify_ic
