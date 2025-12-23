program modify_ic
   use netcdf
implicit none
   character(len=128)::modfile,rstfile
   integer :: ncells, nedges, rc, ncidin, ncidout, varid, dimid, n, idcells, idtime
   real*8, dimension(:), allocatable :: var1a, var1b, area, ssh
   real*8, dimension(:,:), allocatable :: var2
   real*8 :: garea, gssh

   print*,'enter the new ocean restart file name'       
   read(5,'(a128)')rstfile
   print*,'enter the oceanIC mesh file name'       
   read(5,'(a128)')modfile
   !rstfile = 'omipv2.5_015.mpaso.rst.0001-01-01_00000.nc'
   !modfile = 'oceanIC.015km.250908-2000-01-01.nc'

   rc = nf90_open(trim(rstfile),nf90_nowrite,ncidin)
   rc = nf90_open(trim(modfile),nf90_write,ncidout)

   rc = nf90_inq_dimid(ncidout,'nCells',dimid)
   idcells = dimid
   rc = nf90_inquire_dimension(ncidout,dimid,len=ncells)
print*,'ncells ',ncells
   allocate(var1a(ncells))
   allocate(var1b(ncells))
   allocate(area (ncells))
   allocate(ssh (ncells))

   allocate(var2(100,ncells))

   rc = nf90_inq_varid(ncidin,'areaCell',varid)
print*,'thick inqvarid ',rc
   rc = nf90_get_var(ncidin,varid,area)
print*,'thick getvar ',rc
   garea = sum(area)

   rc = nf90_inq_varid(ncidin,'layerThickness',varid)
print*,'thick inqvarid ',rc
   rc = nf90_get_var(ncidin,varid,var2)
print*,'thick getvar ',rc

   rc = nf90_inq_varid(ncidin,'bottomDepth',varid)
print*,'depth inqvarid ',rc
   rc = nf90_get_var(ncidin,varid,var1a)
print*,'depth getvar ',rc
print*,'ncells ',ncells
print*,'size(var1a)',size(var1a)
print*,'size(var1b)',size(var1b)
print*,'size(var2)',size(var2)
   do n = 1,ncells
      ssh(n) = sum(var2(:,n)) - var1a(n)
   enddo

   rc = nf90_inq_varid(ncidin,'atmosphericPressure',varid)
print*,'atmosP inqvarid ',rc
   rc = nf90_get_var(ncidin,varid,var1a)
print*,'atmosP getvar ',rc
   rc = nf90_inq_varid(ncidin,'seaIcePressure',varid)
print*,'seaiceP inqvarid ',rc
   rc = nf90_get_var(ncidin,varid,var1b)
print*,'seaiceP getvar ',rc
   ssh = ssh + var1b / (1027. * 9.806)
   var1a = var1a + var1b
   rc = nf90_inq_varid(ncidout,'seaSurfacePressure',varid)
!print*,'ssP inqvarid ',rc
   rc = nf90_put_var(ncidout,varid,var1a)
print*,'ssP putvar ',rc

   gssh = sum(ssh*area) / garea
   print*,'gssh',gssh
   ssh = ssh - gssh

   ssh = ssh - var1b / (1027. * 9.806)
   rc = nf90_inq_varid(ncidout,'ssh',varid)
   rc = nf90_put_var(ncidout,varid,ssh)
print*,'ssh putvar ',rc

   rc = nf90_close(ncidout)

end program modify_ic
