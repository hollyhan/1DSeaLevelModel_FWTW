!=================================================================================PHYSICSAL & MATHEMATICAL CONSTANTS===!
module constants_mod
!______________________________________________________________________________________________________________________!
   real, parameter :: pi = 3.1415926535898      ! Pi
   complex, parameter :: ii=(0.0,1.0)           ! Square root of -1
   real, parameter :: gravConst = 6.67408E-11   ! Gravitational constant (m^3/kg/s^2)
end module constants_mod

!========================================================================================================PLANETS_MOD===!
module planets_mod
!______________________________________________________________________________________________________________________!
   use constants_mod
   
   real :: radius
   real :: mass
   real :: rhoi
   real :: rhow
   real :: gacc
   real :: omega 
   real :: acoef, ccoef
   real :: moiA, moiC
   real :: kf
   
   contains
   
   !=========================================================================================PLANETS_MOD: EARTH_INIT===!
   subroutine earth_init
   !___________________________________________________________________________________________________________________!
      radius = 6.371E6              ! Radius of the Earth (m)
      mass = 5.976E24               ! Mass of the Earth (kg)
      rhoi = 920.0                  ! Density of ice (kg/m^3)
      rhow = 1000.0                 ! Density of fresh water (kg/m^3)
      gacc = 9.80665                ! Acceleration due to gravity at the Earth's surface (m/s^2)
      omega = 7.292e-5              ! Rotation rate of the Earth (rad/s)
      moiA=0.3296145*mass*radius**2 ! Principal moment of inertia of the Earth
      moiC=0.3307007*mass*radius**2 ! Principal moment of inertia of the Earth
      kf = 0.9342+0.008             ! Fluid (Tidal) Love number
      
   end subroutine earth_init
   
   !==========================================================================================PLANETS_MOD: MARS_INIT===!
   subroutine mars_init
   !___________________________________________________________________________________________________________________!
      radius = 3.3899E6             ! Radius of Mars (m)
      mass = 6.4185E23              ! Mass of Mars (kg)
      rhoi = 1220.0                 ! Density of ice mix on Mars (kg/m^3)
      rhow = 1000.0                 ! Density of fresh water (kg/m^3)
      gacc = 3.713                  ! Acceleration due to gravity at Mars's surface (m/s^2)
      omega = 7.08819118E-5         ! Rotation rate of Mars (rad/s)
      moiA=0.363914*mass*radius**2  ! Principal moments of inertia of Mars
      moiC=0.365905*mass*radius**2  ! Principal moments of inertia of Mars
      !----------------------------------------------------------------------------------!
      ! Lithospheric thickness and corresponding kf (based on Zharkov and Gudkova, 2005)
      !    15    |    48    |    58    |    84    |    110   |    164    |    200    [km]
      ! 1.203959 | 1.127384 | 1.110566 | 1.065899 | 1.023186 | 0.9458358 | 0.8986673
      !----------------------------------------------------------------------------------!
      kf = 0.899                    ! Fluid (Tidal) Love number
      
   end subroutine mars_init
   
end module planets_mod


!! =====================================================================================================netCDF I/O mod===! 
module io_mod
!---------------------------------------------------------------------------------------------------------------------------
   use user_specs_mod, only: nglv, ext, fType, outputfolder, gridfolder, grid_lat, grid_lon
   use netCDF
   implicit none

   contains
	
   subroutine check(status)
    
	 integer, intent ( in) :: status
    
     if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped"
     end if
   
   end subroutine check
   
   subroutine read_sl(data_slm, filename, filepath, suffix, fext)
	   character (len = *), intent (in) :: filename, filepath
	   real, dimension(nglv,2*nglv), intent(inout) :: data_slm 
	 !  integer, intent(in), optional :: indx_time, indx_topoloop
	   character (len = *), optional :: suffix, fext
	 !  character(8) :: suffix
	 

	   ! indentify attributes
	!   select case(filename)
	!      case ("tgrid")
!!		     filepath = outputfolder
!	      case ("beta")
!		     filepath = outputfolder
!	      case ("ocean")
!		     filepath = outputfolder 
 !      end select

	   
	 !  suffix = ''
	   

	!   if (present (indx_time)) then
	!	   write(numstr,'(I4)') indx_time
	!	   numstr = trim(adjustl(numstr))
	!	   suffix == suffix//numstr
	 !  endif
	   
	 !  if (present (indx_topoloop)) then 
	!	   write(numstr,'(I4)') indx_topoloop
	!	   numstr = trim(adjustl(numstr))
	!	   suffix == suffix//'_'//numstr
	!   endif
	   
	   if (fType == 'text') then
		   call read_txt(data_slm, filename, filepath, suffix, fext)
	   elseif (fType == 'binary') then 
		   call read_nf90(data_slm, filename, filepath, suffix, fext)
	   endif
	   
   end subroutine read_sl
   
   subroutine write_sl(data_slm, filename, filepath, suffix, fext)
	   character (len = *), intent (in) :: filename, filepath
	   real, dimension(nglv,2*nglv), intent(in) :: data_slm 
	   character (len = *), optional :: suffix, fext
	   
	   ! indentify attributes
	  ! select case(filename)
	   !   case ("tgrid")
		!     filepath = outputfolder
	     ! case ("beta")
		 !    filepath = outputfolder
	     ! case ("ocean")
		 !    filepath = outputfolder 
       !end select

	   if (fType == 'text') then
		   call write_txt(data_slm, filename, filepath, suffix, fext)
	   elseif (fType == 'binary') then 
		   call write_nf90(data_slm, filename, filepath, suffix, fext)
	   endif
	   
   end subroutine write_sl
    
   
   subroutine write_nf90(data_slm, filename, filepath, suffix, fext)

	   character (len = *), intent(in) :: filename, filepath
	   integer :: ncid, varid, lat_varid, lon_varid, lat_dimid, lon_dimid
	   integer, dimension(2) :: dimids
	   real, dimension(nglv,2*nglv), intent(in) :: data_slm !data in the SLM written to the netCDF file
	   character (len = *), optional ::  suffix 
	   character (len = *), optional :: fext
	   real, dimension(nglv)        :: latgrid
	   real, dimension(2*nglv)      :: longrid
	   
	   ! attribute IDs for I/O in netCDF
	   !character (len = *), parameter :: UNITS = "units"
	   !character (len = *), parameter :: TOPO_UNITS = "meters"
	   !character (len = *), parameter :: SL_UNITS = "meters"
	   !character (len = *), parameter :: LAT_UNITS = "degrees"
	   !character (len = *), parameter :: LON_UNITS = "degrees_east"
	   
	   
	   ! Read in lat-lon grid files
	   open(unit = 1, file = gridfolder//grid_lat, form = 'formatted', access = 'sequential', status = 'old')
	   read(1,*) latgrid
	   close(1)

	   open(unit = 1, file = gridfolder//grid_lon, form = 'formatted', access = 'sequential', status = 'old')
	   read(1,*) longrid
	   close(1)
	   
       !write out data
	   
	   !create the file
   	   if (present (suffix)) then
 		  if (present (fext)) then 
 		     call check( nf90_create(filepath//filename//trim(suffix)//fext, nf90_clobber, ncid) ) 
 		  else
 			 call check( nf90_create(filepath//filename//trim(suffix)//ext, nf90_clobber, ncid) ) 
 		  endif
       else
		  if (present (fext)) then 
		     call check( nf90_create(filepath//filename//fext, nf90_clobber, ncid) ) 
		  else
			 call check( nf90_create(filepath//filename//ext, nf90_clobber, ncid) ) 
		  endif
   	   endif
	   
       call check( nf90_def_dim(ncid, 'lon', nglv*2, lon_dimid)  ) ! Define the dimensions of the griddata
       call check( nf90_def_dim(ncid, 'lat', nglv,   lat_dimid)  )

       dimids =  (/ lon_dimid, lat_dimid /)                        ! Define the dimension of the variable

       call check( nf90_def_var(ncid, filename, nf90_double, dimids, varid)) ! Define variable
       call check( nf90_def_var(ncid, 'lat', nf90_double, lat_dimid, lat_varid))
       call check( nf90_def_var(ncid, 'lon', nf90_double, lon_dimid, lon_varid))
       call check( nf90_enddef(ncid)) ! End definition

       call check( nf90_put_var(ncid, varid, reshape(data_slm,[2*nglv,nglv]))) !write data
   	   call check( nf90_put_var(ncid, lon_varid, longrid))
   	   call check( nf90_put_var(ncid, lat_varid, latgrid))
       call check( nf90_close(ncid))
	   
   end subroutine write_nf90
   
   subroutine read_nf90(data_slm, filename, filepath, suffix, fext)
	     
	   character (len = *), intent(in) :: filename,  filepath !file name and path
	   real, dimension(2*nglv,nglv) :: data_temp !temp. variable name in the SLM in which nc data will be stored
	   real, dimension(nglv,2*nglv), intent(out) :: data_slm
	   character (len = *), optional ::  suffix 
	   character (len = *), optional :: fext
	   integer :: ncid, varid 
	   
	   !This does not work - why? 
       !if (.not.present (fext)) then 
      !	    fext == ext
     !  end if
	   
   	   if (present (suffix)) then
 		  if (present (fext)) then 
 		     call check( nf90_open(filepath//filename//trim(suffix)//fext, nf90_nowrite, ncid) ) !open the file
 		  else
 			 call check( nf90_open(filepath//filename//trim(suffix)//ext, nf90_nowrite, ncid) ) !open the file
 		  endif
       else
		  if (present (fext)) then 
		     call check( nf90_open(filepath//filename//fext, nf90_nowrite, ncid) ) !open the file
		  else
			 call check( nf90_open(filepath//filename//ext, nf90_nowrite, ncid) ) !open the file
		  endif
   	   endif
	   
	   call check( nf90_inq_varid(ncid, filename, varid) ) !get varid of the data variable
	   call check( nf90_get_var(ncid, varid, data_temp) ) ! read the data
	   call check( nf90_close(ncid) ) ! close the file
	   data_slm = reshape(data_temp,[nglv,2*nglv])
	   
   end subroutine read_nf90
   

  ! subroutine read_nf90(data_slm, filename, varname_nc, data_slm)!
	     
!	   character (len = *), intent(in) :: filename, varname_nc !file name and variable name in nc file !
!	   real, dimension(2*nglv,nglv) :: data_temp !temp. variable name in the SLM in which nc data will be stored
!	   real, dimension(nglv,2*nglv), intent(out) :: data_slm
!	   integer :: ncid, varid 
	   
!	   call check( nf90_open(filename, nf90_nowrite, ncid) ) !open the file
!	   call check( nf90_inq_varid(ncid, varname_nc, varid) ) !get varid of the data variable
!	   call check( nf90_get_var(ncid, varid, data_temp) ) ! read the data
!	   call check( nf90_close(ncid) ) ! close the file
!	   data_slm = reshape(data_temp,[nglv,2*nglv])
	   
 !  end subroutine read_nf90
  
   
   subroutine check_ncFile(ncvar, slmvar, varname) 
	   
	   real, dimension(nglv,2*nglv), intent(in) :: ncvar, slmvar 
	   character (len = *), intent(in) :: varname 
	   
	   if (sum(sum(ncvar, dim=1))-sum(sum(slmvar,dim=1)) .NE.0) then
		   write(6,*) 'ncFile benchmark test FAILED. Var name is:                              ', varname
	   else if (sum(sum(ncvar, dim=1))-sum(sum(slmvar, dim=1)).eq.0) then
	       write(6,*) 'ncFile benchmark test PASSED. Var name is:                              ', varname
	   end if
			   
   end subroutine check_ncFile
   
   subroutine write_txt(data_slm, filename, filepath, suffix, fext)
   	
	real, dimension(nglv,2*nglv), intent(in) :: data_slm
   	character (len = *), intent(in) :: filename, filepath 
   	character (len = *), optional :: fext, suffix

    if (.not.present (fext)) then 
   	    fext = ext
    end if

	if (present (suffix)) then
		open(unit = 1, file = filepath//filename//trim(suffix)//fext, form = 'formatted', access = 'sequential', &
        & status = 'replace')
        write(1,'(ES16.9E2)') data_slm
        close(1)
	else
		open(unit = 1, file = filepath//filename//fext, form = 'formatted', access = 'sequential', &
        & status = 'replace')
        write(1,'(ES16.9E2)') data_slm
        close(1)
	endif

   end subroutine write_txt
   
   subroutine read_txt(data_slm, filename, filepath, suffix, fext)
   	
	real, dimension(nglv,2*nglv) :: data_slm
   	character (len = *), intent(in) :: filename, filepath 
   	character (len = *), optional :: fext, suffix

    if (.not.present (fext)) then 
   	    fext = ext
    end if

	if (present (suffix)) then
		open(unit = 1, file = filepath//filename//trim(suffix)//fext, form = 'formatted', access = 'sequential', &
        & status = 'old')
        read(1,*) data_slm
        close(1)
	else
		open(unit = 1, file = filepath//filename//fext, form = 'formatted', access = 'sequential', &
        & status = 'old')
        read(1,*) data_slm
        close(1)
	endif

   end subroutine read_txt
   
end module io_mod

