!===================================================== ==============PHYSICSAL & MATHEMATICAL CONSTANTS===!
module constants_mod
!_________________________________________________________________________________________________________!
   real, parameter :: pi = 3.1415926535898      ! Pi
   complex, parameter :: ii = (0.0,1.0)           ! Square root of -1
   real, parameter :: gravConst = 6.67408E-11   ! Gravitational constant (m^3/kg/s^2)
end module constants_mod

!===========================================================================================PLANETS_MOD===!
module planets_mod
!_________________________________________________________________________________________________________!
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

   ! -------------------------------------------------------------------------
   subroutine earth_init

      radius = 6.371229E6           ! Radius of the Earth (m)
      mass = 5.976E24               ! Mass of the Earth (kg)
      rhoi = 910.0                  ! Density of ice (kg/m^3)
      rhow = 1000.0                 ! Density of fresh water (kg/m^3)
      gacc = 9.80665                ! Acceleration due to gravity at the Earth's surface (m/s^2)
      omega = 7.292e-5              ! Rotation rate of the Earth (rad/s)
      moiA = 0.3296145*mass*radius**2 ! Principal moment of inertia of the Earth
      moiC = 0.3307007*mass*radius**2 ! Principal moment of inertia of the Earth
      kf = 0.9342+0.008             ! Fluid (Tidal) Love number

   end subroutine earth_init


   ! -------------------------------------------------------------------------
   subroutine mars_init

      radius = 3.3899E6             ! Radius of Mars (m)
      mass = 6.4185E23              ! Mass of Mars (kg)
      rhoi = 1220.0                 ! Density of ice mix on Mars (kg/m^3)
      rhow = 1000.0                 ! Density of fresh water (kg/m^3)
      gacc = 3.713                  ! Acceleration due to gravity at Mars's surface (m/s^2)
      omega = 7.08819118E-5         ! Rotation rate of Mars (rad/s)
      moiA = 0.363914*mass*radius**2  ! Principal moments of inertia of Mars
      moiC = 0.365905*mass*radius**2  ! Principal moments of inertia of Mars
      !----------------------------------------------------------------------------------!
      ! Lithospheric thickness and corresponding kf (based on Zharkov and Gudkova, 2005)
      !    15    |    48    |    58    |    84    |    110   |    164    |    200    [km]
      ! 1.203959 | 1.127384 | 1.110566 | 1.065899 | 1.023186 | 0.9458358 | 0.8986673
      !----------------------------------------------------------------------------------!
      kf = 0.899                    ! Fluid (Tidal) Love number

   end subroutine mars_init

end module planets_mod


!=====================================================================================netCDF I/O mod===!
module sl_io_mod
!-------------------------------------------------------------------------------------------------------
   use user_specs_mod, only: nglv, ext, fType_in, fType_out, outputfolder, gridfolder, grid_lat, grid_lon, unit_num, str_len
   use netCDF
   implicit none

   public :: sl_drive_readnl, sl_readnl

   contains

   subroutine check(status)

      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         print *, "NetCDF I/O error occurred in the Sea Level Model."
         print *, "Backtrace of calling routines:"
         !CALL BACKTRACE
         stop "Stopping Sea Level Model."
      endif

   end subroutine check


   ! -------------------------------------------------------------------------
   subroutine read_sl(data_slm, filename, filepath, suffix, fext)
      character (len = *), intent (in) :: filename, filepath
      real, dimension(nglv,2*nglv), intent(inout) :: data_slm
      character (len = *), optional :: suffix, fext

      if (fType_in == 'text') then
         call read_txt(data_slm, filename, filepath, suffix, fext)
      elseif (fType_in == 'netcdf') then
         call read_nf90(data_slm, filename, filepath, suffix, fext)
      else
         write(unit_num,*) "Error: fType_in should be one of 'netcdf' or 'text'"
         stop
      endif

   end subroutine read_sl


   ! -------------------------------------------------------------------------
   subroutine write_sl(data_slm, filename, filepath, suffix, fext)
      character (len = *), intent (in) :: filename, filepath
      real, dimension(nglv,2*nglv), intent(in) :: data_slm
      character (len = *), optional :: suffix, fext

      if ((fType_out == 'text') .or. (fType_out == 'both')) then
         call write_txt(data_slm, filename, filepath, suffix, fext)
      endif
      if ((fType_out == 'netcdf') .or. (fType_out == 'both')) then
         call write_nf90(data_slm, filename, filepath, suffix, fext='.nc')
      endif
      if ((fType_out /= 'netcdf') .and. (fType_out /= 'text') .and. (fType_out /= 'both')) then
         write(unit_num,*) "Error: fType_out should be one of 'netcdf', 'text', or 'both'"
         stop
      endif

   end subroutine write_sl


   ! -------------------------------------------------------------------------
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
      open(unit = 201, file = trim(gridfolder)//trim(grid_lat), form = 'formatted', &
      & access = 'sequential', status = 'old')
      read(201,*) latgrid
      close(201)

      open(unit = 201, file = trim(gridfolder)//trim(grid_lon), form = 'formatted', &
      & access = 'sequential', status = 'old')
      read(201,*) longrid
      close(201)

       !write out data

      !create file
      if (present (suffix)) then
          if (present (fext)) then
             call check( nf90_create(trim(filepath)//trim(filename)//trim(suffix)//trim(fext), &
            & nf90_clobber, ncid) )
          else
             call check( nf90_create(trim(filepath)//trim(filename)//trim(suffix)//trim(ext), &
            nf90_clobber, ncid) )
          endif
      else
         if (present (fext)) then
            call check( nf90_create(trim(filepath)//trim(filename)//trim(fext), nf90_clobber, &
            & ncid) )
         else
            call check( nf90_create(trim(filepath)//trim(filename)//trim(ext), nf90_clobber, &
            & ncid) )
         endif
      endif

      call check( nf90_def_dim(ncid, 'lon', nglv*2, lon_dimid)  ) ! Define the dimensions of the griddata
      call check( nf90_def_dim(ncid, 'lat', nglv,   lat_dimid)  )

      dimids =  (/ lon_dimid, lat_dimid /)                        ! Define the dimension of the variable

      call check( nf90_def_var(ncid, trim(filename), nf90_double, dimids, varid)) ! Define variable
      call check( nf90_def_var(ncid, 'lat', nf90_double, lat_dimid, lat_varid))
      call check( nf90_def_var(ncid, 'lon', nf90_double, lon_dimid, lon_varid))
      call check( nf90_enddef(ncid)) ! End definition

      call check( nf90_put_var(ncid, varid, reshape(data_slm,[2*nglv,nglv]))) !write data
      call check( nf90_put_var(ncid, lon_varid, longrid))
      call check( nf90_put_var(ncid, lat_varid, latgrid(nglv:1:-1)))
      call check( nf90_close(ncid))

   end subroutine write_nf90


   ! -------------------------------------------------------------------------
   subroutine read_nf90(data_slm, filename, filepath, suffix, fext)

      character (len = *), intent(in) :: filename, filepath !file name and path
      real, dimension(2*nglv,nglv) :: data_temp !temp. variable name in the SLM in which nc data will be stored
      real, dimension(nglv,2*nglv), intent(out) :: data_slm
      character (len = *), optional ::  suffix
      character (len = *), optional :: fext
      integer :: ncid, varid
      character(str_len) :: fullfilename

      !open the file
      if (present (suffix)) then
          if (present (fext)) then
             fullfilename = trim(filepath)//trim(filename)//trim(suffix)//trim(fext)
          else
             fullfilename = trim(filepath)//trim(filename)//trim(suffix)//trim(ext)
          endif
      else
         if (present (fext)) then
            fullfilename = trim(filepath)//trim(filename)//trim(fext)
         else
            fullfilename = trim(filepath)//trim(filename)//trim(ext)
         endif
      endif

      write(unit_num,*) 'Reading variable ', trim(filename), ' from netcdf file ', trim(fullfilename)
      call check( nf90_open(trim(fullfilename), nf90_nowrite, ncid) )

      call check( nf90_inq_varid(ncid, trim(filename), varid) ) !get varid of the data variable
      call check( nf90_get_var(ncid, varid, data_temp) ) ! read the data
      call check( nf90_close(ncid) ) ! close the file
      data_slm = reshape(data_temp,[nglv,2*nglv])

   end subroutine read_nf90


   ! -------------------------------------------------------------------------
   subroutine check_ncFile(ncvar, slmvar, varname)

      real, dimension(nglv,2*nglv), intent(in) :: ncvar, slmvar
      character (len = *), intent(in) :: varname

      if (sum(sum(ncvar, dim=1))-sum(sum(slmvar,dim=1)) .NE.0) then
         write(201,*) 'ncFile benchmark test FAILED. Var name is:                              ', varname
      else if (sum(sum(ncvar, dim=1))-sum(sum(slmvar, dim=1)).eq.0) then
         write(201,*) 'ncFile benchmark test PASSED. Var name is:                              ', varname
      endif

   end subroutine check_ncFile


   ! -------------------------------------------------------------------------
   subroutine write_txt(data_slm, filename, filepath, suffix, fext)

   real, dimension(nglv,2*nglv), intent(in) :: data_slm
      character (len = *), intent(in) :: filename, filepath
      character (len = *), optional :: fext, suffix

    if (.not.present (fext)) then
          fext = ext
    end if

   if (present (suffix)) then
      open(unit = 201, file = trim(filepath)//trim(filename)//trim(suffix)//fext, form = 'formatted', &
      & access = 'sequential', status = 'replace')
        write(201,'(ES16.9E2)') data_slm
        close(201)
   else
      open(unit = 201, file = trim(filepath)//trim(filename)//trim(fext), form = 'formatted', &
      & access = 'sequential', status = 'replace')
        write(201,'(ES16.9E2)') data_slm
        close(201)
   endif

   end subroutine write_txt


   ! -------------------------------------------------------------------------
   subroutine read_txt(data_slm, filename, filepath, suffix, fext)

   real, dimension(nglv,2*nglv) :: data_slm
      character (len = *), intent(in) :: filename, filepath
      character (len = *), optional :: fext, suffix

    if (.not.present (fext)) then
          fext = ext
    end if

   if (present (suffix)) then
      open(unit = 201, file = trim(filepath)//trim(filename)//trim(suffix)//trim(fext), &
      & form = 'formatted', access = 'sequential', status = 'old')
      read(201,*) data_slm
      close(201)
   else
      open(unit = 201, file = trim(filepath)//trim(filename)//trim(fext), form = 'formatted', &
      & access = 'sequential', status = 'old')
      read(201,*) data_slm
      close(201)
   endif

   end subroutine read_txt


   ! -------------------------------------------------------------------------
   ! read variables needed to drive the SLM from namelist
   subroutine sl_drive_readnl(itersl, dt1, starttime)

      real, intent(out) :: starttime
      integer, intent(out) :: itersl, dt1

      namelist /time_config/ starttime, itersl, dt1

      open(201, file='namelist.sealevel', status='old')
      read(201, time_config)
      close(201)

   end subroutine sl_drive_readnl


   ! -------------------------------------------------------------------------
   ! read all other variables from namelist
   subroutine sl_readnl(inputfolder_ice, inputfolder, &
                        planetfolder, gridfolder, &
                        outputfolder, outputfolder_ice, &
                        folder_coupled, ext, fType_in, fType_out, &
                        planetmodel, icemodel, icemodel_out, &
                        timearray, topomodel, topo_initial, &
                        grid_lat, grid_lon, checkmarine, &
                        tpw, calcRG, input_times, &
                        initial_topo, iceVolume, coupling, &
                        patch_ice, L_sim, dt1, dt2, &
                        dt3, dt4, Ldt1, Ldt2, &
                        Ldt3, Ldt4, whichplanet)

      character(str_len), intent(out) :: inputfolder_ice
      character(str_len), intent(out) :: inputfolder
      character(str_len), intent(out) :: planetfolder
      character(str_len), intent(out) :: gridfolder
      character(str_len), intent(out) :: outputfolder
      character(str_len), intent(out) :: outputfolder_ice
      character(str_len), intent(out) :: folder_coupled

      character(4), intent(out) :: ext
      character(6), intent(out) :: fType_in
      character(6), intent(out) :: fType_out

      character(str_len), intent(out) :: planetmodel
      character(str_len), intent(out) :: icemodel
      character(str_len), intent(out) :: icemodel_out
      character(str_len), intent(out) :: timearray
      character(str_len), intent(out) :: topomodel
      character(str_len), intent(out) :: topo_initial
      character(str_len), intent(out) :: grid_lat
      character(str_len), intent(out) :: grid_lon

      logical, intent(out)  :: checkmarine
      logical, intent(out)  :: tpw
      logical, intent(out)  :: calcRG
      logical, intent(out)  :: input_times
      logical, intent(out)  :: initial_topo
      logical, intent(out)  :: iceVolume
      logical, intent(out)  :: coupling
      logical, intent(out)  :: patch_ice

      integer, intent(out)  :: L_sim
      integer, intent(out)  :: dt1
      integer, intent(out)  :: dt2
      integer, intent(out)  :: dt3
      integer, intent(out)  :: dt4
      integer, intent(out)  :: Ldt1
      integer, intent(out)  :: Ldt2
      integer, intent(out)  :: Ldt3
      integer, intent(out)  :: Ldt4

      character(5), intent(out) :: whichplanet

      namelist /io_directory/ inputfolder_ice, inputfolder, &
                              planetfolder, gridfolder, &
                              outputfolder, outputfolder_ice, &
                              folder_coupled

      namelist /file_format/ ext, fType_in, fType_out

      namelist   /file_name/ planetmodel, icemodel, icemodel_out, &
                           timearray, topomodel, topo_initial, &
                           grid_lat,grid_lon

      namelist /model_config/ checkmarine, tpw, calcRG, &
                              input_times, initial_topo, iceVolume, &
                              coupling, patch_ice

      namelist /timewindow_config/ L_sim, dt1, dt2, dt3, &
                                   dt4, Ldt1, Ldt2, Ldt3, &
                                   Ldt4

      namelist /others/ whichplanet

      open(201, file='namelist.sealevel', status='old', form='formatted')
      read(201, io_directory)
      read(201, file_format)
      read(201, file_name)
      read(201, model_config)
      read(201, timewindow_config)
      read(201, others)
      close(201)

   end subroutine sl_readnl

end module sl_io_mod
