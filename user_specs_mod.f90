
! For all variables except SL, the prefix d- is equivalent to δ- in Kendall, i.e. incremental change over one timestep.
! The prefix delta- is equivalent to Δ, i.e. total change since first timestep.
! SL is different - dSL refers to spatially heterogeneous change only (the script SL in Kendall), and deltaSL is the 
! total change including eustatic - both since first timestep
! The suffix -xy denotes a quantity in the spatial domain; -lm denotes a quantity in the spectral domain.
! The prefix old- generally refers to the previous timestep, used to calculate increments, although sometimes 
! to previous iteration to check convergence.
! The suffix -0 generally refers to the value at first timestep, used to calculate total changes.
! ninner refers to the inner loop within each timestep to converge the ocean load.
! n refers to the timestep, between 1 and TW_nfiles.
! Other than the above, I have endeavored to be faithful to the notation of Kendall et al.

! The INPUT directory should contain the following files:
!  ********************************************************************************************************************
!  *!!! CAUTION: * Ice files: in this code, ice files are designed to start with a number '0'. eg) 'iceload0',        *
!                  'iceload1'...'iceload#'. Each ice file contains ice thickness in metres on a Gauss-Legendre grid.  *   
!                 and an interval between each ice file should be equal to the value of of the variable 'dt1'         *
! * One needs to be cautious on modifying (if has to) the main code to accomodate the right ice file number           *
!  ******************************************************************************************************************** 
!  
!  - Times: named 'times', containing the times, in years, of each ice files. Present is 0, past is negative, etc.
!      ** Time index starts from '1' (not '0')- i.e. At times(1) there has been no ice melting episode yet. 
!         i.e. The first melting episode from 'iceload0' and 'iceload1' happens between times(1) and times(2).
!  - Present-day bedrock topography: named defined in the "topomodel" variable in the user_specs_mod module, with 
!     values in metres on a Gauss-Legendre grid,
!  - Love numbers in JXM's maxwell.f output format: name is specified in the "planetmodel" variable in user_specs_mod.

! The OUTPUT directory will contain the following files:
!  - OUTPUT FILES 'ocean#', 'beta#' and 'tgrid#' start with a suffix '0'- where '0' means no melting episode 
!    eg) ocean0, beta0 files represent ocean function and beta function (respectively) at initial time times(1). 
!  - 'SL#': total sea-level change (in metres) after '#' number of melting episodes from the beginning of the simulation, 
!         on Gauss-Legendre grid.
!    eg) 'SL1' represents the total sea-level change after 1 melting episode, between times(1) and times(2)
!        'SL37' represents the tital sea-level change after 37 melting episodes, between times(1) and times(38)
!  - 'dS_converged#': total ocean loading changes(in spherical harmonics) after '#' number of melting episodes
!        from the beginning of the simulation.
!  - 'times':  Times in years of each output file,
!  - 'R1','G1','R2','G2',...: Radial and Geoid components of RSL change (only if parameter calcRG is enabled).
! NOTE ON GRID FILES: All of the grid files (i.e., ice, topography, sea-level, radial, geoid) could be in ASCII text or 
!  flat binary formats. For these files, the order of values in read and write operations are as follows: Along lines 
!  of equal latitude, from 0 degrees towards 360 degrees of longitude; each of these lines are then ordered from 0 
!  colatitude towards 180 (i.e., from North to South). 


!================================================================================================USER SPECIFICATIONS===!
module user_specs_mod
!______________________________________________________________________________________________________________________!
  
   ! Directories=======================================================================================================!
   ! 'inputfolder_ice' stores ice history files (if coupled, this folder provides iceloads outside the ice model domain)
   ! 'inputfolder' stores files such as modern observed topography, times array, known initial topography
   ! 'planetfolder': Planetary model directory, input forder for the Earth structure (i.e. PREM files)
   !  The filename of the desired model (i.e., the Love numbers) given by the planetmodel variable below. This is now 
   !  incorporated into the planets_mod module below, since automation limits the freedom in naming, which could be 
   !  complex and highly customized.
   ! 'outputfolder' stores output files from the sea level model (e.g. SL#, tgrid#, beta#, ocean#, dS_converged#, TPW)
   ! 'outputfolder_ice' stores global ice cover files, combining the prescribed ice cover outside the ice domain, 
   !  and the ice cover predicted by the dynamic model. This folder is used only when the SLM is coupled to an ice model.
   ! 'folder_coupled' stores files that are exchanged between the ice (NHiceload) and sea level (bedrock) models. It is
   !  not used if the sea-level model (SLM) is not coupled to an ice sheet model (ISM)

   integer, parameter :: str_len = 200
   integer :: unit_num

   ! Input directory
   character(str_len) :: inputfolder_ice  = 'INPUT_FILES/icemodel/'
   character(str_len) :: inputfolder  = 'INPUT_FILES/others/'
   character(str_len) :: planetfolder = 'INPUT_FILES/earthmodel/'
   character(str_len) :: gridfolder = 'INPUT_FILES/others/'
	  
   ! Output directory
   character(str_len) :: outputfolder = 'OUTPUT_SLM/'
   character(str_len) :: outputfolder_ice = 'ICELOAD_SLM/'

   ! Other directory
   character(str_len) :: folder_coupled = ''
  
   ! Common file extensions============================================================================================!
   ! Since the code does not explicitly specify the precision (single or double) of the input/output files, an 
   !  extension could be used to designate the type. Alternatively, the files could also be delimited text files (to 
   !  be implemented later). Specify the common extension here (applicable to the names of ALL input/output files): 
   character(4) :: ext = ''       ! '.sgl' | '.dbl' | '.txt' | '.nc' if fType == 'netcdf'
   ! ... and their file type:
   character(6) :: fType_in  = 'text'  ! 'netcdf' | 'text'
   character(6) :: fType_out = 'text'  ! 'netcdf' | 'text' | 'both' for both netcdf and text format
   
   ! Various selection ================================================================================================!
   character(5)  :: whichplanet   = 'earth'                    ! e.g. 'earth', 'Mars', etc.
   character(str_len) :: planetmodel   = 'prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22' ! For now, this is generated from maxwell.f by JXM
   character(str_len) :: icemodel      = 'AISload_Run85_0_0_'       ! Common name of ice files in 'inputfolder_ice'
   character(str_len) :: icemodel_out  = 'iceload_out_'       ! Name of ice files in 'outputfolder_ice'
   character(str_len) :: timearray     = 'times'                    ! Name of times array text file
   character(str_len) :: topomodel     = 'etopo2_512_AISbedmap2'       ! Bedrock topography (NO ICE INCLUDED!!) at time = 0ka
   character(str_len) :: topo_initial  = 'etopo2_512_AISbedmap2'
   character(str_len) :: grid_lat       = 'GLlat_512.txt'           ! Grid file for latitude
   character(str_len) :: grid_lon       = 'GLlon_512.txt'           ! Grid file for longitude
   
   ! Model parameters==================================================================================================!
   integer, parameter :: norder = 512           ! Max spherical harmonic degree/order
   integer, parameter :: npam = 500             ! Max relaxation modes
   integer, parameter :: nglv = 512             ! Number of GL points in latitude
   real, parameter :: epsilon1 = 1.0E-5         ! Inner loop convergence criterion
   real, parameter :: epsilon2 = 1.0E-5         ! Outer loop convergence criterion 
                                                !  (if doing a convergence check for outer loop, see below)

   ! CHECK TRUE OR FALSE ==============================================================================================!
   logical :: checkmarine = .false.  ! .true. to check for floating marine-based ice
                                                ! .false. to assume all ice is grounded
   logical :: tpw = .true.           ! .true. to incorporate rotational feedback								                                                                                   ! .false. for non-rotating planet
   logical :: calcRG = .true.       ! .true. to calculate the radial and geoid displacements; note that  
                                                !    the "true" option only works for a fixed number of outer loops
                                                !    (i.e., no convergence checks!).
                                                ! .false. to only calculate RSL.
   logical :: input_times = .false.  ! .true. if time array is provided from an existing text file                                                                    ! .false. if timearray is calculated within the main code                              
   logical :: initial_topo = .true. ! .true. initial topo is known
                                                ! .false. the code assumes initial topography is equal to modern 
                                                !      topography file "truetopo"
   logical :: iceVolume = .true.     ! .true. to output ice volume at each time step
   logical :: coupling = .false.      ! .true. if the SLM is coupled to the ISM
                                                ! .false. if not coupled                                 
   logical :: patch_ice = .false.    ! .true. patch ice data with zeros
                                                ! .false. merge the icemodel files with ice grids provided by the ISM
                                                !. patch_ice is only activated when 'coupling' is .true.

   !Time Window parameters=======================================================================================!

   !                      |--------- total length of a time window------------|
   
   !  schematic diagram   |--------dt4--------|------dt3------|---dt2---|-dt1-|
   !  of a timewindow          past                                        current time step


   ! if you would like a forward simulation WITHOUT a timewindow, simply set 'L_sim' equal to 'Ldt1',
   ! and set Ldt2, Ldt3 and Ldt4 to 0. 

   integer:: L_sim = 20! total length of a simulation, in years
   
   !internal time step intervals (dt's cannot be set as 0 but Ldt's can be)
   !**NOTE** dt# values should be defined such that dt#/dt1 is a positive integer
   integer :: dt1 = 5! the finest time interval in the TW (in years), usually equal to coupling time step
   integer :: dt2 = 10!
   integer :: dt3 = 10!
   integer :: dt4 = 10!

   integer :: Ldt1 = 20! total length of time over which dt1 covers
   integer :: Ldt2 = 0!
   integer :: Ldt3 = 0!
   integer :: Ldt4 = 0!

							
end module user_specs_mod
