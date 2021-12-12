! SL_MODEL_FWTW.90  - Holly Kyeore Han (PhD Student, McGill University 2015-2021, Advisor: by Natalya Gomez),
! Sea Level Model with ForWard and TimeWindow algorithm (FWTW). LAST UPDATE: April 1st, 2021 by Holly Han 

! The is a FORWARD sea-level model with the timewindow algorithm. The code is modified from SL_TPW.f90, a new, 
! benchmarked ice-age sea-level model written by Sam Goldberg, Harvard University EPS '16 (Advised by Jerry Mitrovica) 
! based on pseudo spectral algorithm for ice-age sea-level calculations (Kendall et al., 2005). SL_TPW.f90 was then
! modified to include text format of input and output files and the Structure of Mars by Erik. N.-H. Chan (Post Doctoral
! Fellow, McGill University 2017-2019). The code was then modified the to a forward sea-level model (following the 
! algorithm by Gomez et al., 2010) that can be coupled to a dynamic ice sheet model and with a new time window algorithm 
! by Holly Han during her PhD work at McGill University.   


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


include 'spharmt.f90' ! Spherical harmonic transform module

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

   ! Input directory
   character(*), parameter :: inputfolder_ice  = 'INPUT_NHIS2GC/'
   character(*), parameter :: inputfolder  = '/project/ctb-ng50/Han/INPUT_FILES/TOPOFILES/'
   character(*), parameter :: planetfolder = '/project/ctb-ng50/Han/INPUT_FILES/PREMFILES/'   
   
   ! Output directory
   character(*), parameter :: outputfolder = 'OUTPUT_SLM/' 
   character(*), parameter :: outputfolder_ice = 'ICELOAD_SLM/'

   ! Other directory
   character(*), parameter :: folder_coupled = '' 
  
   
   ! Various selection ================================================================================================!
   character(4), parameter :: ext = ''    ! '.txt' | ''   ! Common file extension
   character(*), parameter :: whichplanet   = 'earth'                  ! e.g. 'earth', 'Mars', etc.
   character(*), parameter :: planetmodel   = 'prem_coll_512.l120C.ump5.lm5' ! For now, this is generated from maxwell.f by JXM
   character(*), parameter :: icemodel      = 'iceload'             ! Common name of ice files in 'inputfolder_ice'
   character(*), parameter :: icemodel_out  = 'iceload'          ! Name of ice files in 'outputfolder_ice'
   character(*), parameter :: timearray     = 'times'                  ! Name of times array text file
   character(*), parameter :: topomodel     = 'etopo2_512_ice6gC '       ! Bedrock topography (NO ICE INCLUDED!!) at time = 0ka       
   character(*), parameter :: topo_initial  = 'etopo2_512_ice6gC' 
   
   ! Model parameters==================================================================================================!
   integer, parameter :: norder = 512           ! Max spherical harmonic degree/order
   integer, parameter :: npam = 500             ! Max relaxation modes
   integer, parameter :: nglv = 512             ! Number of GL points in latitude
   real, parameter :: epsilon1 = 1.0E-5         ! Inner loop convergence criterion
   real, parameter :: epsilon2 = 1.0E-5         ! Outer loop convergence criterion 
                                                !  (if doing a convergence check for outer loop, see below)

   ! CHECK TRUE OR FALSE ==============================================================================================!
   logical, parameter :: checkmarine = .false.  ! .true. to check for floating marine-based ice
                                                ! .false. to assume all ice is grounded
   logical, parameter :: tpw = .true.           ! .true. to incorporate rotational feedback								                                                                                   ! .false. for non-rotating planet												
   logical, parameter :: calcRG = .false.       ! .true. to calculate the radial and geoid displacements; note that  
                                                !    the "true" option only works for a fixed number of outer loops 
                                                !    (i.e., no convergence checks!).
                                                ! .false. to only calculate RSL. 
   logical, parameter :: input_times = .false.  ! .true. if time array is provided from an existing text file                                                                                              ! .false. if timearray is calculated within the main code                              
   logical, parameter :: initial_topo = .true. ! .true. initial topo is known
                                                ! .false. the code assumes initial topography is equal to modern 
                                                !      topography file "truetopo"
   logical, parameter :: iceVolume = .true.     ! .true. to output ice volume at each time step
   logical, parameter :: coupling = .false.      ! .true. if the SLM is coupled to the ISM
                                                ! .false. if not coupled                                 
   logical, parameter :: patch_ice = .false.    ! .true. patch ice data with zeros
                                                ! .false. merge the icemodel files with ice grids provided by the ISM
                                                !. patch_ice is only activated when 'coupling' is .true.
                                
   !Time Window parameters=======================================================================================!

   !                      |--------- total length of a time window------------|
   
   !  schematic diagram   |--------dt4--------|------dt3------|---dt2---|-dt1-|
   !  of a timewindow          past                                        current time step


   ! if you would like a forward simulation WITHOUT a timewindow, simply set 'L_sim' equal to 'Ldt1',
   ! and set Ldt2, Ldt3 and Ldt4 to 0. 

   integer, parameter :: L_sim = 21000! total length of a simulation, in years
   
   !internal time step intervals (dt's cannot be set as 0 but Ldt's can be)
   !**NOTE** dt# values should be defined such that dt#/dt1 is a positive integer
   integer, parameter :: dt1 = 200! the finest time interval in the TW (in years), usually equal to coupling time step
   integer, parameter :: dt2 = 0!  
   integer, parameter :: dt3 = 0!
   integer, parameter :: dt4 = 0! 
   
   integer, parameter :: Ldt1 = 21000! total length of time over which dt1 covers 
   integer, parameter :: Ldt2 = 0! 
   integer, parameter :: Ldt3 = 0!
   integer, parameter :: Ldt4 = 0!


   
end module user_specs_mod

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

!=======================================================================================================================!
!                                                      MAIN BLOCK                                                       !
!_______________________________________________________________________________________________________________________!

!=======================================================================================================MAIN PROGRAM====!
program sl_model
!_______________________________________________________________________________________________________________________!
use spharmt
use user_specs_mod
use planets_mod
implicit none

!=======================================================================================================================!
!                                             VARIABLES                                                                 !
!________________________________________________(Edit with caution)____________________________________________________!

!===============================  Variables for ice sheet - sea level model coupling ===================================|
real, dimension(nglv,2*nglv) :: nh_bedrock        ! Northern Hemispheric bedrock provided by the ice sheet model        |
real, dimension(nglv,2*nglv) :: nh_iceload        ! Northern Hemispheric iceload provided by the ice sheet model        |
!=======================================================================================================================|

!============================================  Variables for the time window============================================|
integer :: L_TW, TW_nfiles                     !    L_TW; length of a full time window                                  |
                                               !    TW_nfiles: total number of ice files in a full time window          |
integer :: ncalls, TW_nmelt, difference, ice   !    TW_nmelt: total number of melting episodes in the TW                |
integer :: iter_TW                             !    nmelt_TW:  index variable for melting episodes in the TW            |
integer :: Travel, Travel_total                !    number of travels that the timewindow has made                      |
                                               !    Travel_total: total number of travelling                            |
integer :: masked_iceload                      !    elements in mask and iceload arrays                                 |
integer, dimension(:), allocatable :: mask, iceload, icefiles!                                                          |
integer, dimension(:), allocatable :: TIMEWINDOW        !   TW array                                                    |
integer :: err                                 !   I/O error                                                            |
integer :: dummy                               ! dummy variable                                                         |
integer, dimension(4) :: int_dt, Rdt, Ndt, Ldt         !   Values of internal timesteps of the time window              |
                                          !    Rdt: ratio between the smallest internal timestep (dt1)                  |
                                          !      to other internal time steps                                           |
                                          !    Ndt: total number of each internal timesteps                             |
                                          !    Ldt: an array for the time covered by each internal timestep             |
                                          !                                                                             |
!=======================================================================================================================|

!======================================== Variables for dynamic arrays==================================================|
real, dimension(:), allocatable :: times        ! Timesteps of ice model (years)                                        |
real, dimension(:,:,:), allocatable :: icexy    ! Spatial inputs of ice model                                           |   
real, dimension(:,:,:), allocatable :: sl       ! Big arrays of sea level change                                        |
complex, dimension(:,:,:), allocatable:: dicestar,dS,deltaicestar,deltaS ! Big arrays of changes in loads               |
                                                                                       !  used in Love number viscous   | 
                                                                                       !  response                      |
real, dimension(:,:,:), allocatable :: rr,gg                         !  R and G (radial displacement and geoid change)  |
complex, dimension(:,:,:), allocatable :: dlambda, deltalambda       ! Big arrays of changes in rotational driving      |
real, dimension(:,:,:), allocatable :: dil                           ! Big array of changes in IL                       |
real, dimension(:,:), allocatable :: dm                              ! Big array of changes in m                        |
real, dimension(:), allocatable :: lovebetatt, lovebetattrr          ! Used in Love number calculations                 |
real, dimension(:,:), allocatable :: lovebeta                        !                                                  |
real, dimension(:,:), allocatable ::   lovebetarr                    ! Used in Love number calculations                 |
!=======================================================================================================================|

!========================= Variables for topography correction (i.e. outer iteration)===================================|
real :: current_time                            ! time (in years) since the start of the simulation                     |
real, dimension(nglv,2*nglv) :: tinit_0         ! topography at the very beginning of the simulation (i.e. time0)       |
real, dimension(nglv,2*nglv) :: tinit_0_last    ! tinit_0 from the previous outer-iteration                             |
real, dimension(nglv,2*nglv) :: pred_pres_topo  ! predicted present topography at current outer-iteration loop          |
real, dimension(nglv,2*nglv) :: init_topo_corr  ! correction applied to compute tinit_0 at current outer-iteration loop |  
character(6) :: iterstr                         ! String for timestep number for reading/writing files                  |
!=======================================================================================================================|
! Inputs
real, dimension(nglv,2*nglv) :: truetopo        ! Present-day topography
real, dimension(npam,norder) :: rprime,r,s      ! Love numbers
real, dimension(norder) :: ke,he                ! Love numbers
real, dimension(npam,norder) :: rprimeT,rT      ! Love numbers (tidal)
real, dimension(norder) :: kTE, hTE             ! Love numbers (tidal)

! Model calculations
real, dimension(nglv,2*nglv) :: glw_matrix              ! weigtht in the Gaussian-Legendre grid 
real, dimension(nglv,2*nglv) :: deltaslxy, dslxy         ! Total sea level change, 
                                                         !  total spatially heterogeneous sea level change
real, dimension(nglv,2*nglv) :: icestarxy                ! Grounded ice thickness
real, dimension(nglv,2*nglv) :: beta, cstarxy, cstar0    ! Grounded ice mask, ice-free ocean function
real, dimension(nglv,2*nglv) :: tOxy, rOxy, tTxy         ! Projections used to calculate loads and shoreline migration
complex, dimension(0:norder,0:norder) :: cstarlm,oldcstarlm,tOlm,rOlm,dSlm,olddSlm,&
                                         icestarlm,dicestarlm,deltaicestarlm,oldicestarlm,icestar0, &
                                         t0lm,oldt0lm,tTlm,oldtTlm,dsllm,deltasllm  ! Above, in spectral domain

real :: conserv                                          ! Uniform geoid shift (ΔΦ/g)
real :: ttl, ekhl                                        ! Used in Love number calculations
complex, dimension(0:norder,0:norder) :: viscous         ! Used in Love number calculations
real :: xi, zeta                                         ! Convergence checks
real :: ice_volume                                ! ice volume 
! For calculating R and G separately
real, dimension(nglv,2*nglv) :: rrxy, drrxy_computed
complex, dimension(0:norder,0:norder) :: rrlm, dgglm, drrlm_computed

complex :: viscousrr

! Rotation calculations
! real, parameter :: --------------------------------->     ! Fluid Love number is defined in the planetary modules
complex, dimension(0:2,0:2) :: lambda, oldlambda, lambda0   ! Rotational driving
real, dimension(3) :: mm, oldm, sum_m                       ! Rotational perturbation vector
real, dimension(3,3) :: il, oldil, sum_il                     ! Load perturbations to MOI
complex, dimension(0:2,0:2) :: dsl_rot, rr_rot              ! Sea level change from rotational feedbacks
real :: ekhTE                                               ! Used in Love number calculations
real :: betatt, betattprime                                 ! Used in Love number calculations
complex :: viscoustt,viscousttrr                            ! Used in Love number calculations

! Miscellaneous variables
integer :: ninner                                           ! Iteration of inner loop
integer :: i,j,k,l,m,n,nn                                   ! Do-loop indices
character(6) :: numstr, numstr2                             ! String for timestep number for reading/writing files
integer :: counti, countf,countrate         ! Computation timing
real :: counti_cpu, countf_cpu
type(sphere) :: spheredat                                   ! SH transform data to be passed to subroutines

! For Jerry's code to read in Love numbers
integer :: legord(norder),nmod(norder),nmodes(norder),ll,nm,np
real :: xn
real, dimension(3,norder) :: elast,asymv,telast,tasymv
real, dimension(npam,norder) :: resh,resl,resk,tresh,tresl,tresk
real :: taurr,taurt,dmx,dmy

real, dimension(nglv,2*nglv) :: beta0, cxy0, cxy
real, dimension(nglv,2*nglv) :: topoxy, topoxy_m1, tinit
                                              ! topoxy_m1: topogramy from the previous timestep (m1: minus one)
                                              ! topoxy: topography at the currect timestep
                                              ! tinit: initial topography within the TimeWindow
                                     
complex, dimension(0:norder,0:norder) :: dS_converged
integer :: nmelt, nfiles, iter, itersl, dtime           ! Number of melting episodes up to the current timestep
real :: starttime                                       ! start time of the simulation 
integer :: iargc, nargs                                 ! Arguments read in from a bash script                              
character(16) :: carg(20)                               ! Arguments from a bash script
character(3) :: skip                                    ! variable used to skip lines in reading TPW file 


! Planetary values
if (whichplanet == 'earth' .or. whichplanet == 'Earth' .or. whichplanet == 'EARTH') then
   call earth_init
elseif (whichplanet == 'mars' .or. whichplanet == 'Mars' .or. whichplanet == 'MARS') then
   call mars_init
else
   write(*,*) 'The parameters for the planet you entered are not built in.' 
   write(*,*) 'Please check your spelling for the variable whichplanet in the user_specs_mod module.' 
   write(*,*) 'If you preferred, you could create a new subroutine for the new planet in the planets_mod module.'
   write(*,*) 'Terminating: program sl_model'
   stop
endif


! Reading in arguments from a bash script
nargs = iargc()
do i=1,nargs
   call getarg(i, carg(i))
enddo
read (carg(1),*) itersl
read (carg(2),*) iter      ! the coupling time step we are on (in years)
read (carg(3),*) dtime     ! coupling time (in years)
read (carg(4),*) starttime ! start time of the simulation (in years)

if (dtime /= dt1) then 
   write(*,*) 'dtime and dt1 should be equal to each other.'
   write(*,*) 'Please check your set up for the variables'
   write(*,*) 'Terminating: program sl_model'
   stop
endif

if (itersl.lt.1) then 
    write(*,*) 'itersl must be equal to or greather than 1'
    write(*,*) 'itersl = 1: No topography correction'
    write(*,*) 'itersl > 1: topography correction'
    write(*,*) 'Terminating: program sl_model'
    stop
endif

!##################################################################################################################
!                                       TIME WINDOW PART                                                          #
!##################################################################################################################

! save internal timestep profiles to a big array 'int_dt'
int_dt(1) = dt4
int_dt(2) = dt3
int_dt(3) = dt2
int_dt(4) = dt1

! Ldt#: total length of time over which each dt# covers 
! save into a big array
Ldt(1) = Ldt4
Ldt(2) = Ldt3
Ldt(3) = Ldt2
Ldt(4) = Ldt1

if (sum(Ldt)>L_sim) then
    write(*,*) 'Total length of simulation CANNOT be smaller than the total lengths of the internal time windows!!'
    write(*,*) 'Make sure the sum(LdT) < L_sim !!'
    write(*,*) 'Terminating: program sl_model'
    stop
endif
! Compute various variables
L_TW = 0                    ! Initializing the length of the timewindow
TW_nmelt = 0                 ! Initializing the total number of melting episodes within a FULL time window

do i=1,size(int_dt)
   Rdt(i) = int_dt(i)/dt1        ! e.g. dt1/dt1, dt2/dt1, dt3/dt1, dt4/dt1      
   Ndt(i) = Ldt(i)/int_dt(i)     ! the number of each internal timestep
   L_TW =  L_TW + Ldt(i)         ! total length of the time window
   TW_nmelt = TW_nmelt + Ndt(i)  ! total number of melting episodes in the full timewindow
enddo

ncalls = L_TW/dt1        ! total number of melting episodes (i.e. iterations) in case of uniform timesteps within the TW 
TW_nfiles = TW_nmelt + 1  ! total number of ice files within in the FULL time window  
Travel_total = (L_sim - L_TW)/dt1 !total number of marching steps that will be taken by the TW to complete the full simulation

if (Travel_total == 0) then 
   write(*,*) 'The length of a total simulation and the length of timewindow are the same.'
   write(*,*)' TW will not march forward!'
endif
! write(*,'(A,I4)') 'ncalls=', ncalls
! write(*,*) 'internal timesteps', int_dt
! write(*,*) 'Rdt', Rdt
! write(*,*) 'Ldt', Ldt
! write(*,*) 'Ndt', Ndt
write(*,*) 'total length of time window (L_TW)', L_TW
! write(*,*) 'ncalls', ncalls
write(*,*) 'total number of melting episode within a Full TW (TW_nmelt)', TW_nmelt
write(*,*) 'total number of icefiles in a  FULL TW (TW_nfiles)', TW_nfiles 



if (iter.LT.ncalls) then !while the time window is growing, the time window has not started travelling when
   Travel = 0 ! 
   if (iter == 0) then 
      !nmelt = 0
      nfiles = 1
   endif
else
   Travel = iter - ncalls ! This is the number of marching steps the full TW has taken upto the current timestep
   !nmelt = TW_nmelt       ! Number of melting episondes in a full time window
endif

! allocate dimensions to arrays
allocate (mask(ncalls+1),iceload(ncalls+1),icefiles(ncalls+1))
allocate (TIMEWINDOW(TW_nfiles))      

! Initialize and grow the time window when t>0 but Travel == 0    
! initialize the mask array and icefile numbers for the timewindow
do i = 1, ncalls+1
    iceload(i) = i-1 ! iceload number calculated from ice sheet model
    mask(i) = 0      ! set all elements of the mask vector to be zero 
end do 
mask(1) = 1          ! the first element of the mask will always be one

! fill in the mask array for the time window
i = 1
k = 1
do while (i.LE.size(int_dt))
   do j = 1, Ndt(i)
      m =  k + Rdt(i)
      mask(m) = 1
      k = m
   enddo
   i=i+1   
enddo

! Multiply the mask created above to an array of iceloads
!    & to find out which ice files to read in
i=1
do while (i.LE.ncalls+1)
    if (i.LE.iter+1) then
       j=1
       do while (j.LT.i+1)
          difference =ncalls+1-(i-j)
          icefiles(j) = iceload(j)*mask(difference)
          j=j+1
       end do     
    elseif (i.GT.iter+1) then
       icefiles(i) = 0
    endif    
    i=i+1
end do

! Output the ice file numbers that are called in a simulation over the FIRST TIMEWINDOW
k = 2
do j=1,ncalls+1 !iter_end is equivalent to ncalls
    ice=icefiles(j)

    if (j.LE.1) then
       TIMEWINDOW(1) = 0
   endif
    if (j.GT.1 .AND. ice.GT.0) then
       TIMEWINDOW(k) = icefiles(j)
       k = k+1
       m = k
       do while (m.LT.TW_nmelt+1) 
          TIMEWINDOW(m) = 0
          m = m + 1
       end do
    endif
end do
nmelt = k - 2    ! Update the number of metling episode  nmelt
nfiles = nmelt+1 ! Number of icefiles read in at each time the sea-level model is called

! Find ice files to read in once the time window has fully grown and started travelling.
if (iter==ncalls .OR. Travel.GT.0) then
   k=2
    do i=2,ncalls+1
       TIMEWINDOW(1) = 0
      masked_iceload = mask(i)*iceload(i)

       if (masked_iceload.GT.0) then
          TIMEWINDOW(k)=masked_iceload
          k=k+1
       endif
    end do

   ! TW is now marching forward
   do i=1,nfiles
      TIMEWINDOW(i) = TIMEWINDOW(i) + Travel
   enddo
    write(*,*) 'Icefiles in the TW at current marching step:', TIMEWINDOW
   write(*,*) 'The TW will stop marching when TRAVEL == Travel_total:'
   write(*,*) ''
   write(*,'(A,I4,A)') '   ', Travel_total - Travel, ' more TW steps to march!'
endif


! Based on the calcualted number of files read in within the time window at each timestep, 
! allocate arrays to the below variables.
allocate (times(nfiles), lovebetatt(nfiles), lovebetattrr(nfiles))
allocate (lovebetarr(nfiles,norder),lovebeta(nfiles,norder))
allocate (icexy(nglv, 2*nglv, nfiles),sl(nglv,2*nglv,nfiles))                     
allocate (dS(0:norder,0:norder,nfiles),deltaS(0:norder,0:norder,nfiles))   
allocate (dicestar(0:norder,0:norder,nfiles), deltaicestar(0:norder,0:norder,nfiles))               
allocate (rr(nglv,2*nglv,nfiles),gg(nglv,2*nglv,nfiles))      
allocate (dil(3,3,nfiles), dlambda(0:2,0:2,TW_nfiles),deltalambda(0:2,0:2,nfiles))
allocate (dm(3,nfiles))   


!!!!!!!!!!!!!!!!!!!ATTENTION!!!!!!!!
! TIMEWINDOW = TIMEWINDOW+1
! TIMEWINDOW = -(TIMEWINDOW*100+starttime)
! CAUTION in modifying the ice file numbers

!write(6,*) 'Timewindow is growing, ice file numbers to read in :', TIMEWINDOW(1:nfiles)
write(6,*) 'number of melting episodes (nmelt):', nmelt
write(*,'(A,I4)') 'Total number of files in the current time window (nfiles) = ', nfiles
write(*,'(A,I4)') ' Number of marching steps the TW has taken =        ', Travel
write(*,*) ''

!===========================================================
!                      BEGIN EXECUTION                      
!___________________________________________________________
call cpu_time(counti_cpu)
call system_clock(count = counti, count_rate = countrate)  ! Total computation time
call spharmt_init(spheredat, 2*nglv, nglv, norder, radius) ! Initialize spheredat (for SH transform subroutines)

!-----------------------------------------------------------
!                    Model initialization & Read input files
!-----------------------------------------------------------
if (coupling) then 
    write(*,*) 'Sea level model is coupled to the ice sheet model, reading in NH_iceload'
    open(unit = 1, file = folder_coupled//'NH_iceload'//ext, form = 'formatted', access = 'sequential', &
    & status = 'old')
    read(1,*) nh_iceload
    close(1)
endif

!========================================================================================================================
!       Initialize the model at times(1) when there has been no melting episodes yet  NMELT = 0 
!========================================================================================================================

if (nmelt==0) then

    write(*,*) 'nmelt=0, INITIALIZING THE SEA LEVEL MODEL..'

    !initialize variables 
    j = TIMEWINDOW(1) ! initial file number (i.e. 0)
    write(*,'(A,I4)') 'initial file number:', j
    write(numstr,'(I4)') j
    numstr = trim(adjustl(numstr))
!     k = -1*starttime
!     write(numstr2,'(I6)') k
!     write(*,'(A,I6)') 'initial icefile suffix:', k
!     numstr2 = trim(adjustl(numstr2))
    
    !====================== topography and ice load========================
    ! read in the initial iceload from the coupled ice input folder
    open(unit = 1, file = inputfolder_ice//icemodel//trim(numstr)//ext, form = 'formatted',  &
    & access = 'sequential', status = 'old')
    read(1,*) icexy(:,:,1)
    close(1)
    
    !  Initialize topography (STEP 1)
    if (initial_topo) then   
       write(*,*) 'Reading in initial topo file'
       open(unit = 1, file = inputfolder//topo_initial//ext, form = 'formatted', access = 'sequential', &
       & status = 'old')
       read(1,*) tinit_0
       close(1)
       
    else  ! if initial topo is unknown
    
       ! Present-day observed topography
       write(*,*) 'Reading in ETOPO2 file'
       open(unit = 1, file = inputfolder//topomodel//ext, form = 'formatted', access = 'sequential', status = 'old')
       read(1,*) truetopo
       close(1)
      
       ! if topography is not iteratively getting improved, or if at the first loop of outer-iteration
       if (itersl == 1) then 
         write(*,*) 'Initial topo is unknown, set the modern-observed topo to be the initial topo'
           ! truetopo(:,:)=truetopo(:,:)-icexy(:,:,TW_nfiles) ! If using ice topography
           tinit_0(:,:) = truetopo(:,:)! (eq. 48) 
     else 
           write(*,*) 'Topographic correction is ON. Updating initial topography'
         ! read in predicted present topo from the previous outer-iteration 'itersl-1'
           write(iterstr,'(I2)') itersl-1
           iterstr = trim(adjustl(iterstr))
        
         open(unit = 1, file = outputfolder//'pred_pres_topo_'//trim(iterstr)//ext, form = 'formatted', &
           & access = 'sequential', status = 'old')
         read(1,*) pred_pres_topo
         close(1)
        
        ! read in tinit_0 from the previous outer-iteration 'itersl-1'
         open(unit = 1, file = outputfolder//'tgrid0_'//trim(iterstr)//ext, form = 'formatted', &
         & access = 'sequential', status = 'old')
         read(1,*) tinit_0_last
         close(1)
        
          ! compute topography correction for the initial topography 
          init_topo_corr(:,:) = truetopo(:,:) - pred_pres_topo(:,:)
          tinit_0(:,:) = tinit_0_last(:,:) + init_topo_corr(:,:)
     endif
    endif

    if (coupling) then  !if coupling the ICE SHEET - SEA LEVEL MODELs
       write(*,*) 'Merge initial topography with NH_bedrock and initial ice load with NH_iceload'
    
       ! Bedrock from the ice sheet model
       open(unit = 1, file = folder_coupled//'NH_bedrock'//ext, form = 'formatted', access = 'sequential', &
       & status = 'old')
       read(1,*) nh_bedrock
       close(1)
       
       ! merge intitial topography with bedrock provided by the ice sheet model.
       do j = 1,2*nglv
          do i = 1,nglv
             if (nh_bedrock(i,j) < 9999) then
                tinit_0(i,j) = nh_bedrock(i,j)
             endif
          enddo
       enddo
    
       ! for the initial ice load file
       if (patch_ice) then 
          ! add zeros to the grids that are not defined in the ice sheet model.
          do j = 1,2*nglv
             do i = 1,nglv
                if (nh_iceload(i,j) == 9999) then 
                   icexy(i,j,nfiles) = 0.0
                endif
             enddo
          enddo 
       else
          ! merge the iceload with that from the ice sheet model.
          do j = 1,2*nglv
             do i = 1,nglv
                if (nh_iceload(i,j) < 9999) then 
                   icexy(i,j,nfiles) = nh_iceload(i,j)
                endif
             enddo
          enddo 
       endif
    
       !write out the current ice load as a new file to the sea-level model ice folder
       open(unit = 1, file = outputfolder_ice//icemodel_out//trim(numstr)//ext, form ='formatted',  &
       & access = 'sequential', status = 'replace')
       write(1,'(ES16.9E2)') icexy(:,:,nfiles)     
       close(1) 
    endif ! end if (coupling)
    
    !write out the initial topo of the simulation, tgrid0 
    open(unit = 1, file = outputfolder//'tgrid'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
    & status = 'replace')
    write(1,'(ES16.9E2)') tinit_0(:,:)
    close(1)
    
    !========================== ocean function =========================
    ! Calculate an initial ocean function based on the present topography as a first guess
    do j = 1,2*nglv
       do i = 1,nglv
          if (tinit_0(i,j)<0) then
             cxy0(i,j) = 1
          else
             cxy0(i,j) = 0
          endif
       enddo
    enddo
    
    !  write out the initial ocean function as a file
    open(unit = 1, file = outputfolder//'ocean'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
    & status = 'replace')
    write(1,'(ES14.4E2)') cxy0(:,:)
    close(1)
    
    !========================== beta function =========================     
    ! calculate initial beta
    do j = 1,2*nglv
       do i = 1,nglv
          if (icexy(i,j,1)==0) then 
             beta0(i,j)=1
          else
             beta0(i,j)=0
          endif
       enddo
    enddo
    
    !  write out the initial beta function as a file
    open(unit = 1, file = outputfolder//'beta'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
    & status = 'replace')
    write(1,'(ES14.4E2)') beta0(:,:)
    close(1)
    
    !================== total ocean loading change =====================
    ! initialize the total ocean loading change and output as a file
    deltaS(:,:,1) = (0.0,0.0) 
    open(unit = 1, file = outputfolder//'dS_converged'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
    & status = 'replace')
    write(1,'(ES16.9E2)') deltaS(:,:,1)
    close(1)

    !========================== computing time =========================
    ! To write out how much time it took to compute sea-level change over one step
    ! Open a new file
    open(unit = 1, file = outputfolder//'elapsed_wall_time'//ext, form = 'formatted', access = 'sequential', &
    & status = 'replace')
    close(1)

    open(unit = 1, file = outputfolder//'elapsed_cpu_time'//ext, form = 'formatted', access = 'sequential', &
    & status = 'replace')
    close(1)
    
    !========================== time array =============================
    if (.not. input_times) then !if time array is not read in from a text file, make a new one       
        ! write a new file
        open(unit = 1, file = outputfolder//timearray//ext, form = 'formatted', access = 'sequential', &
        & status = 'replace')
        write(1,'(ES14.4E2)') starttime
        close(1)
    endif
    
    !=========================== TPW ====================================
    ! initialize the rotational components
    if (tpw) then
    
       !dil(:,:,1) = 0.0
       !dm(:,1) = 0.0
       !dlambda(:,:,1) = (0.0,0.0)
       
       il(:,:) = 0.0
       mm(:) = 0.0
       lambda(:,:) = 0.0
        
       ! write the values (0.0) for the first timestep
       open(unit = 1, file = outputfolder//'TPW'//ext, form = 'formatted', access = 'sequential', &
       & status = 'replace')
       write(1,'(9ES19.8E2/,3ES19.8E2/,18ES19.8E2)') il(:,:), mm(:), lambda(:,:)
       ! write(1,'(9ES19.8E2/,3ES19.8E2/,18ES19.8E2)') dil(:,:,1), dm(:,1), dlambda(:,:,1)
       close(1)
    endif
    !=========================ice volume================================
    if (iceVolume) then
           if (checkmarine) then
              do j = 1,2*nglv
                  do i = 1,nglv
                     if (tinit_0(i,j) > 0) then 
                     ! If not marine...
                        icestarxy(i,j) = icexy(i,j,1)
                     elseif (icexy(i,j,1) > (abs(tinit_0(i,j)) * rhow / rhoi)) then 
                     !...else if marine, but thick enough to be grounded
                        icestarxy(i,j) = icexy(i,j,1)
                     else
                     !...if floating ice
                        icestarxy(i,j) = 0
                     endif
                  enddo
               enddo
               ! Decompose ice field 
               call spat2spec(icestarxy(:,:),icestarlm(:,:),spheredat)
            else ! If not checking for floating ice
               call spat2spec(icexy(:,:,1),icestarlm(:,:),spheredat) ! Decompose ice field
            endif

        ice_volume = icestarlm(0,0)*4*pi*radius**2

        open(unit = 1, file = outputfolder//'ice_volume'//ext, form = 'formatted', access ='sequential', &
        & status = 'replace')
        write(1,'(ES14.4E2)') ice_volume
        close(1)
    endif
   
    !HH: print out the number of iteration it takes for the inner convergence
    open(unit = 1, file = outputfolder//'numiter'//ext, form = 'formatted', access ='sequential', &
    & status = 'replace')
    close(1) 

    !HH: print out the nmelt
    open(unit = 1, file = outputfolder//'nmelt'//ext, form = 'formatted', access ='sequential', &
    & status = 'replace')
    close(1)

    write(*,*) 'DONE INITIALIZATION. EXITING THE PROGRAM'
    call exit

endif !  DONE INITIALIZATION (NMELT=0)

!========================================================================================================================  

!========================================================================================================================
! Compute sea-level change associated with past ice loading changes
if (nmelt.GT.0) then
    
    if (coupling) then 
       ! Ice files
       do n = 1, nfiles-1
          j = TIMEWINDOW(n) ! icefile numbers to read in from the TW array 
          write(*,'(A,I6)') 'ice file read in from the SLM output folder, file number:', j
          write(numstr,'(I6)') j

          k = -1*starttime-j*dt1
          write(*,'(A,I7)') 'ice load, year (ago):',k
          numstr = trim(adjustl(numstr))
    
          ! read in ice files (upto the previous time step) from the sea-level model folder
          open(unit = 1, file = outputfolder_ice//icemodel_out//trim(numstr)//ext, form = 'formatted',  &
          & access = 'sequential', status = 'old')
          read(1,*) icexy(:,:,n)
          close(1)
       enddo
       
       j = TIMEWINDOW(nfiles) ! icefile number to read in from the TW array 
       write(*,'(A,I6)') 'ice file read in from the coupled input ice folder,  file number,:', j
       write(numstr,'(I6)') j
       numstr = trim(adjustl(numstr))
       
!        k = -1*starttime-j*dt1
!        write(*,'(A,I8)') 'ice load, year (ago):',k
!        write(numstr2,'(I6)') k
!        numstr2 = trim(adjustl(numstr2))        
     
       ! for iceload at the current time step, read the corresponding file from 'inputfolder_ice'
       open(unit = 1, file = inputfolder_ice//icemodel//trim(numstr)//ext, form = 'formatted',  &
       & access = 'sequential', status = 'old')
       read(1,*) icexy(:,:,nfiles)
       close(1) 
    
       if (patch_ice) then 
         ! add zeros to the grids that are not defined in the ice sheet model.
          do j = 1,2*nglv
             do i = 1,nglv
                if (nh_iceload(i,j) == 9999) then 
                   icexy(i,j,nfiles) = 0.0
                endif
             enddo
          enddo 
       else
         ! merge the iceload with that from the ice sheet model.
          do j = 1,2*nglv
             do i = 1,nglv
                if (nh_iceload(i,j) < 9999) then 
                   icexy(i,j,nfiles) = nh_iceload(i,j)
                endif
             enddo
          enddo 
       endif
    
    else ! if not coupling
      ! if sea level model is not coupled to an ice sheet model, read in iceloads from directory external_ICE
       do n = 1, nfiles
          j = TIMEWINDOW(n) ! icefile numbers to read in from the TW array 
!          write(*,'(A,I6)') 'iceload file number:', j
          write(numstr,'(I6)') j
          numstr = trim(adjustl(numstr))
!           k = -1*starttime-j*dt1
!           write(*,'(A,I6)') 'ice load, year (ago):',k
!           write(numstr2,'(I6)') k
!           numstr2 = trim(adjustl(numstr2))
          ! read in ice files (upto the previous time step) from the sea-level model folder
          open(unit = 1, file = inputfolder_ice//icemodel//trim(numstr)//ext, form = 'formatted',  &
          & access = 'sequential', status = 'old')
          read(1,*) icexy(:,:,n)
       enddo  
    
    endif !endif coupling
    
    !Time array
    if (input_times) then ! time array is inputted from an existing text file, read in and write out
       open(unit = 1, file = inputfolder//timearray//ext, form = 'formatted', access = 'sequential', status = 'old')
       open(unit = 2, file = outputfolder//timearray//ext, form = 'formatted', access = 'sequential', &
       & status = 'replace')
       read(1,*) times
       write(2,'(ES14.4E2)') times
       close(1)
       close(2)
    else ! if time array is not read in from a text file, calculate the time array 
       ! calculate times that corresponds to ice files that are read in 
       do i = 1, nfiles
          times(i) = starttime + TIMEWINDOW(i)*dt1
       !    times(i) = TIMEWINDOW(i)*100+starttime
       enddo
      ! write(*,*) 'times', times
       
       open(unit = 1, file = outputfolder//timearray//ext, form = 'formatted', access = 'sequential', &
       & status = 'old', position='append')
       write(1,'(ES14.4E2)') times(nfiles)
       close(1)         
    endif
    
    !Read in the initial topography (topo at the beginning of the full simulation)
    !This is used to output the total sea level change from the beginning of the simulation
    open(unit = 1, file = outputfolder//'tgrid0'//ext, form = 'formatted', access = 'sequential', &
    & status = 'old')
    read(1,*) tinit_0
    close(1)
    
    
    ! read in initial (first file within the time window) ocean function
    j = TIMEWINDOW(1) ! first element of the time window as the initial file
    write(*,'(A,I4)') 'file number of the first item in the TW:', j
    write(numstr,'(I4)') j
    numstr = trim(adjustl(numstr))
    
    open(unit = 1, file = outputfolder//'ocean'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
    & status = 'old')
    read(1,*) cxy0(:,:)
    close(1)
    
    if (nmelt == 1) then 
       cxy(:,:) = cxy0(:,:)
    endif
    
    ! read in initial (first file within the time window) beta  
    open(unit = 1, file = outputfolder//'beta'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
    & status = 'old')
    read(1,*) beta0(:,:)
    close(1)

    if (tpw) then
        ! read in variables for the rotation signal 
        open(unit = 1, file = outputfolder//'TPW'//ext, form = 'formatted', access = 'sequential', &
        & status = 'old')

        oldlambda(:,:) = (0.0,0.0)
        oldil(:,:) = 0.0
        oldm(:) = 0.0
       
        do n = 1, nfiles-1
            ! find the number of lines to skip to read in appropriate TPW components
            if (n==1) then
                j = TIMEWINDOW(n)
            else
                j = TIMEWINDOW(n) - TIMEWINDOW(n-1) - 1
            endif
            
            !skip lines to read in the rotational components corresponding to timesteps within the TW 
            do m = 1, j
                read(1,*) !skip line for il
                read(1,*) !skip reading in mm
                read(1,*) !skip reading in lambda
            enddo
            
            !read in TPW components - total rotational change from the beginning of simulation 
            read(1,'(9ES19.8E2)') ((il(i,j), i=1,3), j=1,3)
            read(1,'(3ES19.8E2)') (mm(i), i=1,3)
            read(1,'(18ES19.8E2)') ((lambda(i,j),i=0,2),j=0,2)
    
            !rotational changes between each time step. 
            dm(:,n) = mm(:) - oldm(:)
            oldm(:) = mm(:)
            
            dil(:,:,n) = il(:,:) - oldil(:,:)
            oldil(:,:) = il(:,:)
            
            dlambda(:,:,n) = lambda(:,:) - oldlambda(:,:)
            oldlambda(:,:) = lambda(:,:)
            
            if (n == nfiles-1) then 
                deltalambda(:,:,nfiles-1) = lambda(:,:)
            endif
          
        enddo
        close(1)
    endif !endif (TPW)
    
    ! topography from the previous timestep
    m = TIMEWINDOW(nfiles-1)
    write(numstr2,'(I4)') m
    numstr2 = trim(adjustl(numstr2))
    
    open(unit = 1, file = outputfolder//'tgrid'//trim(numstr2)//ext, form = 'formatted', access = 'sequential', &
    & status = 'old')
    read(1,*) topoxy_m1(:,:)
    close(1)
    
    if (Travel.EQ.0) then 
        !if the TW hasnt started travelling, initial topography of the TW is that of the total simulation
        tinit(:,:) = tinit_0(:,:) 
    else
        !if The TW is moving, initial topography of the TW is read in
        j = TIMEWINDOW(1) 
        write(numstr,'(I4)') j
        numstr = trim(adjustl(numstr))
    
        open(unit = 1, file = outputfolder//'tgrid'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
        & status = 'old')
        read(1,*) tinit(:,:)
        close(1)
    endif
endif ! endif nmelt>0

if (nmelt.GT.1) then
   ! read in converged ocean function from the last timestep 
   j = TIMEWINDOW(nfiles-1)  
   write(*,*) 'Reading in ocean function from previous timestep (file number)', j
   write(numstr,'(I4)') j
   numstr = trim(adjustl(numstr))
   
   open(unit = 1, file = outputfolder//'ocean'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
   & status = 'old')
   read(1,*) cxy(:,:)
   close(1)
   
   ! if there has been more than one melting episode
   ! read in the total ocean loading computed from previous timesteps 0 to nmelt minus 1.
   ! read in the of total ocean loading changes that are needed.
   do n=1, nfiles-1
      
      j = TIMEWINDOW(n)
      ! write(*,*) 'reading in converged ocean loading files'
!      write(*,*) 'sea files, j', j
      write(numstr,'(I4)') j
      numstr = trim(adjustl(numstr))
   
   
      open(unit = 1, file = outputfolder//'dS_converged'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
      & status = 'old')
      read(1,'(ES16.9E2)') dS_converged(:,:)
      close(1)
      
      deltaS(:,:,n) = dS_converged(:,:)  !save the converged total-ocean loading into a big array
      ! dS(:,:,n) = deltaS(:,:,n)-deltaS(:,:,n-1)!calculate the ocean loading changes between every time step
   enddo
endif

!-----------------------------------------------------------
!  Read in Love numbers (Jerry's output from 'maxwell.f')
!-----------------------------------------------------------
open(unit = 2, file = planetfolder//planetmodel, status = 'old')
! Following code borrowed from Jerry
read(2,*) 
do j = 1,norder
   read(2,*) legord(j), nmodes(j)
   nm = nmodes(j)
   xn = real(legord(j))
   ll = legord(j)
   nmod(ll) = nm
   read(2,*) (s(i,ll), i=1, nm)
   read(2,*) elast(1,ll), elast(2,ll), taurr, taurt, elast(3,ll)
   read(2,*) asymv(1,ll), asymv(2,ll), taurr, taurt, asymv(3,ll)
   read(2,*) (resh(i,ll), i=1,nm)
   read(2,*) (resl(i,ll), i=1,nm)
   read(2,*) (resk(i,ll), i=1,nm)
   if (xn .lt. 1.5) cycle
   read(2,*) telast(1,ll), telast(2,ll), dmx, dmy, telast(3,ll)
   read(2,*) tasymv(1,ll), tasymv(2,ll), dmx, dmy, tasymv(3,ll)
   read(2,*) (tresh(i,ll), i=1, nm)
   read(2,*) (tresl(i,ll), i=1, nm)
   read(2,*) (tresk(i,ll), i=1, nm)
enddo
! 1001  format(2i5)
! 1002  format(20a4)
! 1003  format(i10,1x,i5)
! 1020  format(5e16.8)
close(2)

! Divide by l (numbers are multiplied by l in Jerry's output)
do l = 1,norder
   resl(:,l) = resl(:,l) / real(l)
   resk(:,l) = resk(:,l) / real(l)
   tresl(:,l) = tresl(:,l) / real(l)
   tresk(:,l) = tresk(:,l) / real(l)
enddo

! Assign Love number inputs to the letters used in Kendall et al
do l = 1,norder
   he(l) = elast(1,l)
   ke(l) = elast(3,l) / real(l)
   hTE(l) = telast(1,l)
   kTE(l) = telast(3,l) / real(l)
enddo
rprime(:,:) = resk(:,:)
r(:,:) = resh(:,:)
rprimeT(:,:) = tresk(:,:)
rT(:,:) = tresh(:,:)

!===========================================================
!                       CALCULATIONS                        
!___________________________________________________________

write(*,*) ''
write(*,'(A,I4,A,EN15.4E2,A)') '  Timestep',iter,', from ',times(nfiles-1),' years to '
write(*,'(A,EN15.4E2,A)') '                     ',times(nfiles),' years'

! Decompose initial topography (STEP 2) (used to check convergence of outer loop)
call spat2spec(tinit(:,:), t0lm, spheredat)

!=====================================================================================
!                   1. BEGIN ICE PART
!===================================================================================== 

! Read in ice loads and compute the difference in iceload over each time interval
do n=1, nfiles
   ! Calculate icestar (STEP 3) (eq.43)
   if (checkmarine) then
      do j = 1,2*nglv
         do i = 1,nglv
            if (tinit(i,j) > 0) then 
            ! If not marine...
               icestarxy(i,j) = icexy(i,j,n)
            elseif (icexy(i,j,n) > (abs(tinit(i,j)) * rhow / rhoi)) then 
            !...else if marine, but thick enough to be grounded
               icestarxy(i,j) = icexy(i,j,n)
            else
            !...if floating ice
               icestarxy(i,j) = 0
            endif
         enddo
      enddo
      ! Decompose ice field 
      call spat2spec(icestarxy(:,:),icestarlm(:,:),spheredat)
   else ! If not checking for floating ice
      call spat2spec(icexy(:,:,n),icestarlm(:,:),spheredat) ! Decompose ice field
   endif
   
   if (n == 1) then
      dicestarlm(:,:) = 0.0            ! No change at first timestep
      icestar0(:,:) = icestarlm(:,:)   ! define the initial ice field
   else
      dicestarlm(:,:) = icestarlm(:,:) - oldicestarlm(:,:) ! Incremental change
   endif

   
   oldicestarlm(:,:) = icestarlm(:,:)                   ! Save to calculate increment on next time step
   deltaicestarlm(:,:) = icestarlm(:,:) - icestar0(:,:) ! Total change since time0
   dicestar(:,:,n) = dicestarlm(:,:)         ! Save into big matrix (each slice for each time step)
   deltaicestar(:,:,n) = deltaicestarlm(:,:) ! Save into big matrix (each slice for each time step)
enddo

!========================================================================================
!                            BEGIN OCEAN PART 
!========================================================================================

! Calculate beta (STEP 3) (eq. 44)
! Calculate current beta based on iceload at the current timestep
do j = 1,2*nglv
   do i = 1,nglv
       if (icexy(i,j,nfiles) < epsilon(0.0)) then 
          beta(i,j) = 1
       else
          beta(i,j) = 0
       endif
   enddo
enddo 

! calculate initial (initial within the time window) cstar
cstar0(:,:) = cxy0(:,:)*beta0(:,:)  
cstarxy(:,:)= cxy(:,:)*beta(:,:)  ! First guess to the O.F using converged O.F the preivous timestep 

call spat2spec(cstarxy,cstarlm,spheredat) ! Decompose the current cstar

! Calculate the first guess to topography correction 
tOxy(:,:)=tinit(:,:)*(cstarxy(:,:)-cstar0(:,:)) ! (eq. 70)
call spat2spec(tOxy,tOlm,spheredat) ! Decompose the topo correction term 


dS(:,:,1) = deltaS(:,:,1)
if (nmelt > 1) then
   do n = 2, nfiles-1
      dS(:,:,n) = deltaS(:,:,n)-deltaS(:,:,n-1) !Calculate the ocean loading changes between every time step
   enddo
endif

! Initial guess (eustatic) to the ocean loading change at current time step
dS(:,:,nfiles) = (cstarlm(:,:)/cstarlm(0,0)) * ((-rhoi/rhow)*dicestar(0,0,nfiles)) ! Save into a big array
dSlm(:,:) = dS(:,:,nfiles)
!========================================================================================
!                   END OF THE OCEAN PART 
!========================================================================================


!----------------------------------------------------------------------------------------
!>>>>>>>>>>>>>>>>>>> Start of inner loop <<<<<<<<<<<<<<<<<<<
!----------------------------------------------------------------------------------------
ninner = 1
do ! Inner loop
   
   !-----\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    Rotation    \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/----!
  
   if (tpw) then ! If incorporating rotation
     
       ! MOI perturbations (Matsuyama et al, 2006) (only need (1,3),(2,3), and (3,3))
     
      il(3,3) = (-8 * pi * radius**4 / (3 * sqrt(5.0))) * real(rhow * deltaS(2,0,nfiles) + rhoi * deltaicestar(2,0,nfiles))
      il(1,3) = (8 * pi * radius**4 / (sqrt(30.0))) * real(rhow * deltaS(2,1,nfiles) + rhoi * deltaicestar(2,1,nfiles))
      il(2,3) = (-8 * pi * radius**4 / (sqrt(30.0))) &
              & * real(-1.0 * ii * (rhow * deltaS(2,1,nfiles) + rhoi * deltaicestar(2,1,nfiles)))    

      dil(:,:,nfiles) = il(:,:) - oldil(:,:) 

      ! Calculate m (rotation vector perturbation) 
      !  - From Jerry's written notes, also Mitrovica, Wahr, Matsuyama, and Paulson (2005) eq. 2
      sum_il(:,:) = 0.0
      sum_m(:) = 0.0
      do nn = 1,nfiles-1 ! Sum over all previous timesteps
         betatt = 0.0
         betattprime = 0.0
         do k = 1,nmod(2) ! Sum over k=1,K
            betatt = betatt + (rprime(k,2) / s(k,2)) &
                     & * ( 1.0 - exp(-1.0 * s(k,2) * (times(nfiles) - times(nn)) / 1000.0) )
            betattprime = betattprime + &
                        & (rprimeT(k,2) / s(k,2)) * ( 1.0 - exp(-1.0 * s(k,2) * (times(nfiles) - times(nn)) / 1000.0) )
         enddo
         sum_il(:,:) = sum_il(:,:) + dil(:,:,nn) * betatt
         sum_m(:) = sum_m(:) + dm(:,nn) * betattprime
      enddo
    
      do i = 1,2
         mm(i) = (1 / (1 - kTE(2) / kf)) * &
               & ( (1 / (moiC - moiA)) * ((1 + kE(2)) * il(i,3) + sum_il(i,3)) + (1 / kf) * sum_m(i) )
      enddo
      mm(3) = (1.0 / moiC) * ( (1 + kE(2)) * il(3,3) + sum_il(3,3) )
       
      dm(:,nfiles) = mm(:) - oldm(:)

      ! Calculate lambda (rotational driving) from m (Milne and Mitrovica 1998)
      lambda(0,0) = (radius**2 * omega**2 / 3.0) * ((mm(1)**2 + mm(2)**2 + mm(3)**2) + 2.0 * mm(3))
      lambda(2,0) = (radius**2 * omega**2 / (6.0 * sqrt(5.0))) &
                  & * (mm(1)**2 + mm(2)**2 - 2.0 * mm(3)**2 - 4.0 * mm(3))
      lambda(2,1) = (radius**2 * omega**2 / sqrt(30.0)) * ((mm(1) - ii * mm(2)) * (1.0 + mm(3)))
      lambda(2,2) = (radius**2 * omega**2 / sqrt(5.0 * 24.0)) * (mm(2)**2 - mm(1)**2 + 2.0 * ii * mm(1) * mm(2))
       
      dlambda(:,:,nfiles) = lambda(:,:) - deltalambda(:,:,nfiles-1)
     
      ! Calculate effect on sea level (Love numbers) (Kendall)
      dsl_rot(:,:) = (0.0,0.0)
      ekhTE = 1 + kTE(2) - hTE(2)
      do nn = 1,nfiles-1 ! Sum over all previous timesteps
         lovebetatt(nn) = 0.0
         do k = 1,nmod(2) ! Sum over k=1,K
            lovebetatt(nn) = lovebetatt(nn) + ((rprimeT(k,2) - rT(k,2)) / s(k,2)) & 
                             * ( 1 - exp(-1.0 * s(k,2) * (times(nfiles) - times(nn)) / 1000.0) ) ! (eq. B27)
         enddo
      enddo
     
      do m = 0,2
         viscoustt = (0.0,0.0)
         do nn = 1,nfiles-1 ! Sum the loads over all previous timesteps to get the sum on the 4th line of eq. B28
            viscoustt = viscoustt + lovebetatt(nn)*dlambda(2,m,nn)
         enddo       
         dsl_rot(2,m) = (ekhTE * (deltalambda(2,m,nfiles-1) + dlambda(2,m,nfiles)) + viscoustt) / gacc ! (eq. B28/B25)
      enddo         
        
      if (calcRG) then ! For R calculations
         rr_rot(:,:) = (0.0,0.0)
         do nn = 1,nfiles-1 ! Sum over all previous timesteps
            lovebetattrr(nn) = 0.0
            do k = 1,nmod(2) ! Sum over k=1,K
               lovebetattrr(nn) = lovebetattrr(nn) & 
                                 + ((rT(k,2)) / s(k,2)) & 
                                 * (1 - exp(-1.0 * s(k,2) * (times(nfiles) - times(nn)) / 1000.0)) ! (eq. B27)
            enddo
         enddo
         do m = 0,2
            viscousttrr = (0.0,0.0)
            do nn = 1,nfiles-1 ! Sum the loads over all previous timesteps to get the sum on the 4th line of eq. B28
               viscousttrr = viscousttrr + lovebetattrr(nn) * dlambda(2,m,nn)
            enddo
            rr_rot(2,m) = (hTE(2) * (deltalambda(2,m,nfiles-1) + dlambda(2,m,nfiles)) + viscousttrr) / gacc
         enddo
      endif
       
   else ! If not incorporating rotation
      dsl_rot(:,:) = (0.0,0.0)
      rr_rot(:,:) = (0.0,0.0)
   endif ! TPW
          
   !-----/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\    End Rotation    /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\----!
   
   ! Calculate beta (eq. B14) - viscous response factor
   do l = 1,norder
      do nn = 1,nfiles-1 ! Sum over all previous timesteps
         lovebeta(nn,l) = 0.0
         do k = 1,nmod(l) ! Sum over k=1,K
            lovebeta(nn,l) = lovebeta(nn,l) + ((rprime(k,l) - r(k,l)) / s(k,l)) & 
                            * (1 - exp(-1.0 * s(k,l) * (times(nfiles) - times(nn)) / 1000.0)) ! (eq. B14)
         enddo
      enddo
   enddo
   
   if (calcRG) then ! For R calculations
      do l = 1,norder
         do nn = 1,nfiles-1 ! Sum over all previous timesteps
            lovebetarr(nn,l) = 0.0
            do k = 1,nmod(l) ! Sum over k=1,K (modes)
               lovebetarr(nn,l) = lovebetarr(nn,l) & 
                              & + ((r(k,l)) / s(k,l)) * (1 - exp(-1.0 * s(k,l) * (times(nfiles) - times(nn)) / 1000.0))
                              ! (eq. B14)
            enddo
         enddo
      enddo
   endif
   
   ! Compute viscous response outside inner loop, since it doesn't depend on current timestep
   viscous(:,:) = (0.0,0.0)
   do l = 1,norder
      do m = 0,l
         do nn = 1,nfiles-1 ! Sum the loads over all previous timesteps to get the sum on the second line of eq. B18
            viscous(l,m) = viscous(l,m) + lovebeta(nn,l) * (rhoi * dicestar(l,m,nn) + rhow * dS(l,m,nn))
         enddo
      enddo
   enddo

   ! Calculate change in sea level from change in load (STEP 6)  
   dsllm(:,:) = (0.0,0.0) ! No change on first timestep
   do l = 1,norder
      ttl = (4.0 * pi * (radius**3)) / (mass * (2 * real(l) + 1.0)) ! (eq. B16)
      ekhl = 1 + ke(l) - he(l) ! (eq. B13)
      do m = 0,l      
         ! Total  = _____________________________elastic_____________________________ + _viscous_
         dsllm(l,m) = ttl * ekhl * (rhoi * deltaicestar(l,m,nfiles) + rhow * deltaS(l,m,nfiles-1) + rhow * dSlm(l,m)) & 
                       + ttl * viscous(l,m) ! (eq. B18)
      enddo
   enddo

   ! Add rotational effects (calculated above)
   dsllm(0:2,0:2) = dsllm(0:2,0:2) + dsl_rot(0:2,0:2)
   
   ! Convert dSL (total spatially heterogeneous change since time0) to spatial domain
   call spec2spat(dslxy, dsllm, spheredat)
   
   ! Compute r0 (STEP 7)
   rOxy(:,:) = dslxy(:,:) * cstarxy(:,:) ! (eq. 68)
   call spat2spec(rOxy, rOlm, spheredat) ! Decompose r0

   ! Compute conservation term (STEP 8)
   conserv = (1 / cstarlm(0,0)) * (-1.0 * (rhoi / rhow) * deltaicestar(0,0,nfiles) - rOlm(0,0) + tOlm(0,0)) ! (eq. 78)

   ! Compute change in ocean load again (STEP 9)
   olddSlm = dSlm ! Save previous-iterate load to check convergence
   
   dSlm(:,:) = -1.0 * deltaS(:,:,nfiles-1) + rOlm(:,:) + conserv * cstarlm(:,:) - tOlm(:,:) ! (eq. 73)
   deltaS(:,:,nfiles) = deltaS(:,:,nfiles-1) + dSlm(:,:)! Save total load changes, used in Love number calculation

    ! Update the total sea-level change up to the current time step
   deltasllm(:,:) = dsllm(:,:) ! Spatially heterogeneous component
   deltasllm(0,0) = deltasllm(0,0) + conserv ! Add uniform conservation term to (0,0)
   call spec2spat(deltaslxy, deltasllm, spheredat) ! Synthesize deltasl
 
   ! Calculate convergence criterion for inner loop
   if ( abs(sum(abs(dSlm)) - sum(abs(olddSlm))) < epsilon(0.0) .and. abs(sum(abs(olddSlm))) < epsilon(0.0)) then
       xi = 0 ! Otherwise xi = 0 / 0 = NaN at the first loop of the first timestep.
   elseif (abs(sum(abs(olddSlm))) < epsilon(0.0)) then
       xi = abs( (sum(abs(dSlm)) - sum(abs(olddSlm))) / (epsilon(0.0) * 10) ) ! Avoid dividing by 0
   else
       xi = abs( (sum(abs(dSlm)) - sum(abs(olddSlm))) / sum(abs(olddSlm)) ) ! (eq. 83)
   endif
      
   ! If the ocean loading guess has not been converged, 
   if (xi > epsilon1) then      
       ! new guess to the topography correction
       topoxy(:,:) = tinit(:,:) - deltaslxy(:,:) ! (eq. 39)
        
       ! new ocean function
       do j = 1,2*nglv
          do i = 1,nglv
             if (topoxy(i,j) >= 0.0) then
                cxy(i,j) = 0.0
             else
                cxy(i,j) = 1.0
             endif
          enddo
       enddo
       
       !new guess to ocean*beta function
       cstarxy(:,:) = cxy(:,:) * beta(:,:) ! (eq. 65)
       call spat2spec(cstarxy, cstarlm, spheredat) ! Decompose cstar
       
      ! new guess to the topo correction
      tOxy(:,:) = tinit * (cstarxy(:,:) - cstar0(:,:)) ! (eq. 70)
      call spat2spec(tOxy, tOlm, spheredat) ! Decompose tO   
   endif
   
   if (xi <= epsilon1) then ! If converged
      exit
   elseif (ninner == 9999) then ! If no convergence after a huge number of iterations
      write(*,*)
      write(*,'(A,I5,A)') 'WARNING: The inner loop failed to converge after the limit of ', ninner, ' iterations.'
      write(*,'(A,ES15.3,A)') '         The variable xi finished with a value of ', xi, ', resulting from '
      write(*,'(A,ES15.3)') '         sum(abs(olddSlm)) = ', sum(abs(olddSlm))
      write(*,'(A,ES15.3)') '            sum(abs(dSlm)) = ', sum(abs(dSlm))
      write(*,*)
      write(*,'(A)') '!!!---- Program sl_model will now be terminated. ----!!!'
      call abort ! Terminate program
      !exit ! DEBUG line: Continue inner loop despite non-convergence. Normal operation: Enable the 2 lines above.
   endif
   ninner = ninner + 1
      
enddo ! End inner loop
!-----------------------------------------------------------
!<<<<<<<<<<<<<<<<<<<< End of inner loop >>>>>>>>>>>>>>>>>>>>
!-----------------------------------------------------------
write(*,'(A,I4,A)') '  ', ninner, ' inner-loop iterations'

!HH: print out the number of iteration it takes for the inner convergence
open(unit = 1, file = outputfolder//'numiter'//ext, form = 'formatted', access ='sequential', &
& status = 'old',position='append')
write(1,'(I5)') ninner
close(1)

! Write out the converged rotation-related quantities 
if (tpw) then
   open(unit = 1, file = outputFolder//'TPW'//ext, &
   & form = 'formatted', access = 'sequential', status = 'old', position='append')
!        write(1,'(9ES19.8E2/,3ES19.8E2/,18ES19.8E2)') dil(:,:,nfiles), dm(:,nfiles), dlambda(:,:,nfiles)
   write(1,'(9ES19.8E2/,3ES19.8E2/,18ES19.8E2)') il(:,:), mm(:), lambda(:,:)
   close(1)
endif
 
!      write(*,*) 'dil', dil(:,:,nfiles)
!      write(*,*) 'dm', dm(:,nfiles)
!      write(*,*) 'dlambda', dlambda(:,:,nfiles)
 
if (calcRG) then ! For R calculations
   if (nmelt == 0) then
      rrlm(:,:) = (0.0,0.0)! No change on first timestep
   else
      do l = 1,norder
         ttl = (4.0 * pi * (radius**3)) / (mass * (2.0 * real(l) + 1.0)) ! (eq. B16)
         do m = 0,l      
            ! Viscous response
            viscousrr = (0.0,0.0)
            do nn = 1,nfiles-1 ! Sum the loads over all previous timesteps to get the sum on the second line of eq. B18
               viscousrr = viscousrr + lovebetarr(nn,l) * (rhoi * dicestar(l,m,nn) + rhow * dS(l,m,nn))
            enddo
            rrlm(l,m) = ttl * he(l) * (rhoi * deltaicestar(l,m,nfiles) + rhow * deltaS(l,m,nfiles-1) + rhow * dSlm(l,m)) &
                     & + ttl * viscousrr
         enddo
      enddo
   endif
   rrlm(0:2,0:2) = rrlm(0:2,0:2) + rr_rot(0:2,0:2)
   call spec2spat(rrxy, rrlm, spheredat)
   rr(:,:,n) = rrxy(:,:)
endif
 

! Update topography fieldS, ocean functions (STEP 11)
oldt0lm = t0lm ! Save old initial topography to check convergence

! Update the topography at the current time step 
topoxy(:,:) = tinit(:,:) - deltaslxy(:,:) ! (eq. 12)

! Update ocean function
do j = 1,2*nglv
   do i = 1,nglv
      if (topoxy(i,j) >= 0.0) then
         cxy(i,j) = 0.0
      else
         cxy(i,j) = 1.0
      endif
   enddo
enddo

!=========================================================================================
!                          OUTPUT                           
!_________________________________________________________________________________________

if (nmelt.GT.0) then

   write(*,'(A)') 'Writing output files...'
   j = TIMEWINDOW(nfiles)
   write(*,*) 'FILENUMBER of new outputs : ',j
   write(numstr,'(I4)') j
   numstr = trim(adjustl(numstr))
  
   !Total sea level change from the beginning of the simulation
!   open(unit = 1, file = outputfolder//'SL'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
!   & status = 'replace')
!   write(1,'(ES16.9E2)') tinit_0(:,:)-topoxy(:,:)
!   close(1)

   !HH: print out the nmelt
   open(unit = 1, file = outputfolder//'nmelt'//ext, form = 'formatted', access ='sequential', &
   & status = 'old',position='append')
   write(1,'(I4)') nmelt
   close(1)

   ! topography at the current timestep
   open(unit = 1, file = outputfolder//'tgrid'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
   & status = 'replace')
   write(1,'(ES16.9E2)') topoxy(:,:)
   close(1)

   ! converged ocean function at the current timestep
   open(unit = 1, file = outputfolder//'ocean'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
   & status = 'replace')
   write(1,'(ES14.4E2)') cxy(:,:)
   close(1)
   
   ! converged beta function at the current timestpe
   open(unit = 1, file = outputfolder//'beta'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
   & status = 'replace')
   write(1,'(ES14.4E2)') beta(:,:)
   close(1)   

   ! output converged total ocean loading changes 
   open(unit = 1, file = outputfolder//'dS_converged'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
   & status = 'replace')
   write(1,'(ES14.7E2)') deltaS(:,:,nfiles)
   close(1)
   
   if (iceVolume) then
	  ice_volume = icestarlm(0,0)*4*pi*radius**2 !multiply the (0,0) component of ice to the area of a sphere
      open(unit = 1, file = outputfolder//'ice_volume'//ext, form = 'formatted', access = 'sequential', &
      & status = 'old', position = 'append')
      write(1,'(ES14.4E2)') ice_volume
      close(1)
   endif
   
   current_time = iter * dt1     !time passed since the start of the simulation  
   if (current_time == L_sim) then !if we are at the last time step of simulation
     write(iterstr,'(I2)') itersl
     iterstr = trim(adjustl(iterstr))
     write(*,*) 'Last time step of the simulation! writing out files for next outer-iteration loop'
     ! write out the predicted present day topography into a file so it can be used in the next outer-iteration

     open(unit = 1, file = outputfolder//'pred_pres_topo_'//trim(iterstr)//ext, form = 'formatted',  &
     & access = 'sequential', status = 'replace')
     write(1,'(ES16.9E2)') topoxy(:,:)
     close(1)
     
     ! write out the initial topography of the simulation at the currect outer-loop into a file
     open(unit = 1, file = outputfolder//'tgrid0_'//trim(iterstr)//ext, form = 'formatted',  &
     & access = 'sequential', status = 'replace')
     write(1,'(ES16.9E2)') tinit_0(:,:)
     close(1)
   endif 
   
   if (calcRG) then
      open(unit = 1, file = outputfolder//'R'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
      & status = 'replace')
      write(1,'(ES14.4E2)') rr(:,:,n)
      close(1)
      
      ! Compute geoid displacement
      gg(:,:,n) = deltaslxy(:,:)+rr(:,:,n)
    
      open(unit = 1, file = outputfolder//'G'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
      & status = 'replace')
      write(1,'(ES14.4E2)') gg(:,:,n)
      close(1)
   endif

   !   open(unit = 1, file = outputfolder//'Tice'//numstr, form = 'formatted', access = 'sequential', &
   !   & status = 'replace')
   !   write(1,'(ES14.4E2)') T(:,:,n) + icexy(:,:,n)
   !   close(1)

   !   open(unit = 1, file = outputfolder//'Tbr'//numstr, form = 'formatted', access = 'sequential', &
   !   & status = 'replace')
   !   write(1,'(ES14.4E2)') T(:,:,n)
   !   close(1)
   
   if (coupling) then 
      ! topography change between the previous and the current timestep 
      ! this is the information passed to the ice sheet model

      open(unit = 1, file = folder_coupled//'bedrock'//ext, form = 'formatted', access = 'sequential', &
      & status = 'replace')
      write(1,'(ES16.9E2)') topoxy_m1(:,:)-topoxy(:,:)
      close(1)
      
      !write out the current ice load as a new file
      open(unit = 1, file = outputfolder_ice//icemodel_out//trim(numstr)//ext, form ='formatted', access = 'sequential', &
      & status = 'replace')
      write(1,'(ES16.9E2)') icexy(:,:,nfiles)
      close(1)
   endif !endif coupling

endif !endif nmaelt>0

call system_clock(countf) ! Total time
call cpu_time(countf_cpu)
if (nmelt .GT. 0) then 
   ! Write out total compuatation time of sea level change over current timestep
   open(unit = 1, file = outputfolder//'elapsed_wall_time'//ext, form = 'formatted', access = 'sequential', &
   & status = 'old', position='append')
   write(1,'(ES14.4E2)') float(countf-counti)/float(countrate)
   close(1)

   open(unit = 1, file = outputfolder//'elapsed_cpu_time'//ext, form = 'formatted', access = 'sequential', &
   & status = 'old', position='append')
   write(1,'(ES14.4E2)') countf_cpu-counti_cpu
   close(1)
endif

write(*,'(A,F7.2,A)') 'Done! Total time ', real(countf - counti) / real(countrate), ' seconds'
write(*,*) ''
write(*,*) ''

if (Travel_total > 0 .and. Travel == Travel_total) then
   write(*,*) ' THE TIMEWINDOW HAS REACHED THE END OF THE SIMULATION!'
   write(*,*) ' GREAT JOB TW!'
endif
write(*,*) ''
end program sl_model
