!=======================================================================================================MAIN MODULE=====!
!_______________________________________________________________________________________________________________________!
module sl_model_mod

   use spharmt
   use user_specs_mod
   use planets_mod
   use sl_io_mod

   implicit none
   private

   public :: sl_timewindow, set_planet, sl_solver_checkpoint, sl_solver_init, &
             sl_solver, sl_deallocate_array, sl_set_unit_num , sl_call_readnl! public module
   public :: iterstr ! public variables


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

   ! grid lat-lon
   real, dimension(nglv)        :: latgrid
   real, dimension(2*nglv)      :: longrid

   ! Model calculations
   real, dimension(nglv,2*nglv) :: glw_matrix               ! weigtht in the Gaussian-Legendre grid
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
   real :: ice_volume, mean_slc                             ! ice volume and mean sea level change
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
   integer :: counti, countf,countrate                         ! Computation timing
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
   integer :: nmelt, nfiles                                ! Number of melting episodes up to the current timestep
!   real :: starttime                                       ! start time of the simulation
!   integer :: nmelt, nfiles, iter, itersl, dtime           ! Number of melting episodes up to the current timestep
   integer :: iargc, nargs                                 ! Arguments read in from a bash script
   character(16) :: carg(20)                               ! Arguments from a bash script
   character(3) :: skip                                    ! variable used to skip lines in reading TPW file

   !HH variables for debugging
   real, dimension(nglv,2*nglv) :: deltaslxy_outside_maliMesh, deltaslxy_inside_maliMesh
   complex, dimension(0:norder,0:norder) :: deltasllm_outside_maliMesh, deltasllm_inside_maliMesh
   complex, dimension(0:norder,0:norder) :: maliMask_lm00, kFactor2_from_mali_lm00
   real, dimension(nglv,2*nglv) :: mask_maliMesh
   real :: mean_slc_outside_maliMesh, mean_slc_inside_maliMesh, mean_slc_sum, area_maliMesh
   contains

   !======================================================================================================================!
   subroutine set_planet
      ! Planetary values
      if (whichplanet == 'earth' .or. whichplanet == 'Earth' .or. whichplanet == 'EARTH') then
         call earth_init
      elseif (whichplanet == 'mars' .or. whichplanet == 'Mars' .or. whichplanet == 'MARS') then
         call mars_init
      else
         write(unit_num,*) 'The parameters for the planet you entered are not built in.'
         write(unit_num,*) 'Please check your spelling for the variable whichplanet in the user_specs_mod module.'
         write(unit_num,*) 'If you preferred, you could create a new subroutine for the new planet in the planets_mod module.'
         write(unit_num,*) 'Terminating: program sl_model'
         stop
      endif
   end subroutine set_planet
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================!
   subroutine sl_call_readnl

      call sl_readnl(inputfolder_ice, inputfolder, planetfolder, gridfolder, &
                        outputfolder, outputfolder_ice, folder_coupled, ext, fType_in, fType_out, &
                        planetmodel, icemodel, icemodel_out, timearray, &
                        topomodel, topo_initial, grid_lat, grid_lon, &
                        checkmarine, tpw, calcRG, input_times, &
                        initial_topo, iceVolume, coupling, patch_ice, &
                        L_sim, dt1, dt2, dt3, dt4, Ldt1, Ldt2, Ldt3, Ldt4, whichplanet)

   end subroutine sl_call_readnl
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================!
   subroutine sl_timewindow(iter) !create a time window that points to which ice history to read in

      integer :: iter

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

      if (sum(Ldt) > L_sim) then
          write(unit_num,*) 'Total length of simulation CANNOT be smaller than the total lengths of the internal time windows!!'
          write(unit_num,*) 'Make sure the sum(LdT) ≤ L_sim !!'
          write(unit_num,*) 'Terminating: program sl_model'
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
      Travel_total = (L_sim - L_TW)/dt1 !total number of marching steps taken by the TW

      if (Travel_total == 0) then
         write(unit_num,*) 'The length of a total simulation and the length of timewindow are the same.'
         write(unit_num,*)' TW will not march forward.'
      endif

      write(unit_num,*) 'Total length of time window in years (L_TW)', L_TW
      write(unit_num,*) 'Total number of melting episode within a Full TW (TW_nmelt)', TW_nmelt
      write(unit_num,*) 'Total number of icefiles in a FULL TW (TW_nfiles)', TW_nfiles

      if (iter.LT.ncalls) then !while the time window is growing, the time window has not started travelling when
         Travel = 0
         if (iter == 0) then
            nfiles = 1
         endif
      else
         Travel = iter - ncalls ! This is the number of marching steps the full TW has taken upto the current timestep
      endif

      ! allocate dimensions to arrays
      allocate (mask(ncalls+1),iceload(ncalls+1),icefiles(ncalls+1))
      mask(:) = 0
      iceload(:) = 0
      icefiles(:) = 0
      allocate (TIMEWINDOW(TW_nfiles))
      TIMEWINDOW(:)=0

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
         write(unit_num,*) 'Icefiles in the TW at current marching step:', TIMEWINDOW
         write(unit_num,*) 'The TW will stop marching when TRAVEL == Travel_total:'
         write(unit_num,*) ''
         write(unit_num,'(A,I4,A)') '   ', Travel_total - Travel, ' more TW steps to march!'
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

      times = 0.0
      lovebetatt = 0.0
      lovebetattrr = 0.0
      lovebetarr = 0.0
      lovebeta = 0.0
      icexy = 0.0
      sl = 0.0
      dS = 0.0
      deltaS = 0.0
      dicestar = 0.0
      deltaicestar = 0.0
      rr = 0.0
      gg = 0.0
      dil = 0.0
      dlambda = 0.0
      deltalambda = 0.0
      dm = 0.0

      !!!!!!!!!!!!!!!!!!!ATTENTION!!!!!!!!
      ! TIMEWINDOW = TIMEWINDOW+1
      ! TIMEWINDOW = -(TIMEWINDOW*100+starttime)
      ! CAUTION in modifying the ice file numbers

      !write(unit_num,*) 'Timewindow is growing, ice file numbers to read in :', TIMEWINDOW(1:nfiles)
      write(unit_num,*) 'number of melting episodes (nmelt):', nmelt
      write(unit_num,'(A,I4)') 'Total number of files in the current time window (nfiles) = ', nfiles
      write(unit_num,'(A,I4)') ' Number of marching steps the TW has taken =        ', Travel
      write(unit_num,*) ''

   end subroutine sl_timewindow
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================!
   subroutine sl_deallocate_array

      deallocate(mask, iceload, icefiles, TIMEWINDOW)
      deallocate(times, lovebetatt, lovebetattrr)
      deallocate(lovebetarr, lovebeta)
      deallocate(icexy, sl)
      deallocate(dS, deltaS)
      deallocate(dicestar, deltaicestar)
      deallocate(rr, gg)
      deallocate(dil, dlambda ,deltalambda)
      deallocate(dm)

   end subroutine sl_deallocate_array
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================!
   subroutine sl_set_unit_num(unit_num_in)

      integer :: unit_num_in

         unit_num = unit_num_in

   end subroutine sl_set_unit_num
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================!
   subroutine sl_solver_checkpoint(itersl, dtime)

      integer :: itersl, dtime

      if (dtime /= dt1) then
         write(unit_num,*) 'dtime and dt1 should be equal to each other.'
         write(unit_num,*) 'Please check your set up for the variables'
         write(unit_num,*) 'Terminating: program sl_model'
         stop
      endif

      if (itersl.lt.1) then
          write(unit_num,*) 'itersl must be equal to or greather than 1'
          write(unit_num,*) 'itersl = 1: No topography correction'
          write(unit_num,*) 'itersl > 1: topography correction'
          write(unit_num,*) 'Terminating: program sl_model'
          stop
      endif

      if (coupling) then ! Set a default unit number for I/O
         write(unit_num,*) 'Sea level model is coupled to the ice sheet model. Unit number is:', unit_num
      else
         call sl_set_unit_num(6)! Set a default unit number for I/O
         write(unit_num,*) 'Sea level model is running standalone. Unit number is:', unit_num
      endif

      ! set up the planet profile
      call set_planet

      ! initalize spherical harmonics
      call spharmt_init(spheredat, 2*nglv, nglv, norder, radius) ! Initialize spheredat (for SH transform subroutines)

   end subroutine sl_solver_checkpoint
   !_______________________________________________________________________________________________________________________!


   !=======================================================================================================================
   subroutine sl_solver_init(itersl, starttime, mali_iceload, mali_bedrock, mali_mask)

      integer :: itersl
      real :: starttime
      real, dimension(nglv,2*nglv), optional :: mali_iceload, mali_bedrock, mali_mask

      !================================================================================================
      !       Initialize the model at times(1) when there has been no melting episodes yet  NMELT = 0
      !================================================================================================

      write(unit_num,*) 'nmelt=0, INITIALIZING THE SEA LEVEL MODEL..'
      flush(unit_num)

      !initialize variables
      j = TIMEWINDOW(1) ! initial file number (i.e. 0)
      write(unit_num,'(A,I4)') 'initial file number:', j
      write(numstr,'(I4)') j
      numstr = trim(adjustl(numstr))
      !     k = -1*starttime
      !     write(numstr2,'(I6)') k
      !     write(unit_num,'(A,I6)') 'initial icefile suffix:', k
      !     numstr2 = trim(adjustl(numstr2))

      !====================== topography and ice load========================
      ! read in the initial iceload from the coupled ice input folder
      call read_sl(icexy(:,:,1), icemodel, inputfolder_ice, suffix=numstr)

      !  Initialize topography (STEP 1)
      if (initial_topo) then
         write(unit_num,*) 'Reading in initial topo file'
         call read_sl(tinit_0, topo_initial, inputfolder)
      else  ! if initial topo is unknown

         ! Present-day observed topography
         write(unit_num,*) 'Reading in ETOPO2 file'
         call read_sl(truetopo, topomodel, inputfolder)

         ! if topography is not iteratively getting improved, or if at the first loop of outer-iteration
         if (itersl == 1) then
            write(unit_num,*) 'Initial topo is unknown, set the modern-observed topo to be the initial topo'
            ! truetopo(:,:)=truetopo(:,:)-icexy(:,:,TW_nfiles) ! If using ice topography
            tinit_0(:,:) = truetopo(:,:)! (eq. 48)
          else

             write(unit_num,*) 'Topographic correction is ON. Updating initial topography'
             ! read in predicted present topo from the previous outer-iteration 'itersl-1'
             write(iterstr,'(I2)') itersl-1
             iterstr = trim(adjustl(iterstr))

             call read_sl(pred_pres_topo, 'pred_pres_topo_', outputfolder, suffix=iterstr)

             ! read in tinit_0 from the previous outer-iteration 'itersl-1'
             call read_sl(tinit_0_last, 'tgrid0_', outputfolder, suffix=iterstr)

             ! compute topography correction for the initial topography
             init_topo_corr(:,:) = truetopo(:,:) - pred_pres_topo(:,:)
             tinit_0(:,:) = tinit_0_last(:,:) + init_topo_corr(:,:)
          endif
       endif

       if (coupling) then  !if coupling the ICE SHEET - SEA LEVEL MODELs
          call spat2spec(mali_mask(:,:),maliMask_lm00(:,:),spheredat)
          area_maliMesh = maliMask_lm00(0,0)*4*pi*radius**2
          open(unit = 201, file = trim(outputfolder)//'area_mali_domain', form = 'formatted', access ='sequential', &
          & status = 'replace')
          write(201,'(ES14.4E2)') area_maliMesh
          close(201)

          write(unit_num,*) 'Merge initial topography and iceload with the ISM data'

          ! merge intitial topography with bedrock provided by the ice sheet model.
          do j = 1,2*nglv
             do i = 1,nglv
                tinit_0(i,j) = mali_bedrock(i,j) + tinit_0(i,j)*(1 - mali_mask(i,j))
             enddo
          enddo

          ! merge intitial topography with bedrock provided by the ice sheet model.
          if (patch_ice) then
             do j = 1,2*nglv
                do i = 1,nglv
                   icexy(i,j,nfiles) = 0.0
                enddo
             enddo
          endif

          do j = 1,2*nglv
             do i = 1,nglv
                icexy(i,j,nfiles) = mali_iceload(i,j) + icexy(i,j,nfiles)*(1 - mali_mask(i,j))
             enddo
          enddo

          !write out the current ice load as a new file to the sea-level model ice folder
          call write_sl(icexy(:,:,nfiles), icemodel_out, outputfolder_ice, suffix=numstr)
          !HH write out the mali mask file for debugging
          call write_sl(mali_mask(:,:), "mali_mask_init", outputfolder)
      endif ! end if (coupling)

      !write out the initial topo of the simulation, tinit_0
      call write_sl(tinit_0(:,:), 'tgrid', outputfolder, suffix=numstr)


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
      call write_sl(cxy0(:,:), 'ocean', outputfolder, suffix=numstr)


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
      call write_sl(beta0(:,:), 'beta', outputfolder, suffix=numstr)


      !================== total ocean loading change =====================
      ! initialize the total ocean loading change and output as a file
      deltaS(:,:,1) = (0.0,0.0)
      open(unit = 201, file = trim(outputfolder)//'dS_converged'//trim(numstr), form = 'formatted', access = 'sequential', &
      & status = 'replace')
      write(201,'(ES16.9E2)') deltaS(:,:,1)
      close(201)


      !========================== computing time =========================
      ! To write out how much time it took to compute sea-level change over one step
      ! Open a new file
      open(unit = 201, file = trim(outputfolder)//'elapsed_wall_time', form = 'formatted', access = 'sequential', &
      & status = 'replace')
      close(201)

      open(unit = 201, file = trim(outputfolder)//'elapsed_cpu_time', form = 'formatted', access = 'sequential', &
      & status = 'replace')
      close(201)


      !========================== time array =============================
      if (.not. input_times) then !if time array is not read in from a text file, make a new one
         ! write a new file
         open(unit = 201, file = trim(outputfolder)//trim(timearray), form = 'formatted', access = 'sequential', &
         & status = 'replace')
         write(201,'(ES14.4E2)') starttime
         close(201)
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
      open(unit = 201, file = trim(outputfolder)//'TPW', form = 'formatted', access = 'sequential', &
      & status = 'replace')
      write(201,'(9ES19.8E2/,3ES19.8E2/,18ES19.8E2)') il(:,:), mm(:), lambda(:,:)
      ! write(unit_num,'(9ES19.8E2/,3ES19.8E2/,18ES19.8E2)') dil(:,:,1), dm(:,1), dlambda(:,:,1)
      close(201)
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

         open(unit = 201, file = trim(outputfolder)//'ice_volume', form = 'formatted', access ='sequential', &
         & status = 'replace')
         write(201,'(ES14.4E2)') ice_volume
         close(201)

         mean_slc = 0.0
         ! for mean sea level directly calculated based on the deltasl variable
         open(unit = 201, file = trim(outputfolder)//'gmsle_change', form = 'formatted', access ='sequential', &
         & status = 'replace')
         write(201,'(ES14.4E2)') mean_slc
         close(201)

         open(unit = 201, file = trim(outputfolder)//'gmsle_change_outside_maliMesh', form = 'formatted', access ='sequential', &
         & status = 'replace')
         write(201,'(ES14.4E2)') mean_slc
         close(201)

         open(unit = 201, file = trim(outputfolder)//'gmsle_change_inside_maliMesh', form = 'formatted', access ='sequential', &
         & status = 'replace')
         write(201,'(ES14.4E2)') mean_slc
         close(201)

         open(unit = 201, file = trim(outputfolder)//'gmsle_change_sumMesh', form = 'formatted', access ='sequential', &
         & status = 'replace')
         write(201,'(ES14.4E2)') mean_slc
         close(201)
      endif

      !print out the number of iteration it takes for the inner convergence
      open(unit = 201, file = trim(outputfolder)//'numiter', form = 'formatted', access ='sequential', &
      & status = 'replace')
      close(201)

      !print out the number of melting episodes taken into account at the current timestep (i.e., nfiles-1)
      open(unit = 201, file = trim(outputfolder)//'nmelt', form = 'formatted', access ='sequential', &
      & status = 'replace')
      close(201)

      write(unit_num,*) 'DONE INITIALIZATION.'
   !   call exit
   end subroutine sl_solver_init

   !========================================================================================================================
   subroutine sl_solver(itersl, iter, dtime, starttime, mali_iceload, mali_mask, kFactor2_from_mali, slchange, kFactor2Inv_to_mali, mali_mask_ones_only)
      ! Compute sea-level change associated with past ice loading changes

      real :: starttime, k0, lat0, area_maliMesh_kFactor2
      integer :: iter, itersl, dtime
      real, dimension(nglv) :: lat512
      real, dimension(2*nglv) :: lon512
      real, dimension(nglv,2*nglv) :: kFactor
      real, dimension(2*nglv,nglv) :: kFactor_test
      real, dimension(nglv,2*nglv), optional :: mali_iceload, mali_mask, kFactor2_from_mali ! variables for coupled ISM-SLM simulations
      real, dimension(nglv,2*nglv), intent(out), optional :: slchange, kFactor2Inv_to_mali ! variable exchanged with the ISM
      integer, dimension(nglv,2*nglv), intent(out), optional :: mali_mask_ones_only

      !===========================================================
      !                   BEGIN TIMING & EXECUTION
      !___________________________________________________________
      call cpu_time(counti_cpu)
      call system_clock(count = counti, count_rate = countrate)  ! Total computation time

      if (coupling) then
         ! calculate area of the mali domain based on the interpolated kFactor^2
         call spat2spec(kFactor2_from_mali(:,:),kFactor2_from_mali_lm00(:,:),spheredat)
         area_maliMesh_kFactor2 = kFactor2_from_mali_lm00(0,0)*4*pi*radius**2
         open(unit = 201, file = trim(outputfolder)//'area_mali_domain_kFactor2', form = 'formatted', access ='sequential', &
         & status = 'replace')
         write(201,'(ES14.4E2)') area_maliMesh_kFactor2
         close(201)

         ! create mali masks of 1 on the sea-level model
         mali_mask_ones_only(:,:) = 0.0
         do j = 1,2*nglv
            do i = 1,nglv
                if (mali_mask(i,j) .NE. 0) then
                   mali_mask_ones_only(i,j) = 1
                endif
            enddo
         enddo

         open(unit = 201, file = trim(gridfolder)//trim(grid_lat), form = 'formatted', &
         & access = 'sequential', status = 'old')
         read(201,*) lat512
         close(201)
         open(unit = 201, file = trim(gridfolder)//trim(grid_lon), form = 'formatted', &
         & access = 'sequential', status = 'old')
         read(201,*) lon512
         close(201)

         do j = 1,2*nglv
            do i = 1,nglv
                if (mali_mask(i,j) .NE. 0) then
                   mali_mask_ones_only(i,j) = 1
                endif
            enddo
         enddo

         open(unit = 201, file = trim(outputfolder)//trim('mali_mask_ones_only_on_GLgrid'), form = 'formatted', &
             & access = 'sequential', status = 'replace')
             write(201,'(I2)') mali_mask_ones_only
         close(201)

         k0 = (1 - sin(-71 * 3.14159265359 / 180)) / 2 ! in radian
         lat0 = -90 * 3.14159265359 / 180 ! center of the latitude in radian
         !lon0 = 0 ! center of the longitude
         kFactor= 0.0
         do j = 1,2*nglv
            do i = 1,nglv
                !kFactor(i,j) =  2 * k0 / ( 1 + (sin(lat0) * sin(lat512(i) * 3.14159265359 / 180)) + (cos(lat0) * cos(lat512(i)* 3.14159265359 / 180) * cos((lon512(j)-180)* 3.14159265359 / 180)) )
                kFactor_test(j,i) =  2 * k0 / ( 1 + (sin(lat0) * sin(lat512(i) * 3.14159265359 / 180)) + (cos(lat0) * cos(lat512(i)* 3.14159265359 / 180) * cos((lon512(j)-180)* 3.14159265359 / 180)) )
            enddo
         enddo
         !call write_txt(kFactor, 'kFactor_on_GLgrid', outputfolder)
         call write_txt(reshape(kFactor_test,[nglv,2*nglv]), 'kFactor_on_GLgrid_reshape', outputfolder)
         !kFactor2Inv_to_mali = 1 / kFactor**2
         kFactor2Inv_to_mali = reshape(1 / kFactor_test**2,[nglv,2*nglv])
         call write_txt(kFactor2Inv_to_mali, 'kFactor2Inv_on_GLgrid', outputfolder)

         if (iter*dtime .GT. L_sim) then
            write(unit_num,*) 'ERROR: The current SLM simulation time exceeds the prescribed simulation length'
            write(unit_num,*) 'Current simulation time (in yr) is: ', iter*dtime
            write(unit_num,*) 'Prescribed simulation length (in yr) is: ', L_sim
            write(unit_num,*) 'Check the variable L_sim in namelist.sealevel file.'
            write(unit_num,*) 'Terminating the program.'
            stop
         else
            ! Ice files
            do n = 1, nfiles-1
               j = TIMEWINDOW(n) ! icefile numbers to read in from the TW array
               !write(unit_num,'(A,I6)') 'ice file read in from the SLM output folder, file number:', j
               write(numstr,'(I6)') j

               k = -1 * starttime - j * dt1
               !write(unit_num,'(A,I7)') 'ice load, year (ago):',k ! year 'relative to present'
               numstr = trim(adjustl(numstr))

              ! read in ice files (upto the previous time step) from the sea-level model folder
              call read_sl(icexy(:,:,n), icemodel_out, outputfolder_ice, suffix=numstr)

            enddo

            j = TIMEWINDOW(nfiles) ! icefile number to read in from the TW array
            !write(unit_num,'(A,I6)') 'ice file read in from the coupled input ice folder,  file number,:', j
            write(numstr,'(I6)') j
            numstr = trim(adjustl(numstr))

            ! read in ice thickness at the current time step outside the ISM domain from 'inputfolder_ice'
            call read_sl(icexy(:,:,nfiles), icemodel, inputfolder_ice, suffix=numstr)

            ! ice thickness at the current time step inside the ISM domain is provided by the ISM
            ! merge iceload configuration
            if (patch_ice) then
               do j = 1,2*nglv
                  do i = 1,nglv
                     icexy(i,j,nfiles) = 0.0
                  enddo
               enddo
            endif

            do j = 1,2*nglv
               do i = 1,nglv
                  icexy(i,j,nfiles) = mali_iceload(i,j) + icexy(i,j,nfiles)*(1 - mali_mask(i,j))
               enddo
            enddo
         endif

      else ! if not coupling
         ! if sea level model is not coupled to an ice sheet model, read in iceloads from the folder inputfolder_ice
         do n = 1, nfiles
            j = TIMEWINDOW(n) ! icefile numbers to read in from the TW array
            write(numstr,'(I6)') j
            numstr = trim(adjustl(numstr))
            ! read in ice files (upto the previous time step) from the sea-level model folder
            call read_sl(icexy(:,:,n), icemodel, inputfolder_ice, suffix=numstr)
         enddo

      endif !endif coupling

      !Time array
      if (input_times) then ! time array is inputted from an existing text file, read in and write out
         open(unit = 201, file = trim(inputfolder)//trim(timearray), form = 'formatted', access = 'sequential', status = 'old')
         open(unit = 202, file = trim(outputfolder)//trim(timearray), form = 'formatted', access = 'sequential', &
         & status = 'replace')
         read(201,*) times
         write(202,'(ES14.4E2)') times
         close(201)
         close(202)
      else ! if time array is not read in from a text file, calculate the time array
         ! calculate times that corresponds to ice files that are read in
         do i = 1, nfiles
            times(i) = starttime + TIMEWINDOW(i)*dt1
         enddo
          open(unit = 201, file = trim(outputfolder)//trim(timearray), form = 'formatted', access = 'sequential', &
         & status = 'old', position='append')
         write(201,'(ES14.4E2)') times(nfiles)
         close(201)
      endif

      !Read in the initial topography (topo at the beginning of the full simulation)
      !This is used to output the total sea level change from the beginning of the simulation
      call read_sl(tinit_0, 'tgrid', outputfolder, suffix='0')

      j = TIMEWINDOW(1) ! first element of the time window as the initial file
      write(unit_num,'(A,I4)') 'file number of the first item in the TW:', j
      write(numstr,'(I4)') j
      numstr = trim(adjustl(numstr))

      ! read in initial (first file within the time window) ocean function
      call read_sl(cxy0(:,:), 'ocean', outputfolder, suffix=numstr)

      if (nmelt == 1) then
         cxy(:,:) = cxy0(:,:)
      endif

      ! read in initial (first file within the time window) beta
      call read_sl(beta0(:,:), 'beta', outputfolder, suffix=numstr)

      if (tpw) then
         ! read in variables for the rotation signal
         open(unit = 201, file = trim(outputfolder)//'TPW', form = 'formatted', access = 'sequential', &
         & status = 'old')

         oldlambda(:,:) = (0.0,0.0)
         oldil(:,:) = 0.0
         oldm(:) = 0.0

         do n = 1, nfiles-1
            ! find the number of lines to skip to read in appropriate TPW components
            if (n == 1) then
               j = TIMEWINDOW(n)
            else
               j = TIMEWINDOW(n) - TIMEWINDOW(n-1) - 1
            endif

            !skip lines to read in the rotational components corresponding to timesteps within the TW
            do m = 1, j
               read(201,*) !skip line for il
               read(201,*) !skip reading in mm
               read(201,*) !skip reading in lambda
            enddo

            !read in TPW components - total rotational change from the beginning of simulation
            read(201,'(9ES19.8E2)') ((il(i,j), i = 1,3), j = 1, 3)
            read(201,'(3ES19.8E2)') (mm(i), i = 1, 3)
            read(201,'(18ES19.8E2)') ((lambda(i, j), i = 0, 2), j = 0, 2)

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
         close(201)
      endif !endif (TPW)

      ! find index of the previous timestep
      m = TIMEWINDOW(nfiles-1)
      write(numstr2,'(I4)') m
      numstr2 = trim(adjustl(numstr2))

      ! read in topography of the previous timestep
      call read_sl(topoxy_m1(:,:), 'tgrid', outputfolder, suffix=numstr2)

      ! condition for time window travel
      if (Travel.EQ.0) then
         !if the TW hasnt started travelling, initial topography of the TW is that of the total simulation
         tinit(:,:) = tinit_0(:,:)
      else
         !if The TW is moving, initial topography of the TW is read in
         j = TIMEWINDOW(1)
         write(numstr,'(I4)') j
         numstr = trim(adjustl(numstr))

        call read_sl(tinit(:,:), 'tgrid', outputfolder, suffix=numstr)
      endif
      !endif nmelt>0

      if (nmelt.GT.1) then

         j = TIMEWINDOW(nfiles-1)
         write(unit_num,*) 'Reading in ocean function from previous timestep (file number)', j
         write(numstr,'(I4)') j
         numstr = trim(adjustl(numstr))

         ! read in converged ocean function from the last timestep
         call read_sl(cxy(:,:), 'ocean', outputfolder, suffix=numstr)

         ! if there has been more than one melting episode
         ! read in the total ocean loading computed from previous timesteps 0 to nmelt minus 1.
         ! read in the of total ocean loading changes that are needed.
         do n=1, nfiles-1

            j = TIMEWINDOW(n)
            ! write(unit_num,*) 'reading in converged ocean loading files'
            ! write(unit_num,*) 'sea files, j', j
            write(numstr,'(I4)') j
            numstr = trim(adjustl(numstr))

            open(unit = 201, file = trim(outputfolder)//'dS_converged'//trim(numstr), form = 'formatted', access = 'sequential', &
            & status = 'old')
            read(201,'(ES16.9E2)') dS_converged(:,:)
            close(201)

            deltaS(:,:,n) = dS_converged(:,:)  !save the converged total-ocean loading into a big array

         enddo
      endif

      !-----------------------------------------------------------
      !  Read in Love numbers (Jerry's output from 'maxwell.f')
      !-----------------------------------------------------------
      open(unit = 201, file = trim(planetfolder)//trim(planetmodel), status = 'old')
      ! Following code borrowed from Jerry
      read(201,*)
      do j = 1,norder
         read(201,*) legord(j), nmodes(j)
         nm = nmodes(j)
         xn = real(legord(j))
         ll = legord(j)
         nmod(ll) = nm
         read(201,*) (s(i,ll), i=1, nm)
         read(201,*) elast(1,ll), elast(2,ll), taurr, taurt, elast(3,ll)
         read(201,*) asymv(1,ll), asymv(2,ll), taurr, taurt, asymv(3,ll)
         read(201,*) (resh(i,ll), i=1,nm)
         read(201,*) (resl(i,ll), i=1,nm)
         read(201,*) (resk(i,ll), i=1,nm)
         if (xn .lt. 1.5) cycle
         read(201,*) telast(1,ll), telast(2,ll), dmx, dmy, telast(3,ll)
         read(201,*) tasymv(1,ll), tasymv(2,ll), dmx, dmy, tasymv(3,ll)
         read(201,*) (tresh(i,ll), i=1, nm)
         read(201,*) (tresl(i,ll), i=1, nm)
         read(201,*) (tresk(i,ll), i=1, nm)
      enddo
      close(201)

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

      write(unit_num,*) ''
      write(unit_num,'(A,I4,A,EN15.4E2,A)') '  Timestep',iter,', from ',times(nfiles-1),' years to '
      write(unit_num,'(A,EN15.4E2,A)') '                     ',times(nfiles),' years'

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
      do j = 1, 2 * nglv
         do i = 1, nglv
             if (icexy(i, j, nfiles) < epsilon(0.0)) then
                beta(i,j) = 1
             else
                beta(i,j) = 0
             endif
         enddo
      enddo

      ! calculate initial (initial within the time window) cstar
      cstar0(:,:)  = cxy0(:,:) * beta0(:,:)
      cstarxy(:,:) = cxy(:,:) * beta(:,:)  ! First guess to the O.F using converged O.F the preivous timestep

      call spat2spec(cstarxy,cstarlm,spheredat) ! Decompose the current cstar

      ! Calculate the first guess to topography correction
      tOxy(:,:) = tinit(:,:) * (cstarxy(:,:) - cstar0(:,:)) ! (eq. 70)
      call spat2spec(tOxy, tOlm, spheredat) ! Decompose the topo correction term


      dS(:,:,1) = deltaS(:,:,1)
      if (nmelt > 1) then
         do n = 2, nfiles - 1
            dS(:,:,n) = deltaS(:,:,n) - deltaS(:,:,n-1) !Calculate the ocean loading changes between every time step
         enddo
      endif

      ! Initial guess (eustatic) to the ocean loading change at current time step
      dS(:,:,nfiles) = (cstarlm(:,:) / cstarlm(0,0)) * ((-rhoi / rhow) * dicestar(0,0,nfiles)) ! Save into a big array
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
            write(unit_num,*)
            write(unit_num,'(A,I5,A)') 'WARNING: The inner loop failed to converge after the limit of ', ninner, ' iterations.'
            write(unit_num,'(A,ES15.3,A)') '         The variable xi finished with a value of ', xi, ', resulting from '
            write(unit_num,'(A,ES15.3)') '         sum(abs(olddSlm)) = ', sum(abs(olddSlm))
            write(unit_num,'(A,ES15.3)') '            sum(abs(dSlm)) = ', sum(abs(dSlm))
            write(unit_num,*)
            write(unit_num,'(A)') '!!!---- Program sl_model will now be terminated. ----!!!'
            call abort ! Terminate program
            !exit ! DEBUG line: Continue inner loop despite non-convergence. Normal operation: Enable the 2 lines above.
         endif
         ninner = ninner + 1

      enddo ! End inner loop
      !-----------------------------------------------------------
      !<<<<<<<<<<<<<<<<<<<< End of inner loop >>>>>>>>>>>>>>>>>>>>
      !-----------------------------------------------------------
      write(unit_num,'(A,I4,A)') '  ', ninner, ' inner-loop iterations'

      !HH: print out the number of iteration it takes for the inner convergence
      open(unit = 201, file = trim(outputfolder)//'numiter', form = 'formatted', access ='sequential', &
      & status = 'old', position='append')
      write(201,'(I5)') ninner
      close(201)

      ! Write out the converged rotation-related quantities
      if (tpw) then
         open(unit = 201, file = trim(outputfolder)//'TPW', form = 'formatted', access = 'sequential', &
         & status = 'old', position='append')
         write(201,'(9ES19.8E2/,3ES19.8E2/,18ES19.8E2)') il(:,:), mm(:), lambda(:,:)
         close(201)
      endif

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
         rr(:,:,nfiles) = rrxy(:,:)
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

      write(unit_num,'(A)') 'Writing output files...'
      j = TIMEWINDOW(nfiles)
      write(unit_num,*) 'FILENUMBER of new outputs : ',j
      write(numstr,'(I4)') j
      numstr = trim(adjustl(numstr))

      ! nmelt
      open(unit = 201, file = trim(outputfolder)//'nmelt', form = 'formatted', access ='sequential', &
      & status = 'old',position='append')
      write(201,'(I4)') nmelt
      close(201)

      ! topography at the current timestep
      call write_sl(topoxy, 'tgrid', outputfolder, suffix=numstr)

      ! converged ocean function at the current timestep
      call write_sl(cxy, 'ocean', outputfolder, suffix=numstr)

      ! converged beta function at the current timestpe
      call write_sl(beta, 'beta', outputfolder, suffix=numstr)

      ! output converged total ocean loading changes
      open(unit = 201, file = trim(outputfolder)//'dS_converged'//trim(numstr), form = 'formatted', access = 'sequential', &
      & status = 'replace')
      write(201,'(ES14.7E2)') deltaS(:,:,nfiles)
      close(201)

      current_time = iter * dt1     !time passed since the start of the simulation
      if (current_time == L_sim) then !if we are at the last time step of simulation
         write(iterstr,'(I2)') itersl
           iterstr = trim(adjustl(iterstr))
         write(unit_num,*) 'Last time step of the simulation! writing out files for next outer-iteration loop'

         ! write out the predicted present day topography into a file so it can be used in the next outer-iteration
         call write_sl(topoxy, 'pred_pres_topo_', outputfolder, suffix=iterstr)

         ! write out the initial topography of the simulation at the currect outer-loop into a file
         call write_sl(tinit_0, 'tgrid0_', outputfolder, suffix=iterstr)
       endif

      if (calcRG) then
         call write_sl(rr(:,:,nfiles), 'R', outputfolder, suffix=numstr)

         ! Compute geoid displacement and write out
         gg(:,:,nfiles) = deltaslxy(:,:)+rr(:,:,nfiles)

         call write_sl(gg(:,:,nfiles), 'G', outputfolder, suffix=numstr)
      endif

      if (iceVolume) then
         ice_volume = icestarlm(0,0)*4*pi*radius**2 !multiply the (0,0) component of ice to the area of a sphere
         open(unit = 201, file = trim(outputfolder)//'ice_volume', form = 'formatted', access = 'sequential', &
         & status = 'old', position = 'append')
         write(201,'(ES14.4E2)') ice_volume
         close(201)

         ! calculate the mean SLC
         mean_slc = deltasllm(0,0)
         open(unit = 201, file = trim(outputfolder)//'gmsle_change', form = 'formatted', access ='sequential', &
         & status = 'old', position = 'append')
         write(201,'(ES14.4E2)') mean_slc
         close(201)
      endif

      ! HH: outputs for debugging
      if (coupling) then
      call read_sl(mask_maliMesh, 'mali_mask_init', outputfolder)
      deltaslxy_outside_maliMesh(:,:) = deltaslxy(:,:) * (1 - mask_maliMesh(:,:))
      deltaslxy_inside_maliMesh(:,:) = deltaslxy(:,:) * mask_maliMesh(:,:)
      call write_sl(deltaslxy_outside_maliMesh(:,:), 'deltaslxy_outside_maliMesh', outputfolder, suffix=numstr)
      call write_sl(deltaslxy_inside_maliMesh(:,:), 'deltaslxy_inside_maliMesh', outputfolder, suffix=numstr)
      call spat2spec(deltaslxy_outside_maliMesh(:,:), deltasllm_outside_maliMesh(:,:),spheredat)
      call spat2spec(deltaslxy_inside_maliMesh(:,:), deltasllm_inside_maliMesh(:,:),spheredat)
      mean_slc_outside_maliMesh = deltasllm_outside_maliMesh(0,0)
      mean_slc_inside_maliMesh = deltasllm_inside_maliMesh(0,0)
      mean_slc_sum = mean_slc_outside_maliMesh + mean_slc_inside_maliMesh
      endif

      open(unit = 201, file = trim(outputfolder)//'gmsle_change_outside_maliMesh', form = 'formatted', access ='sequential', &
      & status = 'old', position = 'append')
      write(201,'(ES14.4E2)') mean_slc_outside_maliMesh
      close(201)

      open(unit = 201, file = trim(outputfolder)//'gmsle_change_inside_maliMesh', form = 'formatted', access ='sequential', &
      & status = 'old', position = 'append')
      write(201,'(ES14.4E2)') mean_slc_inside_maliMesh
      close(201)

      open(unit = 201, file = trim(outputfolder)//'gmsle_change_sumMesh', form = 'formatted', access ='sequential', &
      & status = 'old', position = 'append')
      write(201,'(ES14.4E2)') mean_slc_sum
      close(201)

      if (coupling) then
         ! topography change between the previous and the current timestep
         ! this is the information passed to the ice sheet model
         !call write_sl(topoxy_m1(:,:)-topoxy(:,:), 'bedrock', folder_coupled)
          slchange = topoxy_m1(:,:)-topoxy(:,:)
         !write out the current ice load as a new file
         call write_sl(icexy(:,:,nfiles), icemodel_out, outputfolder_ice, suffix=numstr)
      endif !endif coupling

      call system_clock(countf) ! Total time
      call cpu_time(countf_cpu)

      ! Write out total compuatation time of sea level change over current timestep
      open(unit = 201, file = trim(outputfolder)//'elapsed_wall_time', form = 'formatted', access = 'sequential', &
      & status = 'old', position='append')
      write(201,'(ES14.4E2)') float(countf-counti)/float(countrate)
      close(201)

      open(unit = 201, file = trim(outputfolder)//'elapsed_cpu_time', form = 'formatted', access = 'sequential', &
      & status = 'old', position='append')
      write(201,'(ES14.4E2)') countf_cpu-counti_cpu
      close(201)


      write(unit_num,'(A,F7.2,A)') 'Done! Total time ', real(countf - counti) / real(countrate), ' seconds'
      write(unit_num,*) ''
      write(unit_num,*) ''

      if (Travel_total > 0 .and. Travel == Travel_total) then
         write(unit_num,*) ' THE TIMEWINDOW HAS REACHED THE END OF THE SIMULATION!'
         write(unit_num,*) ' GREAT JOB TW!'
      endif
      write(unit_num,*) ''
   end subroutine sl_solver


!----------------------------------------------------------------------------------------------------------------------!



end module sl_model_mod
