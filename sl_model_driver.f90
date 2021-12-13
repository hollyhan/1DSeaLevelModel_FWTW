program sl_model_driver
	use sl_model_mod 
	implicit none


integer :: i, itersl_sh, iter_sh, dtime_sh
real :: starttime_sh                                  ! start time of the simulation 
integer :: iargc, nargs                                 ! Arguments read in from a bash script                              
character(16) :: carg(20)                               ! Arguments from a bash script

! Reading in arguments from a bash script
nargs = iargc()
do i=1,nargs
   call getarg(i, carg(i))
enddo
read (carg(1),*) itersl_sh
read (carg(2),*) iter_sh      ! the coupling time step we are on (in years)
read (carg(3),*) dtime_sh     ! coupling time (in years)
read (carg(4),*) starttime_sh ! start time of the simulation (in years)

! set up output unit number
call sl_set_unit_num(6)

! check point for time array and coupling
call sl_solver_checkpoint(itersl_sh, dtime_sh)

! set up the temporal resolution 
call sl_timewindow(iter_sh)

! intialize and execute the sea-level solver
if (iter_sh .eq. 0) then 
	call sl_solver_init(itersl_sh, starttime_sh)
elseif (iter_sh .gt. 0) then 
	call sl_solver(itersl_sh, iter_sh, dtime_sh, starttime_sh)
endif

end program sl_model_driver