!include 'spharmt.f90' ! Spherical harmonic transform module
!include 'init_slm_mod.f90'
!include 'user_specs_mod.f90'
!include 'sl_model_mod.f90' 

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

call sl_solver(itersl_sh, iter_sh, dtime_sh, starttime_sh)

end program sl_model_driver