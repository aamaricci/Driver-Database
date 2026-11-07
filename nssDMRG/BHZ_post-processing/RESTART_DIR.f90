!This is just of a convenience for the post-processing driver.
!Extract and remove the `restart` directory automatically from the
!post-processing drivers.
!It is essentially a wrapper around to bash commands.
MODULE RESTART_DIR
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  logical :: temporary_restart=.false.
  integer :: restart_status,restart_ierr


contains


  ! Only the master mutates files. Existing restart directories are user-owned.
  subroutine prepare_restart()
    logical :: master=.true.
    logical :: exists
    !
    restart_status=0
    !
#ifdef _MPI
    if(check_MPI())then
      master=get_Master_MPI(MPI_COMM_WORLD)
    endif
#endif
    if(master)then
      inquire(file='restart',exist=exists)
      if(.not.exists)then
        ! finalize_dmrg archives a single top-level restart_<timestamp>/ directory.
        ! Strip that component on extraction, leaving the original archive intact.
        call execute_command_line( &
          'set -- restart_*.tgz; ' // &
          'if [ -f "$1" ]; then ' // &
          '  [ "$#" -eq 1 ] || { echo "Multiple restart archives; prepare restart/ explicitly" >&2; exit 1; }; ' // &
          '  mkdir restart || exit 1; ' // &
          '  tar -xzf "$1" -C restart --strip-components=1 || { rm -rf -- restart; exit 1; }; ' // &
          'else ' // &
          '  set -- restart_*/; ' // &
          '  [ "$#" -eq 1 ] && [ -d "$1" ] || ' // &
          '    { echo "Expected one uncompressed restart_<timestamp>/ directory" >&2; exit 1; }; ' // &
          '  mkdir restart || exit 1; ' // &
          '  cp -R "$1". restart/ || { rm -rf -- restart; exit 1; }; ' // &
          'fi',wait=.true.,exitstat=restart_status,cmdstat=restart_ierr)
        if(restart_ierr/=0)restart_status=1
        temporary_restart=restart_status==0
      endif
    endif
#ifdef _MPI
    call Bcast_MPI(MPI_COMM_WORLD,restart_status)
    call Bcast_MPI(MPI_COMM_WORLD,temporary_restart)
#endif
    if(restart_status/=0)error stop 'Cannot prepare restart; original checkpoint retained'
  end subroutine prepare_restart







  subroutine cleanup_restart()
    logical :: master=.true.
    !
#ifdef _MPI
    if(check_MPI())then
      master=get_Master_MPI(MPI_COMM_WORLD)
    endif
    call Barrier_MPI(MPI_COMM_WORLD)
#endif
    !
    restart_status=0
    if(master.and.temporary_restart)then
      call execute_command_line('rm -rf -- restart',wait=.true., &
        exitstat=restart_status,cmdstat=restart_ierr)
      if(restart_ierr/=0)restart_status=1
    endif
#ifdef _MPI
    call Bcast_MPI(MPI_COMM_WORLD,restart_status)
#endif
    if(restart_status/=0)error stop 'Cannot remove temporary restart directory'
  end subroutine cleanup_restart





end module RESTART_DIR
