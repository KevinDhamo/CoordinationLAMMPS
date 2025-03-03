module io_mod
    use types_mod
    use config_mod, only: VERBOSE, INPUT_FILE, OUTPUT_FILE  ! Add VERBOSE to the import list
    use trajectory_io_mod, only: init_trajectory_reader, read_trajectory_frame, &
                                count_trajectory_frames, cleanup_trajectory_reader, &
                                seek_to_frame
    use data_file_io_mod, only: read_data_file, read_data_file_dims_only
    use output_io_mod, only: init_output, write_output_frame, &
                           write_final_report, cleanup_output, disable_output
    use progress_mod, only: init_progress, update_progress, &
                          finalize_progress
    implicit none
    private

    ! Re-export only necessary procedures
    public :: read_data_file          ! From data_file_io_mod
    public :: read_data_file_dims_only ! From data_file_io_mod
    public :: read_trajectory_frame   ! From trajectory_io_mod
    public :: seek_to_frame           ! From trajectory_io_mod
    public :: write_output_frame      ! From output_io_mod
    public :: write_final_report      ! From output_io_mod
    public :: count_trajectory_frames ! From trajectory_io_mod
    public :: update_progress         ! From progress_mod
    public :: disable_output          ! From output_io_mod
    
    ! Initialize and cleanup procedures
    public :: initialize_io
    public :: cleanup_io

contains
    subroutine initialize_io(trajectory_file, output_file, n_frames, start_at_frame)
        character(len=*), intent(in) :: trajectory_file, output_file
        integer, intent(in) :: n_frames
        integer, intent(in), optional :: start_at_frame
        integer :: frame_to_seek = 0
        
        ! Initialize trajectory reader
        call init_trajectory_reader(trajectory_file)
        
        ! Initialize output file
        call init_output(output_file)
        
        ! Initialize progress bar
        call init_progress(n_frames)
        
        ! Seek to starting frame if specified
        if (present(start_at_frame)) then
            frame_to_seek = start_at_frame
            if (frame_to_seek > 0) then
                call seek_to_frame(frame_to_seek)
                if (VERBOSE) then
                    write(*,'(A,I0)') ' Seeking to frame ', frame_to_seek
                end if
            end if
        end if
    end subroutine initialize_io

    subroutine cleanup_io()
        call cleanup_trajectory_reader()
        call cleanup_output()
        call finalize_progress()
    end subroutine cleanup_io

end module io_mod
