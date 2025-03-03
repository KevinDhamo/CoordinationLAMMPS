program coordination_analysis
    use types_mod
    use config_mod
    use cell_list_mod
    use coordination_mod
    use io_mod
    use error_mod
    use benchmark_mod
    use omp_lib
    implicit none

    ! System variables
    integer :: n_atoms = 0
    integer :: n_types = 0
    integer :: n_types_actual = 0  ! For verification
    type(atom_type_info), allocatable :: atom_info(:)
    real(dp) :: box_length(3)
    
    ! Atomic data
    real(dp), allocatable :: coords(:,:)
    integer, allocatable :: atom_types(:)
    character(len=8), allocatable :: elements(:)
    
    ! Analysis variables
    integer :: frame = 0
    integer :: total_frames = 0
    logical :: eof = .false.
    real(dp) :: start_time, end_time
    logical, allocatable :: include_mask(:)
    
    ! Local variables
    integer :: io_stat, num_threads
    logical :: benchmark_mode = .false.
    character(len=32) :: arg
    logical :: setup_file_read = .false.
    
    ! Check command line arguments for benchmark mode
    if (command_argument_count() > 0) then
        call get_command_argument(1, value=arg)
        benchmark_mode = (trim(arg) == '--benchmark')
    end if
    
    ! Initialize configuration
    call initialize_config()
    call read_config()
    
    ! First, read DATA file to get atom count regardless of setup file 
    ! This ensures we have the correct atom count to allocate arrays 
    call read_data_file_dims_only(DATA_FILE, n_atoms, n_types_actual, box_length)
    
    if (VERBOSE) then
        write(*,'(A)') " Data file dimensions loaded:"
        write(*,'(A,I0)') "   Actual atoms: ", n_atoms 
        write(*,'(A,I0)') "   Actual atom types: ", n_types_actual
    end if
    
    ! Try to read setup file for additional configuration
    if (.not. benchmark_mode) then
        setup_file_read = read_setup_file(atom_info, n_types, pairs, n_pairs)
        
        if (setup_file_read) then
            ! Check if number of types matches
            if (n_types /= n_types_actual) then
                write(*,'(A,I0,A,I0)') " WARNING: Number of atom types in setup file (", &
                                       n_types, ") doesn't match data file (", n_types_actual, ")"
            end if
            
            if (requested_cores > 0) then
                ! Use number of cores from setup file
                num_threads = min(requested_cores, omp_get_max_threads())
            else
                ! Use all available cores
                num_threads = omp_get_max_threads()
            end if
            
            ! Set OpenMP parameters
            call omp_set_num_threads(num_threads)
            call omp_set_schedule(omp_sched_guided, 0)  ! Use guided scheduling
            if (VERBOSE) write(*,'(A,I4,A)') ' Using', num_threads, ' OpenMP threads'
            
            ! Set environment variables for better thread affinity
            if (num_threads > 1) then
                call execute_command_line("export OMP_PROC_BIND=close", wait=.false.)
                call execute_command_line("export OMP_PLACES=cores", wait=.false.)
            end if
        end if
    end if
    
    ! If setup file was not read, or in benchmark mode, use traditional approach for atom types
    if (.not. setup_file_read) then
        ! Read atom types from data file
        call read_data_file(DATA_FILE, atom_info, n_atoms, n_types, box_length)
        n_types = n_types_actual  ! Ensure consistency
        
        ! Set OpenMP threads (not in benchmark mode)
        if (.not. benchmark_mode) then
            num_threads = omp_get_max_threads()
            call omp_set_num_threads(num_threads)
            call omp_set_schedule(omp_sched_guided, 0)  ! Use guided scheduling
            if (VERBOSE) write(*,'(A,I4,A)') ' Using', num_threads, ' OpenMP threads'
        end if
    end if
    
    ! Allocate arrays with error checking - use n_atoms from data file
    allocate(coords(n_atoms,3), atom_types(n_atoms), elements(n_atoms), &
             include_mask(n_types), stat=io_stat)
    if (io_stat /= 0) call handle_error("Failed to allocate arrays", ERR_ALLOCATION)
    
    ! Initialize arrays
    coords = 0.0_dp
    atom_types = 0
    elements = ''
    include_mask = atom_info%include
    
    ! Count frames in trajectory file
    if (VERBOSE) write(*,'(A)') ' Counting frames in trajectory file...'
    
    ! Use file from setup if available
    if (setup_file_read) then
        total_frames = count_trajectory_frames(input_trajectory)
    else 
        total_frames = count_trajectory_frames(INPUT_FILE)
    end if
    
    ! Apply frame selection based on setup.txt
    if (setup_file_read) then
        ! Apply start_frame and end_frame parameters from setup file
        if (start_frame > 0) then
            if (VERBOSE) write(*,'(A,I0)') ' Starting at frame: ', start_frame
        end if
        
        if (end_frame > 0) then
            total_frames = min(total_frames, end_frame)
            if (VERBOSE) write(*,'(A,I0)') ' Will process up to frame: ', end_frame
        end if
    end if
    
    if (VERBOSE) write(*,'(A,I0,A)') ' Will process ', total_frames, ' frames'
    
    if (benchmark_mode) then
        ! Run benchmark mode
        call initialize_benchmark()
        call run_benchmark(coords, atom_types, elements, box_length, &
                         n_atoms, n_types, atom_info, total_frames)
        call cleanup_benchmark()
    else
        ! Normal analysis mode
        ! Set up analysis
        call initialize_coordination(n_atoms, n_types, atom_info)
        call initialize_cell_list(box_length, get_cutoffs(), n_atoms)
        
        ! Initialize I/O
        if (setup_file_read) then
            call initialize_io(input_trajectory, output_data_file, total_frames, start_frame)
            frame = start_frame  ! Set starting frame if specified
        else
            call initialize_io(INPUT_FILE, OUTPUT_FILE, total_frames)
        end if
        
        ! Process trajectory
        start_time = omp_get_wtime()
        
        do while (.not. eof)
            call read_trajectory_frame(coords, atom_types, elements, box_length, &
                                    frame, eof, atom_info)
            if (eof) exit
            
            call calculate_coordination(coords, atom_types, n_atoms, box_length, &
                                     frame, include_mask)
            
            call write_output_frame(frame, elements, atom_types, coord_numbers, n_atoms, atom_info)
            call update_progress(frame, total_frames)
            
            ! Check if we've reached the end_frame limit
            if (setup_file_read .and. end_frame > 0 .and. frame >= end_frame) exit
        end do
        
        ! Finalize
        end_time = omp_get_wtime()
        
        ! Write final report before cleanup
        call write_final_report(frame, n_atoms, end_time - start_time, atom_info, &
                              n_types, coord_numbers, atom_types)
        
        ! Cleanup analysis
        call cleanup_io()
        call cleanup_coordination()
        call cleanup_cell_list()
    end if
    
    ! Final cleanup
    call cleanup_config()
    
    if (allocated(coords)) deallocate(coords)
    if (allocated(atom_types)) deallocate(atom_types)
    if (allocated(elements)) deallocate(elements)
    if (allocated(include_mask)) deallocate(include_mask)
    if (allocated(atom_info)) deallocate(atom_info)

contains
    function get_cutoffs() result(cutoffs)
        real(dp), allocatable :: cutoffs(:)
        integer :: i
        
        allocate(cutoffs(n_pairs))
        do i = 1, n_pairs
            cutoffs(i) = pairs(i)%cutoff
        end do
    end function get_cutoffs

end program coordination_analysis
