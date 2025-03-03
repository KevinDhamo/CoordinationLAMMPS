module benchmark_mod
    use types_mod
    use config_mod
    use cell_list_mod
    use coordination_mod
    use io_mod
    use error_mod
    use omp_lib
    implicit none
    private

    ! Public procedures
    public :: run_benchmark
    public :: initialize_benchmark
    public :: cleanup_benchmark

    ! Benchmark parameters
    integer, parameter :: BENCHMARK_TIME_LIMIT = 30  ! seconds
    character(len=*), parameter :: BENCHMARK_OUTPUT = 'benchmark_results.txt'
    
    ! Module variables
    integer :: benchmark_unit = -1
    logical :: benchmark_file_open = .false.
    real(dp) :: single_thread_fps = 0.0_dp
    
    ! Type for storing benchmark results
    type :: benchmark_result_type
        integer :: n_threads
        real(dp) :: total_time
        real(dp) :: frames_per_second
        integer :: frames_processed
        real(dp) :: speedup
        real(dp) :: efficiency
    end type benchmark_result_type
    
contains
    subroutine initialize_benchmark()
        integer :: io_stat
        
        ! Initialize module variables
        single_thread_fps = 0.0_dp
        
        ! Open benchmark output file
        open(newunit=benchmark_unit, file=BENCHMARK_OUTPUT, status='replace', &
             action='write', iostat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to open benchmark file", ERR_FILE_IO)
        end if
        benchmark_file_open = .true.
        
        ! Write header
        write(benchmark_unit,'(A)') "# Coordination Analysis Benchmark Results"
        write(benchmark_unit,'(A)') "# Threads, Time(s), Frames/s, Frames, Speedup, Efficiency"
        
        ! Disable output file writing during benchmarks
        disable_output = .true.
    end subroutine initialize_benchmark

    subroutine init_coordination_defaults(n_atoms, n_types, atom_info)
        integer, intent(in) :: n_atoms, n_types
        type(atom_type_info), intent(in) :: atom_info(:)
        integer :: i, j, pair_idx
        
        ! Calculate number of unique pairs
        n_pairs = 0
        do i = 1, n_types
            if (.not. atom_info(i)%include) cycle
            do j = i, n_types
                if (.not. atom_info(j)%include) cycle
                n_pairs = n_pairs + 1
            end do
        end do
        
        ! Allocate arrays
        if (allocated(pairs)) deallocate(pairs)
        if (allocated(coord_numbers)) deallocate(coord_numbers)
        allocate(pairs(n_pairs), coord_numbers(n_atoms, n_pairs))
        
        ! Initialize pair information with default cutoffs
        pair_idx = 1
        do i = 1, n_types
            if (.not. atom_info(i)%include) cycle
            do j = i, n_types
                if (.not. atom_info(j)%include) cycle
                
                pairs(pair_idx)%type1_id = i
                pairs(pair_idx)%type2_id = j
                pairs(pair_idx)%cutoff = DEFAULT_CUTOFF
                pairs(pair_idx)%cutoff_sq = DEFAULT_CUTOFF**2
                pair_idx = pair_idx + 1
            end do
        end do
        
        ! Initialize coordination numbers
        coord_numbers = 0
    end subroutine init_coordination_defaults

    subroutine run_benchmark(coords, atom_types, elements, box_length, n_atoms, &
                            n_types, atom_info, total_frames)
        real(dp), intent(inout) :: coords(:,:)
        integer, intent(inout) :: atom_types(:)
        character(len=*), intent(inout) :: elements(:)
        real(dp), intent(inout) :: box_length(3)
        integer, intent(in) :: n_atoms, n_types, total_frames
        type(atom_type_info), intent(in) :: atom_info(:)
        
        integer :: max_threads, n_threads, frame
        logical :: eof
        real(dp) :: start_time, current_time
        type(benchmark_result_type) :: result
        logical, allocatable :: include_mask(:)
        
        ! Get maximum available threads
        max_threads = omp_get_max_threads()
        
        ! Allocate mask array
        allocate(include_mask(n_types))
        include_mask = atom_info%include
        
        write(*,'(A)') "=== Starting Benchmark ==="
        write(*,'(A,I0,A)') "Maximum available threads: ", max_threads
        
        ! Test different thread counts (powers of 2)
        n_threads = 1
        do while (n_threads <= max_threads)
            write(*,'(A,I0,A)') "Testing with ", n_threads, " threads..."
            
            ! Set OpenMP parameters
            call omp_set_num_threads(n_threads)
            call omp_set_schedule(omp_sched_guided, 0)  ! Use guided scheduling
            
            ! Set OpenMP environment variables for better performance
            if (n_threads > 1) then
                call execute_command_line("export OMP_PROC_BIND=close", wait=.false.)
                call execute_command_line("export OMP_PLACES=cores", wait=.false.)
            end if
            
            ! Initialize for this run
            frame = 0
            eof = .false.
            call init_coordination_defaults(n_atoms, n_types, atom_info)
            call initialize_cell_list(box_length, get_cutoffs(), n_atoms)
            
            ! Initialize IO for this run
            call cleanup_io()
            call initialize_io(INPUT_FILE, OUTPUT_FILE, total_frames)
            
            ! Start timing
            start_time = omp_get_wtime()
            
            ! Main processing loop
            do while (.not. eof)
                call read_trajectory_frame(coords, atom_types, elements, &
                                        box_length, frame, eof, atom_info)
                if (eof) exit
                
                call calculate_coordination(coords, atom_types, n_atoms, &
                                         box_length, frame, include_mask)
                
                ! Update progress (with reduced frequency for benchmark)
                if (mod(frame, 10) == 0) then
                    write(*,'(a1,1x,A,I0,A,I0,A)',advance='no') achar(13), &
                        'Processing frame ', frame, ' of ', total_frames, '    '
                end if
                
                ! Check if we've reached time limit and processed enough frames
                current_time = omp_get_wtime()
                if ((current_time - start_time) >= BENCHMARK_TIME_LIMIT .and. frame >= 10) exit
            end do
            
            ! Clear progress line
            write(*,'(a1,A)',advance='yes') achar(13), repeat(' ', 80)
            
            ! Record results
            result%n_threads = n_threads
            result%total_time = omp_get_wtime() - start_time
            result%frames_processed = frame
            
            ! Calculate performance metrics
            if (result%total_time > 0.0_dp) then
                result%frames_per_second = real(result%frames_processed) / result%total_time
            else
                result%frames_per_second = 0.0_dp
            end if
            
            ! Calculate speedup and efficiency
            if (n_threads == 1) then
                single_thread_fps = result%frames_per_second
                result%speedup = 1.0_dp
                result%efficiency = 1.0_dp
            else if (single_thread_fps > 0.0_dp) then
                result%speedup = result%frames_per_second / single_thread_fps
                result%efficiency = result%speedup / real(n_threads, dp)
            else
                result%speedup = 0.0_dp
                result%efficiency = 0.0_dp
            end if
            
            ! Write results
            call write_benchmark_result(result)
            
            ! Cleanup for this run
            call cleanup_coordination()
            call cleanup_cell_list()
            
            ! Double the thread count for next iteration
            if (n_threads == max_threads) exit
            n_threads = min(n_threads * 2, max_threads)
        end do
        
        ! Final cleanup
        call cleanup_io()
        deallocate(include_mask)
        write(*,'(A)') "Benchmark complete. Results written to: " // BENCHMARK_OUTPUT
        write(*,'(A)') "Note: Performance difference between benchmark and normal run is expected"
        write(*,'(A)') "      due to output file writing being disabled during benchmark."
    end subroutine run_benchmark

    subroutine write_benchmark_result(result)
        type(benchmark_result_type), intent(in) :: result
        
        ! Write to screen
        write(*,'(A,I0,A)') "Results for ", result%n_threads, " threads:"
        write(*,'(A,F10.3,A)') "  Time: ", result%total_time, " seconds"
        write(*,'(A,I0)') "  Frames processed: ", result%frames_processed
        write(*,'(A,F10.3,A)') "  Performance: ", result%frames_per_second, " frames/second"
        write(*,'(A,F10.3)') "  Speedup: ", result%speedup
        write(*,'(A,F10.3)') "  Efficiency: ", result%efficiency
        write(*,*)
        
        ! Write to file
        write(benchmark_unit,'(I6,5(F12.3))') result%n_threads, result%total_time, &
            result%frames_per_second, real(result%frames_processed), &
            result%speedup, result%efficiency
        call flush(benchmark_unit)
    end subroutine write_benchmark_result

    subroutine cleanup_benchmark()
        if (benchmark_file_open) then
            close(benchmark_unit)
            benchmark_file_open = .false.
        end if
        disable_output = .false.  ! Re-enable output for normal operation
    end subroutine cleanup_benchmark

    function get_cutoffs() result(cutoffs)
        real(dp), allocatable :: cutoffs(:)
        integer :: i
        
        allocate(cutoffs(n_pairs))
        do i = 1, n_pairs
            cutoffs(i) = pairs(i)%cutoff
        end do
    end function get_cutoffs

end module benchmark_mod
