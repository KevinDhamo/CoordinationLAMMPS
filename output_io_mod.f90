module output_io_mod
    use types_mod
    use config_mod
    use error_mod
    use coordination_mod, only: pairs, n_pairs
    use cell_list_mod, only: num_cell_resets, num_cell_updates, num_selective_updates, &
                            get_cell_statistics, total_empty_cells, &
                            max_atoms_in_cell, min_atoms_in_cell, &
                            avg_atoms_per_cell
    use omp_lib
    implicit none
    private

    ! Public procedures
    public :: init_output
    public :: write_output_frame
    public :: write_final_report
    public :: cleanup_output
    
    ! Output control parameters
    logical, public :: disable_output = .false.  ! New: Option to disable output file writing

    ! Module variables
    integer :: output_unit = -1
    logical :: output_file_open = .false.

contains
    subroutine init_output(filename)
        character(len=*), intent(in) :: filename
        integer :: io_stat, i
        
        ! If output is disabled, do nothing
        if (disable_output) return
        
        ! Open output file
        open(newunit=output_unit, file=filename, status='replace', &
             action='write', iostat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to open output file: " // trim(filename), ERR_FILE_IO)
        end if
        output_file_open = .true.
        
        ! Write header
        write(output_unit,'(A)') "# Coordination Number Analysis"
        write(output_unit,'(A)') "# Column format:"
        write(output_unit,'(A)') "#   1. Atom ID"
        write(output_unit,'(A)') "#   2. Atom Type"
        do i = 1, n_pairs
            write(output_unit,'(A,I0,A,I0,A,I0)') "#   ", i+2, &
                ". CN(Type ", pairs(i)%type1_id, "-", pairs(i)%type2_id, ")"
        end do
        write(output_unit,'(A)') "# Frame data follows"
    end subroutine init_output

    subroutine write_output_frame(frame_number, elements, atom_types, coord_numbers_local, n_atoms, atom_info)
        integer, intent(in) :: frame_number
        character(len=*), intent(in) :: elements(:)
        integer, intent(in) :: atom_types(:)
        integer, intent(in) :: coord_numbers_local(:,:)
        integer, intent(in) :: n_atoms
        type(atom_type_info), intent(in) :: atom_info(:)
        
        integer :: i, j
        character(len=32) :: fmt_str
        
        ! Skip if output is disabled
        if (disable_output .or. .not. output_file_open) return
        
        ! Write frame header
        write(output_unit,'(A,I0)') "# Frame: ", frame_number
        
        ! Construct format string based on number of pairs
        write(fmt_str,'(A,I0,A)') '(I6,A15,', n_pairs, 'I10)'
        
        ! Write atom data with type names
        do i = 1, n_atoms
            if (atom_types(i) >= 1 .and. atom_types(i) <= size(atom_info)) then
                write(output_unit, fmt_str) i, trim(atom_info(atom_types(i))%name), &
                                            (coord_numbers_local(i,j), j=1,n_pairs)
            else
                write(output_unit, fmt_str) i, "Unknown", (coord_numbers_local(i,j), j=1,n_pairs)
            end if
        end do
    end subroutine write_output_frame

    subroutine write_final_report(n_frames, n_atoms, total_time, atom_info, &
                                 n_types, coord_numbers_local, atom_types_local)
        integer, intent(in) :: n_frames, n_atoms, n_types
        real(dp), intent(in) :: total_time
        type(atom_type_info), intent(in) :: atom_info(:)
        integer, intent(in) :: coord_numbers_local(:,:)
        integer, intent(in) :: atom_types_local(:)
        
        real(dp), allocatable :: type_avg_coord(:,:)
        integer, allocatable :: type_counts(:)
        real(dp) :: frames_per_second
        integer :: i, j, alloc_stat
        real(dp) :: cell_stats(4)
        
        ! Allocate arrays for statistics with error checking
        allocate(type_avg_coord(n_types, n_pairs), type_counts(n_types), stat=alloc_stat)
        if (alloc_stat /= 0) then
            call handle_error("Failed to allocate statistics arrays", ERR_ALLOCATION)
            return
        end if
        
        ! Initialize arrays
        type_avg_coord = 0.0_dp
        type_counts = 0
        
        ! Calculate statistics safely
        !$OMP PARALLEL DO REDUCTION(+:type_counts,type_avg_coord) PRIVATE(i,j) SCHEDULE(guided)
        do i = 1, n_atoms
            if (atom_types_local(i) > 0 .and. atom_types_local(i) <= n_types) then
                if (atom_info(atom_types_local(i))%include) then
                    type_counts(atom_types_local(i)) = type_counts(atom_types_local(i)) + 1
                    do j = 1, n_pairs
                        type_avg_coord(atom_types_local(i),j) = &
                            type_avg_coord(atom_types_local(i),j) + coord_numbers_local(i,j)
                    end do
                end if
            end if
        end do
        !$OMP END PARALLEL DO
        
        ! Calculate averages
        do i = 1, n_types
            if (type_counts(i) > 0) then
                type_avg_coord(i,:) = type_avg_coord(i,:) / type_counts(i)
            end if
        end do
        
        ! Get cell list statistics safely
        cell_stats = get_cell_statistics()
        
        ! Write performance report
        write(*,*)
        write(*,'(A)') '=== Analysis Complete ==='
        write(*,'(A,I0,A,I0,A)') ' Processed ', n_frames, ' frames with ', &
                                 n_atoms, ' atoms each'
        write(*,'(A,F10.3,A)') ' Total time: ', total_time, ' seconds'
        frames_per_second = real(n_frames)/max(total_time, 1.0e-6_dp)
        write(*,'(A,F10.3,A)') ' Performance: ', frames_per_second, ' frames/second'
        
        ! Write cell list statistics
        write(*,*)
        write(*,'(A)') 'Cell List Performance:'
        write(*,'(A,I0)') '  Total cell resets: ', num_cell_resets
        write(*,'(A,I0)') '  Total cell updates: ', num_cell_updates
        if (selective_rebuild) then
            write(*,'(A,I0)') '  Selective cell updates: ', num_selective_updates
        end if
        write(*,'(A,I0)') '  Average empty cells: ', total_empty_cells
        write(*,'(A,F8.2)') '  Average atoms per non-empty cell: ', avg_atoms_per_cell
        write(*,'(A,I0,A,I0)') '  Min/Max atoms per cell: ', &
                               min_atoms_in_cell, '/', max_atoms_in_cell
        
        ! Write coordination statistics
        write(*,*)
        write(*,'(A)') 'Coordination Statistics by Type:'
        do i = 1, n_types
            if (atom_info(i)%include) then
                write(*,'(A,A,A,I0,A,I0,A)') ' Type ', trim(atom_info(i)%name), &
                    ' (', i, ', ', type_counts(i), ' atoms):'
                do j = 1, n_pairs
                    if (pairs(j)%type1_id == i .or. pairs(j)%type2_id == i) then
                        write(*,'(A,A,A,I0,A,A,A,I0,A,F10.3)') &
                            '   With ', trim(atom_info(pairs(j)%type1_id)%name), &
                            '-', pairs(j)%type1_id, ' to ', &
                            trim(atom_info(pairs(j)%type2_id)%name), '-', &
                            pairs(j)%type2_id, ': ', type_avg_coord(i,j)
                    end if
                end do
            end if
        end do
        
        ! Write composition summary
        write(*,*)
        write(*,'(A)') 'System Composition:'
        do i = 1, n_types
            if (atom_info(i)%include) then
                write(*,'(A,A,A,I0,A,I0,A,F6.2,A)') ' ', &
                    trim(atom_info(i)%name), ' (Type-', i, '): ', &
                    type_counts(i), ' atoms (', &
                    (real(type_counts(i), dp) / real(n_atoms, dp)) * 100.0_dp, '%)'
            end if
        end do
        
        ! Write output status
        if (disable_output) then
            write(*,*)
            write(*,'(A)') 'Note: Output file writing was disabled (e.g., for benchmarking)'
        end if
        
        ! Clean up arrays
        if (allocated(type_avg_coord)) deallocate(type_avg_coord)
        if (allocated(type_counts)) deallocate(type_counts)
    end subroutine write_final_report

    subroutine cleanup_output()
        if (output_file_open) then
            close(output_unit)
            output_file_open = .false.
        end if
    end subroutine cleanup_output

end module output_io_mod
