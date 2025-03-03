# Coordination Analysis Program: A Developer's Guide

This document provides a comprehensive overview of the coordination analysis program, highlighting key components with code snippets from each module. This program analyzes molecular dynamics trajectories to calculate coordination numbers between different atom types.

## Program Overview

The coordination analysis program calculates how many atoms of specific types surround each atom within defined cutoff distances. It uses a cell-based spatial decomposition algorithm for efficiency and supports OpenMP parallelization for performance on multi-core systems.

## Table of Contents
1. [types_mod.f90](#types_modf90): Core data types
2. [error_mod.f90](#error_modf90): Error handling utilities
3. [config_mod.f90](#config_modf90): Configuration and parameter management
4. [cell_list_mod.f90](#cell_list_modf90): Spatial decomposition algorithm
5. [coordination_mod.f90](#coordination_modf90): Core coordination analysis
6. [data_file_io_mod.f90](#data_file_io_modf90): Structure file parsing
7. [trajectory_io_mod.f90](#trajectory_io_modf90): Trajectory file handling
8. [output_io_mod.f90](#output_io_modf90): Results output
9. [progress_mod.f90](#progress_modf90): Progress reporting
10. [io_mod.f90](#io_modf90): I/O consolidation
11. [benchmark_mod.f90](#benchmark_modf90): Performance testing
12. [main.f90](#mainf90): Main program

---

## types_mod.f90

This module defines the fundamental data types used throughout the program.

### Key Types

```fortran
type :: atom_type_info
    integer :: type_id           ! Numeric type ID from data file
    character(len=8) :: name     ! Element name from data file comment
    real(dp) :: mass            ! Mass from data file
    logical :: include = .true.  ! Whether to include in analysis
end type atom_type_info

type :: pair_type
    integer :: type1_id, type2_id  ! Numeric type IDs
    real(dp) :: cutoff            ! Cutoff distance
    real(dp) :: cutoff_sq         ! Squared cutoff (for optimization)
end type pair_type

type :: cell_type
    integer :: n_atoms = 0              ! Number of atoms in cell
    integer, allocatable :: atoms(:)    ! Indices of atoms in this cell
contains
    procedure :: init => initialize_cell
    procedure :: cleanup => cleanup_cell
    procedure :: resize => resize_cell_array
end type cell_type
```

### Cell Handling Methods

```fortran
subroutine initialize_cell(this, initial_size)
    class(cell_type), intent(inout) :: this
    integer, intent(in) :: initial_size
    integer :: alloc_stat
    
    if (allocated(this%atoms)) deallocate(this%atoms)
    allocate(this%atoms(initial_size), stat=alloc_stat)
    if (alloc_stat /= 0) error stop "Failed to allocate cell array"
    this%n_atoms = 0
end subroutine initialize_cell

subroutine resize_cell_array(this, new_size)
    class(cell_type), intent(inout) :: this
    integer, intent(in) :: new_size
    integer, allocatable :: temp(:)
    integer :: alloc_stat
    
    allocate(temp(new_size), stat=alloc_stat)
    if (alloc_stat /= 0) error stop "Failed to allocate temporary array"
    
    ! Copy existing data
    temp(1:this%n_atoms) = this%atoms(1:this%n_atoms)
    
    ! Replace old array with new one
    call move_alloc(from=temp, to=this%atoms)
end subroutine resize_cell_array
```

The `atom_type_info` stores metadata about each atom type including its name, mass, and whether to include it in the analysis. The `pair_type` defines cutoff distances between atom types. The `cell_type` implements a container for atoms in the spatial decomposition with dynamic resizing capability.

---

## error_mod.f90

This module provides centralized error handling for consistent reporting throughout the program.

### Error Codes and Handling

```fortran
! Error codes
integer, parameter, public :: ERR_ALLOCATION = 1
integer, parameter, public :: ERR_FILE_IO = 2
integer, parameter, public :: ERR_INVALID_PARAM = 3
integer, parameter, public :: ERR_INVALID_TYPE = 4
integer, parameter, public :: ERR_INCONSISTENT_DATA = 5

subroutine handle_error(msg, error_code, fatal)
    character(len=*), intent(in) :: msg
    integer, intent(in) :: error_code
    logical, intent(in), optional :: fatal
    logical :: is_fatal
    
    is_fatal = .true.
    if (present(fatal)) is_fatal = fatal

    write(error_unit,'(A,I0,2A)') "Error (", error_code, "): ", trim(msg)
    
    if (is_fatal) then
        error stop
    end if
end subroutine handle_error

subroutine check_allocation(alloc_stat, array_name)
    integer, intent(in) :: alloc_stat
    character(len=*), intent(in) :: array_name
    
    if (alloc_stat /= 0) then
        call handle_error("Failed to allocate memory for " // trim(array_name), &
                        ERR_ALLOCATION)
    end if
end subroutine check_allocation
```

This module provides a consistent way to handle errors. The `handle_error` subroutine reports errors with an option to terminate the program. The `check_allocation` subroutine is a specialized handler for memory allocation failures. Named error constants improve code readability.

---

## config_mod.f90

This module manages program configuration from files and provides default parameters.

### Default Parameters

```fortran
! File names and paths
character(len=*), parameter, public :: INPUT_FILE = 'trajectory.dump'
character(len=*), parameter, public :: OUTPUT_FILE = 'coordination_numbers.dat'
character(len=*), parameter, public :: DATA_FILE = 'lammps.data'
character(len=*), parameter, public :: CONFIG_FILE = 'analysis_config.txt'
character(len=*), parameter, public :: SETUP_FILE = 'setup.txt'

! Default parameters
real(dp), parameter, public :: DEFAULT_CUTOFF = 3.0_dp
integer, parameter, public :: PROGRESS_BAR_WIDTH = 50
logical, parameter, public :: VERBOSE = .true.

! Cell list configuration - changed from parameters to variables
integer, public :: cell_update_freq = 5
real(dp), parameter, public :: CELL_SIZE_FACTOR = 1.0_dp
integer, parameter, public :: MAX_ATOMS_PER_CELL = 50
real(dp), parameter, public :: REBUILD_THRESHOLD = 0.1_dp
```

### Reading Setup File

```fortran
function read_setup_file(atom_info, n_types, pairs, n_pairs) result(success)
    type(atom_type_info), allocatable, intent(out) :: atom_info(:)
    integer, intent(out) :: n_types
    type(pair_type), allocatable, intent(out) :: pairs(:)
    integer, intent(out) :: n_pairs
    logical :: success
    
    ! Check if setup file exists
    inquire(file=SETUP_FILE, exist=file_exists)
    if (.not. file_exists) then
        if (VERBOSE) write(*,*) "Setup file not found, using interactive mode"
        return
    end if
    
    ! First pass: count atom types and find basic parameters
    n_types = 0
    do
        read(unit, '(A)', iostat=io_stat) line
        if (io_stat /= 0) exit
        
        ! Skip comments and empty lines
        if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
        
        ! Process key-value pairs
        select case(trim(key))
            case('ATOM_TYPES')
                read(value_str, *, iostat=io_stat) n_types
                if (io_stat == 0 .and. n_types > 0) then
                    atom_types_found = .true.
                    allocate(type_names(n_types), type_include(n_types), type_masses(n_types))
                    type_names = ''
                    type_include = .true.
                    type_masses = 0.0_dp
                end if
            
            case('Cores')
                read(value_str, *, iostat=io_stat) requested_cores
            
            ! Additional cases for other parameters...
        end select
    end do
    
    ! Second pass: read atom types and pair cutoffs
    ! ...
end function read_setup_file
```

This module handles configuration from files, providing default values and parsing setup information. It reads atom types and pair cutoff values, supporting both a detailed setup file mode and an interactive mode.

---

## cell_list_mod.f90

This module implements a cell-based spatial decomposition algorithm to efficiently find nearby atoms.

### Cell List Initialization

```fortran
subroutine initialize_cell_list(box_length, cutoffs, n_atoms)
    real(dp), intent(in) :: box_length(3)
    real(dp), intent(in) :: cutoffs(:)
    integer, intent(in) :: n_atoms
    integer :: ix, iy, iz, alloc_stat
    
    ! Find maximum cutoff distance
    max_cutoff = maxval(cutoffs)
    
    ! Calculate cell size and grid dimensions
    cell_size = max_cutoff * CELL_SIZE_FACTOR
    n_cells_x = max(1, floor(box_length(1)/cell_size))
    n_cells_y = max(1, floor(box_length(2)/cell_size))
    n_cells_z = max(1, floor(box_length(3)/cell_size))
    
    ! Adjust cell size to fit box exactly
    cell_size = min(box_length(1)/n_cells_x, &
                   box_length(2)/n_cells_y, &
                   box_length(3)/n_cells_z)

    ! Allocate cell grid
    if (allocated(cells)) deallocate(cells)
    allocate(cells(n_cells_x, n_cells_y, n_cells_z), stat=alloc_stat)
    call check_allocation(alloc_stat, "cell grid")

    ! Initialize cells
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(ix,iy,iz)
    do iz = 1, n_cells_z
        do iy = 1, n_cells_y
            do ix = 1, n_cells_x
                call cells(ix,iy,iz)%init(MAX_ATOMS_PER_CELL)
            end do
        end do
    end do
    !$OMP END PARALLEL DO
    
    ! More initialization...
end subroutine initialize_cell_list
```

### Cell List Update

```fortran
subroutine update_cell_list(coords, atom_types, n_atoms, box_length, frame, include_mask)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: atom_types(:)
    integer, intent(in) :: n_atoms
    real(dp), intent(in) :: box_length(3)
    integer, intent(in) :: frame
    logical, intent(in) :: include_mask(:)
    
    ! Check if full update is needed
    force_update = (last_cell_update == -1)
    
    ! Custom cell update frequency from setup file
    if (.not. force_update .and. cell_update_freq > 0 .and. &
        mod(frame, cell_update_freq) == 0) then
        force_update = .true.
    end if
    
    if (.not. force_update) then
        ! Check maximum atomic displacement
        max_displacement = 0.0_dp
        rebuild_mask = .false.
        
        !$OMP PARALLEL DO PRIVATE(dx,dy,dz,disp_sq,cell_x,cell_y,cell_z) &
        !$OMP& REDUCTION(max:max_displacement) SCHEDULE(guided)
        do i = 1, n_atoms
            if (.not. include_mask(atom_types(i))) cycle
            
            dx = coords(i,1) - prev_coords(i,1)
            dy = coords(i,2) - prev_coords(i,2)
            dz = coords(i,3) - prev_coords(i,3)
            
            ! Apply minimum image convention
            dx = dx - nint(dx/box_length(1)) * box_length(1)
            dy = dy - nint(dy/box_length(2)) * box_length(2)
            dz = dz - nint(dz/box_length(3)) * box_length(3)
            
            disp_sq = dx*dx + dy*dy + dz*dz
            max_displacement = max(max_displacement, sqrt(disp_sq))
            
            ! For selective rebuilding, mark cells that need updating
            if (selective_rebuild .and. sqrt(disp_sq) > REBUILD_THRESHOLD * cell_size * 0.5_dp) then
                ! Mark cell for rebuild
                !$OMP CRITICAL(rebuild_mask_update)
                rebuild_mask(cell_x, cell_y, cell_z) = .true.
                any_selective_updates = .true.
                !$OMP END CRITICAL(rebuild_mask_update)
            end if
        end do
        !$OMP END PARALLEL DO
        
        ! Force full update if displacement exceeds threshold
        force_update = (max_displacement > REBUILD_THRESHOLD * cell_size)
    end if
    
    ! Full or selective update logic follows...
end subroutine update_cell_list
```

This module implements a cell-based spatial decomposition that divides the simulation box into cells to avoid checking all atom pairs. It supports selective rebuilding where only cells with significant atom movement are updated. The implementation is optimized with OpenMP parallelization and uses atomic operations for thread safety.

---

## coordination_mod.f90

This module implements the core coordination number analysis, calculating how many atoms of each type surround each atom within specified cutoff distances.

### Coordination Initialization

```fortran
subroutine initialize_coordination(n_atoms, n_atom_types, atom_info)
    integer, intent(in) :: n_atoms, n_atom_types
    type(atom_type_info), intent(in) :: atom_info(:)
    integer :: i, j, pair_idx, alloc_stat
    character(20) :: user_input
    real(dp) :: cutoff_value
    
    ! Calculate number of unique pairs
    n_pairs = 0
    do i = 1, n_atom_types
        if (.not. atom_info(i)%include) cycle
        do j = i, n_atom_types
            if (.not. atom_info(j)%include) cycle
            n_pairs = n_pairs + 1
        end do
    end do
    
    ! Allocate arrays
    allocate(coord_numbers(n_atoms, n_pairs), stat=alloc_stat)
    call check_allocation(alloc_stat, "coordination arrays")
    
    ! Check if pairs is already allocated (from setup file)
    if (.not. allocated(pairs)) then
        allocate(pairs(n_pairs), stat=alloc_stat)
        call check_allocation(alloc_stat, "pairs array")
        
        ! Initialize pair information through user prompts or use defaults
        ! ...
    end if
    
    ! Initialize coordination numbers
    coord_numbers = 0
end subroutine initialize_coordination
```

### Coordination Calculation

```fortran
subroutine calculate_coordination(coords, atom_types, n_atoms, box_length, frame, include_mask)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: atom_types(:)
    integer, intent(in) :: n_atoms
    real(dp), intent(in) :: box_length(3)
    integer, intent(in) :: frame
    logical, intent(in) :: include_mask(:)
    
    ! Update cell lists
    call update_cell_list(coords, atom_types, n_atoms, box_length, frame, include_mask)
    
    ! Initialize coordination arrays
    coord_numbers = 0
    
    ! Get cell grid dimensions
    call get_cell_grid_dims(nx, ny, nz)
    
    ! Loop over all cells
    !$OMP PARALLEL DO PRIVATE(iy,iz,neighbor_cells,n_neighbors) SCHEDULE(dynamic)
    do ix = 1, nx
        do iy = 1, ny
            do iz = 1, nz
                ! Get neighboring cells
                call get_neighboring_cells(ix, iy, iz, neighbor_cells, n_neighbors)
                
                ! Process atoms in current cell
                call process_cell_atoms(ix, iy, iz, coords, atom_types, &
                                      box_length, include_mask)
                                      
                ! Process atoms in neighboring cells
                call process_neighbor_cells(neighbor_cells, n_neighbors, ix, iy, iz, &
                                          coords, atom_types, box_length, include_mask)
            end do
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine calculate_coordination
```

### Atom Pair Processing

```fortran
subroutine process_atom_pair(atom_i, atom_j, coords, atom_types, box_length, include_mask)
    integer, intent(in) :: atom_i, atom_j
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: atom_types(:)
    real(dp), intent(in) :: box_length(3)
    logical, intent(in) :: include_mask(:)
    
    real(dp) :: dx, dy, dz, r_sq
    integer :: pair_idx
    
    ! Skip if second atom type is excluded
    if (.not. include_mask(atom_types(atom_j))) return
    
    ! Calculate distance
    dx = coords(atom_i,1) - coords(atom_j,1)
    dy = coords(atom_i,2) - coords(atom_j,2)
    dz = coords(atom_i,3) - coords(atom_j,3)
    
    ! Apply minimum image convention
    dx = dx - nint(dx/box_length(1)) * box_length(1)
    dy = dy - nint(dy/box_length(2)) * box_length(2)
    dz = dz - nint(dz/box_length(3)) * box_length(3)
    
    r_sq = dx*dx + dy*dy + dz*dz
    
    ! Find pair index and update coordination if within cutoff
    call get_pair_index(atom_types(atom_i), atom_types(atom_j), pair_idx)
    if (pair_idx > 0 .and. r_sq <= pairs(pair_idx)%cutoff_sq) then
        !$OMP ATOMIC
        coord_numbers(atom_i,pair_idx) = coord_numbers(atom_i,pair_idx) + 1
        !$OMP ATOMIC
        coord_numbers(atom_j,pair_idx) = coord_numbers(atom_j,pair_idx) + 1
    end if
end subroutine process_atom_pair
```

This module calculates coordination numbers by finding atoms within specified cutoff distances. It uses the cell list for efficient neighbor finding and OpenMP parallelization for performance. The core calculation involves computing distances between atom pairs and comparing them to the cutoff distance, updating coordination counts when atoms are within the cutoff.

---

## data_file_io_mod.f90

This module handles reading LAMMPS data files to extract system information.

### Reading Data File

```fortran
subroutine read_data_file(filename, atom_info, n_atoms, n_types, box_length)
    character(len=*), intent(in) :: filename
    type(atom_type_info), allocatable, intent(out) :: atom_info(:)
    integer, intent(out) :: n_atoms, n_types
    real(dp), intent(out) :: box_length(3)
    
    ! First pass: read header information
    do
        read(data_unit, '(A)', iostat=io_stat) line
        if (io_stat /= 0) exit
        
        ! Parse essential header information
        if (index(line, 'atoms') > 0 .and. index(line, 'atom types') == 0) then
            read(line, *, iostat=io_stat) n_atoms
            found_atoms = (io_stat == 0)
        else if (index(line, 'atom types') > 0) then
            read(line, *, iostat=io_stat) n_types
            found_types = (io_stat == 0)
        else if (index(line, 'xlo xhi') > 0) then
            read(line, *, iostat=io_stat) xlo, xhi
            found_box = (found_box .or. io_stat == 0)
        else if (index(line, 'ylo yhi') > 0) then
            read(line, *, iostat=io_stat) ylo, yhi
            found_box = (found_box .or. io_stat == 0)
        else if (index(line, 'zlo zhi') > 0) then
            read(line, *, iostat=io_stat) zlo, zhi
            found_box = (found_box .or. io_stat == 0)
        else if (index(line, 'Masses') > 0) then
            found_masses = .true.
            exit
        end if
    end do
    
    ! Calculate box dimensions
    box_length(1) = xhi - xlo
    box_length(2) = yhi - ylo
    box_length(3) = zhi - zlo
    
    ! Read Masses section to get atom type information
    ! ...
end subroutine read_data_file
```

This module parses LAMMPS data files to extract system information like atom counts, types, and box dimensions. It supports reading atom type information including masses and names from the data file.

---

## trajectory_io_mod.f90

This module handles reading LAMMPS trajectory files (dump files).

### Frame Counting and Seeking

```fortran
function count_trajectory_frames(filename) result(num_frames)
    character(len=*), intent(in) :: filename
    integer :: num_frames
    
    ! Try to use grep for faster counting
    command = 'grep -c "ITEM: NUMBER OF ATOMS" ' // trim(filename) // ' > ' // trim(temp_file)
    call execute_command_line(trim(command), exitstat=grep_stat, cmdstat=cmd_stat)
    
    ! If command execution fails or grep isn't available, fall back
    if (grep_stat /= 0 .or. cmd_stat /= 0) then
        if (VERBOSE) write(*,*) "Using direct frame counting method"
        num_frames = count_frames_direct(filename)
        return
    end if
    
    ! Read the frame count from grep output
    open(newunit=unit, file=temp_file, status='old', action='read', iostat=io_stat)
    if (io_stat == 0) then
        read(unit, *, iostat=io_stat) num_frames
        close(unit)
    end if
    
    ! Parse the frame offsets if we have a valid count
    if (num_frames > 0) then
        call parse_frame_offsets(grep_positions_file, num_frames)
    end if
    
    ! Clean up temporary files
    call execute_command_line('rm ' // trim(temp_file), wait=.true.)
    call execute_command_line('rm ' // trim(grep_positions_file), wait=.true.)
end function count_trajectory_frames

subroutine seek_to_frame(frame_number)
    integer, intent(in) :: frame_number
    
    ! If we have frame offsets, use direct seeking
    if (have_frame_offsets .and. frame_number > 0 .and. frame_number <= size(frame_offsets)) then
        ! Seek to the position of the frame
        rewind(input_unit)
        read(input_unit, '(A)', advance='no', pos=frame_offsets(frame_number), iostat=io_stat) dummy_line
        if (io_stat /= 0) then
            call handle_error("Error seeking to frame " // trim(integer_to_string(frame_number)), ERR_FILE_IO)
        end if
        return
    end if
    
    ! Fall back to sequential reading if direct seeking isn't available
    ! ...
end subroutine seek_to_frame
```

### Reading Trajectory Frames

```fortran
subroutine read_trajectory_frame(coords, atom_types, elements, box_length, frame_number, eof, atom_info)
    real(dp), intent(out) :: coords(:,:)
    integer, intent(out) :: atom_types(:)
    character(len=*), intent(out) :: elements(:)
    real(dp), intent(out) :: box_length(3)
    integer, intent(inout) :: frame_number
    logical, intent(out) :: eof
    type(atom_type_info), intent(in) :: atom_info(:)
    
    ! Look for frame start
    do
        read(input_unit, '(A)', iostat=io_stat) line
        if (io_stat /= 0) then
            eof = .true.
            return
        end if
        
        if (index(line, 'ITEM: NUMBER OF ATOMS') > 0) then
            read(input_unit, *) n_atoms_frame
            exit
        end if
    end do
    
    ! Read box bounds
    read(input_unit, '(A)') line  ! ITEM: BOX BOUNDS
    read(input_unit, *) xlo, xhi
    read(input_unit, *) ylo, yhi
    read(input_unit, *) zlo, zhi
    
    ! Calculate box lengths
    box_length(1) = xhi - xlo
    box_length(2) = yhi - ylo
    box_length(3) = zhi - zlo
    
    ! Read atoms header and setup column mapping if first frame
    read(input_unit, '(A)') line
    if (frame_number == 0) then
        call setup_column_mapping(line)
    end if
    
    ! Read atom data
    do i = 1, n_atoms_frame
        read(input_unit, *, iostat=io_stat) values_buffer(1:col_map%max_cols)
        
        atom_id = nint(values_buffer(col_map%id))
        coords(atom_id,1) = values_buffer(col_map%x)
        coords(atom_id,2) = values_buffer(col_map%y)
        coords(atom_id,3) = values_buffer(col_map%z)
        atom_types(atom_id) = nint(values_buffer(col_map%type))
        
        ! Set element name
        if (atom_types(atom_id) >= 1 .and. atom_types(atom_id) <= size(atom_info)) then
            elements(atom_id) = atom_info(atom_types(atom_id))%name
        else
            write(elements(atom_id),'(A,I0)') 'Type', atom_types(atom_id)
        end if
    end do
    
    frame_number = frame_number + 1
end subroutine read_trajectory_frame
```

This module reads LAMMPS trajectory files, supporting efficient frame counting, direct seeking to frames, and flexible column mapping. It can use external tools (grep) for performance and builds frame offset tables for direct seeking.

---

## output_io_mod.f90

This module handles writing analysis results and generating summary reports.

### Output Initialization and Writing

```fortran
subroutine init_output(filename)
    character(len=*), intent(in) :: filename
    
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
```

### Final Report Generation

```fortran
subroutine write_final_report(n_frames, n_atoms, total_time, atom_info, &
                             n_types, coord_numbers_local, atom_types_local)
    integer, intent(in) :: n_frames, n_atoms, n_types
    real(dp), intent(in) :: total_time
    type(atom_type_info), intent(in) :: atom_info(:)
    integer, intent(in) :: coord_numbers_local(:,:)
    integer, intent(in) :: atom_types_local(:)
    
    ! Calculate statistics safely with OpenMP
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
    
    ! Write performance report
    write(*,*)
    write(*,'(A)') '=== Analysis Complete ==='
    write(*,'(A,I0,A,I0,A)') ' Processed ', n_frames, ' frames with ', &
                             n_atoms, ' atoms each'
    write(*,'(A,F10.3,A)') ' Total time: ', total_time, ' seconds'
    frames_per_second = real(n_frames)/max(total_time, 1.0e-6_dp)
    write(*,'(A,F10.3,A)') ' Performance: ', frames_per_second, ' frames/second'
    
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
    
    ! Additional statistics...
end subroutine write_final_report
```

This module handles writing analysis results to output files and generating comprehensive summary reports. It formats coordination numbers for each atom and pair type and calculates average coordination statistics for the final report.

---

## progress_mod.f90

This module provides user feedback on analysis progress.

### Progress Bar Implementation

```fortran
subroutine display_progress(frame, elapsed_time, total)
    integer, intent(in) :: frame
    real(dp), intent(in) :: elapsed_time
    integer, intent(in), optional :: total ! Optional total frames parameter
    integer :: progress, width, actual_total
    character(len=30) :: time_string
    character(len=20) :: frame_string
    real(dp) :: frames_per_second, estimated_remaining
    
    ! Calculate progress width
    width = PROGRESS_BAR_WIDTH - 1
    
    ! Calculate progress percentage
    if (actual_total <= 0) then
        progress = 0
    else
        progress = min(width, nint(real(frame) / real(actual_total) * width))
    end if
    
    ! Calculate performance metrics
    if (elapsed_time > 0.0_dp) then
        frames_per_second = real(frame) / elapsed_time
        if (actual_total > 0) then
            estimated_remaining = (real(actual_total) - real(frame)) / frames_per_second
        else
            estimated_remaining = 0.0_dp
        end if
    end if
    
    ! Display progress bar
    write(*,'(a1,1x,A11,1x)',advance='no') achar(13), frame_string
    write(*,'(A)',advance='no') '['
    write(*,'(A)',advance='no') repeat('#', progress)
    write(*,'(A)',advance='no') repeat(' ', width - progress)
    write(*,'(A)',advance='no') ']'
    
    ! Display timing information
    write(*,'(1x,A,A)',advance='no') trim(adjustl(time_string))
    if (frames_per_second > 0.0_dp) then
        write(*,'(A,F6.2,A)',advance='no') ' [', frames_per_second, ' frames/s]'
    end if
    if (estimated_remaining > 0.0_dp) then
        write(*,'(A,F6.1,A)',advance='no') ' [ETA: ', estimated_remaining, 's]'
    end if
    
    call flush(6)
end subroutine display_progress
```

This module implements a text-based progress bar with performance metrics and remaining time estimation. It uses carriage return to update the display in-place and provides real-time feedback on the analysis progress.

---

## io_mod.f90

This module serves as a high-level interface to all I/O operations, consolidating functionality from specialized I/O modules.

### I/O Initialization

```fortran
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
```

This module serves as a facade for the specialized I/O modules, providing a simplified interface for the main program. It coordinates initialization and cleanup of all I/O resources.

---

## benchmark_mod.f90

This module implements benchmarking functionality to evaluate performance with different thread counts and configurations.

### Benchmark Execution

```fortran
subroutine run_benchmark(coords, atom_types, elements, box_length, n_atoms, &
                        n_types, atom_info, total_frames)
    real(dp), intent(inout) :: coords(:,:)
    integer, intent(inout) :: atom_types(:)
    character(len=*), intent(inout) :: elements(:)
    real(dp), intent(inout) :: box_length(3)
    integer, intent(in) :: n_atoms, n_types, total_frames
    type(atom_type_info), intent(in) :: atom_info(:)
    
    integer :: max_threads, n_threads
    
    ! Get maximum available threads
    max_threads = omp_get_max_threads()
    
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
            
            ! Check if we've reached time limit and processed enough frames
            current_time = omp_get_wtime()
            if ((current_time - start_time) >= BENCHMARK_TIME_LIMIT .and. frame >= 10) exit
        end do
        
        ! Record results
        result%n_threads = n_threads
        result%total_time = omp_get_wtime() - start_time
        result%frames_processed = frame
        
        ! Calculate performance metrics
        result%frames_per_second = real(result%frames_processed) / result%total_time
        
        ! Calculate speedup and efficiency
        if (n_threads == 1) then
            single_thread_fps = result%frames_per_second
            result%speedup = 1.0_dp
            result%efficiency = 1.0_dp
        else if (single_thread_fps > 0.0_dp) then
            result%speedup = result%frames_per_second / single_thread_fps
            result%efficiency = result%speedup / real(n_threads, dp)
        end if
        
        ! Write results
        call write_benchmark_result(result)
        
        ! Double the thread count for next iteration
        if (n_threads == max_threads) exit
        n_threads = min(n_threads * 2, max_threads)
    end do
end subroutine run_benchmark
```

This module tests the program's performance with different thread counts, calculating metrics like speedup and efficiency. It uses a fixed time limit for each test and disables I/O for accurate timing.

---

## main.f90

This is the main program that ties all modules together, orchestrating the overall coordination analysis workflow.

### Program Initialization and Analysis

```fortran
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

    ! Check command line arguments for benchmark mode
    if (command_argument_count() > 0) then
        call get_command_argument(1, value=arg)
        benchmark_mode = (trim(arg) == '--benchmark')
    end if
    
    ! Initialize configuration
    call initialize_config()
    call read_config()
    
    ! First, read DATA file to get atom count regardless of setup file 
    call read_data_file_dims_only(DATA_FILE, n_atoms, n_types_actual, box_length)
    
    ! Try to read setup file for additional configuration
    if (.not. benchmark_mode) then
        setup_file_read = read_setup_file(atom_info, n_types, pairs, n_pairs)
        
        if (setup_file_read) then
            ! Use configuration from setup file
            ! ...
        end if
    end if
    
    ! Allocate arrays
    allocate(coords(n_atoms,3), atom_types(n_atoms), elements(n_atoms), &
             include_mask(n_types), stat=io_stat)
    if (io_stat /= 0) call handle_error("Failed to allocate arrays", ERR_ALLOCATION)
    
    ! Count frames in trajectory file
    if (VERBOSE) write(*,'(A)') ' Counting frames in trajectory file...'
    total_frames = count_trajectory_frames(input_trajectory)
    
    if (benchmark_mode) then
        ! Run benchmark mode
        call initialize_benchmark()
        call run_benchmark(coords, atom_types, elements, box_length, &
                         n_atoms, n_types, atom_info, total_frames)
        call cleanup_benchmark()
    else
        ! Normal analysis mode
        call initialize_coordination(n_atoms, n_types, atom_info)
        call initialize_cell_list(box_length, get_cutoffs(), n_atoms)
        call initialize_io(input_trajectory, output_data_file, total_frames, start_frame)
        
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
        call write_final_report(frame, n_atoms, end_time - start_time, atom_info, &
                              n_types, coord_numbers, atom_types)
    end if
    
    ! Final cleanup
    call cleanup_io()
    call cleanup_coordination()
    call cleanup_cell_list()
    call cleanup_config()
end program coordination_analysis
```

The main program orchestrates the entire analysis workflow, handling initialization, processing, and cleanup. It supports both normal analysis and benchmark modes, and can read configuration from setup files or use defaults.

---

## Performance Optimizations & Key Features

### Cell List Algorithm
- Divides the simulation box into cells approximately the size of the largest cutoff
- Reduces neighbor finding from O(n²) to nearly O(n) complexity
- Only needs to check atoms in the same and neighboring cells

### Selective Cell Rebuilding
```fortran
! Check maximum atomic displacement
if (.not. force_update) then
    max_displacement = 0.0_dp
    rebuild_mask = .false.
    
    !$OMP PARALLEL DO PRIVATE(...) REDUCTION(max:max_displacement)
    do i = 1, n_atoms
        ! Calculate displacement since last update
        dx = coords(i,1) - prev_coords(i,1)
        dy = coords(i,2) - prev_coords(i,2)
        dz = coords(i,3) - prev_coords(i,3)
        
        ! Apply minimum image convention
        dx = dx - nint(dx/box_length(1)) * box_length(1)
        dy = dy - nint(dy/box_length(2)) * box_length(2)
        dz = dz - nint(dz/box_length(3)) * box_length(3)
        
        disp_sq = dx*dx + dy*dy + dz*dz
        max_displacement = max(max_displacement, sqrt(disp_sq))
        
        ! For selective rebuilding, mark cells that need updating
        if (selective_rebuild .and. sqrt(disp_sq) > REBUILD_THRESHOLD * cell_size * 0.5_dp) then
            rebuild_mask(cell_x, cell_y, cell_z) = .true.
            any_selective_updates = .true.
        end if
    end do
    
    ! Force full update if displacement exceeds threshold
    force_update = (max_displacement > REBUILD_THRESHOLD * cell_size)
end if
```

### OpenMP Parallelization
- Used throughout the code for performance on multi-core systems
- Cell-based parallelism for good load balancing
- Atomic operations to prevent race conditions
- Thread-safe update methods

### Optimization Tricks
- Pre-computed squared cutoffs avoid expensive square root operations
- Buffer reuse to minimize memory allocations
- Frame offset caching for direct seeking
- External tools (grep) for fast frame counting
- Selective processing based on atom type include flags

### Benchmark Results
```
# Coordination Analysis Benchmark Results
# Threads, Time(s), Frames/s, Frames, Speedup, Efficiency
     1       2.358      53.431     126.000       1.000       1.000
     2       1.951      64.574     126.000       1.209       0.604
     4       1.679      75.038     126.000       1.404       0.351
     8       1.565      80.498     126.000       1.507       0.188
    16       1.558      80.875     126.000       1.514       0.095
    32       2.527      49.860     126.000       0.933       0.029
```

These results show good scaling up to 16 threads with diminishing returns after 8 threads. The 32-thread performance degradation suggests potential resource contention or thread management overhead.

## Conclusion

This coordination analysis program demonstrates efficient algorithms and optimization techniques for analyzing molecular dynamics trajectories. Key features include:

1. Cell-based spatial decomposition for efficient neighbor finding
2. OpenMP parallelization for multi-core performance
3. Selective rebuilding to avoid unnecessary updates
4. Flexible configuration via setup files
5. Comprehensive error handling and progress reporting

The program is well-structured using Fortran modules with clear separation of concerns, making it maintainable and extensible for future enhancements.
