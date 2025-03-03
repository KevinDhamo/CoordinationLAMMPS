module coordination_mod
    use types_mod
    use config_mod
    use cell_list_mod
    use error_mod
    use omp_lib
    implicit none
    private

    ! Public procedures
    public :: initialize_coordination
    public :: calculate_coordination
    public :: cleanup_coordination

    ! Public variables
    public :: coord_numbers, pairs, n_pairs

    ! Module variables
    integer, allocatable :: coord_numbers(:,:)
    type(pair_type), allocatable :: pairs(:)
    integer :: n_pairs

contains
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
            
            ! Initialize pair information through user prompts
            write(*,'(A)') ' Setting up pair cutoff distances:'
            pair_idx = 1
            do i = 1, n_atom_types
                if (.not. atom_info(i)%include) cycle
                do j = i, n_atom_types
                    if (.not. atom_info(j)%include) cycle
                    
                    pairs(pair_idx)%type1_id = i
                    pairs(pair_idx)%type2_id = j
                    
                    write(*,'(A,I0,A,A,A,I0,A,A,A,F5.2,A)', advance='no') &
                        '   Enter cutoff for type ', i, ' (', trim(atom_info(i)%name), ') - type ', j, &
                        ' (', trim(atom_info(j)%name), ') (Å) [', DEFAULT_CUTOFF, ']: '
                    
                    read(*,'(A)') user_input
                    if (len_trim(user_input) == 0) then
                        cutoff_value = DEFAULT_CUTOFF
                    else
                        read(user_input,*,iostat=alloc_stat) cutoff_value
                        if (alloc_stat /= 0) then
                            call handle_error("Invalid cutoff value", ERR_INVALID_PARAM)
                        end if
                    end if
                    
                    ! Validate cutoff value
                    if (cutoff_value <= 0.0_dp) then
                        call handle_error("Cutoff value must be positive", ERR_INVALID_PARAM)
                    end if
                    
                    pairs(pair_idx)%cutoff = cutoff_value
                    pairs(pair_idx)%cutoff_sq = cutoff_value**2
                    pair_idx = pair_idx + 1
                end do
            end do
        else
            ! Pairs already allocated (from setup file)
            ! Just print out the cutoff values for reference
            write(*,'(A)') ' Using pair cutoffs from setup file:'
            do pair_idx = 1, n_pairs
                i = pairs(pair_idx)%type1_id
                j = pairs(pair_idx)%type2_id
                write(*,'(A,I0,A,A,A,I0,A,A,A,F5.2)') &
                    '   Type ', i, ' (', trim(atom_info(i)%name), ') - type ', j, &
                    ' (', trim(atom_info(j)%name), '): ', pairs(pair_idx)%cutoff
            end do
        end if
        
        ! Initialize coordination numbers
        coord_numbers = 0
    end subroutine initialize_coordination

    subroutine calculate_coordination(coords, atom_types, n_atoms, box_length, frame, include_mask)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        integer, intent(in) :: n_atoms
        real(dp), intent(in) :: box_length(3)
        integer, intent(in) :: frame
        logical, intent(in) :: include_mask(:)
        
        integer :: ix, iy, iz, nx, ny, nz
        integer :: neighbor_cells(27,3), n_neighbors
        
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

    subroutine process_cell_atoms(ix, iy, iz, coords, atom_types, box_length, include_mask)
        integer, intent(in) :: ix, iy, iz
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        real(dp), intent(in) :: box_length(3)
        logical, intent(in) :: include_mask(:)
        
        integer :: i, j, atom_i, atom_j
        type(cell_type) :: current_cell
        
        ! Get current cell
        call get_cell(ix, iy, iz, current_cell)
        
        ! Process atoms within the cell
        do i = 1, current_cell%n_atoms
            atom_i = current_cell%atoms(i)
            if (.not. include_mask(atom_types(atom_i))) cycle
            
            do j = i+1, current_cell%n_atoms
                atom_j = current_cell%atoms(j)
                call process_atom_pair(atom_i, atom_j, coords, atom_types, &
                                     box_length, include_mask)
            end do
        end do
    end subroutine process_cell_atoms

    subroutine process_neighbor_cells(neighbor_cells, n_neighbors, ix, iy, iz, &
                                    coords, atom_types, box_length, include_mask)
        integer, intent(in) :: neighbor_cells(:,:), n_neighbors, ix, iy, iz
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        real(dp), intent(in) :: box_length(3)
        logical, intent(in) :: include_mask(:)
        
        integer :: k
        
        do k = 1, n_neighbors
            call process_cell_pair(ix, iy, iz, &
                                 neighbor_cells(k,1), neighbor_cells(k,2), neighbor_cells(k,3), &
                                 coords, atom_types, box_length, include_mask)
        end do
    end subroutine process_neighbor_cells

    subroutine process_cell_pair(ix1, iy1, iz1, ix2, iy2, iz2, &
                                coords, atom_types, box_length, include_mask)
        integer, intent(in) :: ix1, iy1, iz1, ix2, iy2, iz2
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        real(dp), intent(in) :: box_length(3)
        logical, intent(in) :: include_mask(:)
        
        integer :: i, j, atom_i, atom_j
        type(cell_type) :: cell1, cell2
        
        ! Skip if same cell
        if (ix1 == ix2 .and. iy1 == iy2 .and. iz1 == iz2) return
        
        ! Get cells
        call get_cell(ix1, iy1, iz1, cell1)
        call get_cell(ix2, iy2, iz2, cell2)
        
        ! Process atoms between cells
        do i = 1, cell1%n_atoms
            atom_i = cell1%atoms(i)
            if (.not. include_mask(atom_types(atom_i))) cycle
            
            do j = 1, cell2%n_atoms
                atom_j = cell2%atoms(j)
                call process_atom_pair(atom_i, atom_j, coords, atom_types, &
                                     box_length, include_mask)
            end do
        end do
    end subroutine process_cell_pair

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

    subroutine get_pair_index(type1, type2, idx)
        integer, intent(in) :: type1, type2
        integer, intent(out) :: idx
        integer :: t1, t2
        
        ! Order types for consistent pair lookup
        t1 = min(type1, type2)
        t2 = max(type1, type2)
        
        do idx = 1, n_pairs
            if (pairs(idx)%type1_id == t1 .and. pairs(idx)%type2_id == t2) return
        end do
        
        ! Pair not found
        idx = -1
    end subroutine get_pair_index

    subroutine cleanup_coordination()
        if (allocated(coord_numbers)) deallocate(coord_numbers)
        if (allocated(pairs)) deallocate(pairs)
    end subroutine cleanup_coordination

end module coordination_mod