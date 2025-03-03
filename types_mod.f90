module types_mod
    use iso_fortran_env, only: dp => real64
    implicit none
    private
    
    ! Public type definitions
    public :: atom_type_info, pair_type, cell_type
    public :: dp

    !===============================================================================
    ! Type definitions
    !===============================================================================
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

contains
    
    subroutine initialize_cell(this, initial_size)
        class(cell_type), intent(inout) :: this
        integer, intent(in) :: initial_size
        integer :: alloc_stat
        
        if (allocated(this%atoms)) deallocate(this%atoms)
        allocate(this%atoms(initial_size), stat=alloc_stat)
        if (alloc_stat /= 0) error stop "Failed to allocate cell array"
        this%n_atoms = 0
    end subroutine initialize_cell
    
    subroutine cleanup_cell(this)
        class(cell_type), intent(inout) :: this
        if (allocated(this%atoms)) deallocate(this%atoms)
        this%n_atoms = 0
    end subroutine cleanup_cell
    
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

end module types_mod