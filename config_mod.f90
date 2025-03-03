module config_mod
    use types_mod
    use error_mod
    implicit none
    private

    ! Public parameters and procedures
    public :: initialize_config
    public :: read_config
    public :: read_setup_file
    public :: get_config_value
    public :: cleanup_config

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

    ! Setup file parameters
    logical, public :: setup_file_found = .false.
    integer, public :: requested_cores = 0
    integer, public :: start_frame = 0       ! Start processing from this frame
    integer, public :: end_frame = 0         ! End processing at this frame (0 = all frames)
    logical, public :: selective_rebuild = .false. ! Enable selective cell rebuilding
    character(len=256), public :: input_data_file = ''
    character(len=256), public :: input_trajectory = ''
    character(len=256), public :: output_data_file = ''

    ! Type definition for configuration values
    type :: config_value_type
        character(len=32) :: key = ''
        character(len=128) :: value = ''
    end type config_value_type

    ! Module variables
    type(config_value_type), allocatable :: config_data(:)

contains

    subroutine initialize_config()
        integer :: n_config_items = 20  ! Initial size
        integer :: alloc_stat
        
        if (allocated(config_data)) deallocate(config_data)
        allocate(config_data(n_config_items), stat=alloc_stat)
        call check_allocation(alloc_stat, "config_data")
        
        ! Initialize all entries to empty strings
        config_data%key = ''
        config_data%value = ''
        
        ! Set defaults for setup parameters
        setup_file_found = .false.
        requested_cores = 0
        start_frame = 0
        end_frame = 0
        cell_update_freq = 5  ! Default value of 5, can be overridden in setup.txt
        selective_rebuild = .false.
        input_data_file = DATA_FILE
        input_trajectory = INPUT_FILE
        output_data_file = OUTPUT_FILE
    end subroutine initialize_config

    subroutine read_config()
        integer :: io_stat, unit
        character(len=256) :: line
        character(len=32) :: key
        character(len=128) :: value
        logical :: file_exists
        
        inquire(file=CONFIG_FILE, exist=file_exists)
        if (.not. file_exists) return
        
        open(newunit=unit, file=CONFIG_FILE, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to open configuration file", ERR_FILE_IO)
            return
        end if
        
        do
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            
            ! Skip comments and empty lines
            if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
            
            ! Parse key-value pair
            call parse_config_line(line, key, value)
            if (len_trim(key) > 0) then
                call store_config_value(key, value)
            end if
        end do
        
        close(unit)
    end subroutine read_config
    
    ! Fixed function to parse pair cutoffs properly
    subroutine parse_cutoff_key(key_str, type1, type2, success)
        character(len=*), intent(in) :: key_str
        integer, intent(out) :: type1, type2
        logical, intent(out) :: success
        
        integer :: underscore_pos, io_stat
        character(len=32) :: prefix, type1_str, type2_str
        
        success = .false.
        type1 = -1
        type2 = -1
        
        ! Check if it's a cutoff key
        prefix = "PAIR_CUTOFF_"
        if (len_trim(key_str) <= len_trim(prefix)) return
        if (key_str(1:len_trim(prefix)) /= prefix) return
        
        ! Extract the part after the prefix
        type1_str = key_str(len_trim(prefix)+1:)
        
        ! Find the underscore
        underscore_pos = index(type1_str, "_")
        if (underscore_pos <= 1) return
        
        ! Split at the underscore
        type2_str = type1_str(underscore_pos+1:)
        type1_str = type1_str(1:underscore_pos-1)
        
        ! Convert to integers
        read(type1_str, *, iostat=io_stat) type1
        if (io_stat /= 0) return
        
        read(type2_str, *, iostat=io_stat) type2
        if (io_stat /= 0) return
        
        success = .true.
    end subroutine parse_cutoff_key
 
    function read_setup_file(atom_info, n_types, pairs, n_pairs) result(success)
        type(atom_type_info), allocatable, intent(out) :: atom_info(:)
        integer, intent(out) :: n_types
        type(pair_type), allocatable, intent(out) :: pairs(:)
        integer, intent(out) :: n_pairs
        logical :: success
        
        integer :: io_stat, unit, i, j, pair_idx
        character(len=256) :: line, value_str
        character(len=32) :: key
        logical :: file_exists
        logical :: atom_types_found = .false.
        character(len=8), allocatable :: type_names(:)
        logical, allocatable :: type_include(:)
        real(dp), allocatable :: type_masses(:)
        character(len=256) :: comment
        integer :: comment_pos
        logical :: success_parse
        
        success = .false.
        
        ! Check if setup file exists
        inquire(file=SETUP_FILE, exist=file_exists)
        if (.not. file_exists) then
            if (VERBOSE) write(*,*) "Setup file not found, using interactive mode"
            return
        end if
        
        ! Open setup file
        open(newunit=unit, file=SETUP_FILE, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to open setup file", ERR_FILE_IO, fatal=.false.)
            return
        end if
        
        ! First pass: count atom types and find basic parameters
        n_types = 0
        do
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            
            ! Skip comments and empty lines
            if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
            
            ! Parse line
            key = adjustl(line)
            value_str = ''
            if (index(line, '=') > 0) then
                key = adjustl(line(:index(line, '=')-1))
                value_str = adjustl(line(index(line, '=')+1:))
                
                ! Check for comments
                comment_pos = index(value_str, '#')
                if (comment_pos > 0) then
                    comment = adjustl(value_str(comment_pos+1:))
                    value_str = adjustl(value_str(:comment_pos-1))
                end if
            end if
            
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
                
                case('START_FRAME')
                    read(value_str, *, iostat=io_stat) start_frame
                    if (io_stat /= 0) start_frame = 0
                
                case('END_FRAME')
                    if (trim(value_str) == 'unlimited' .or. &
                        trim(value_str) == 'all') then
                        end_frame = 0  ! Indicating unlimited
                    else
                        read(value_str, *, iostat=io_stat) end_frame
                        if (io_stat /= 0) end_frame = 0
                    end if
                    
                case('CELL_UPDATE_FREQ')
                    read(value_str, *, iostat=io_stat) cell_update_freq
                    if (io_stat /= 0) cell_update_freq = 5  ! Default to 5 if invalid input
                
                case('SELECTIVE_REBUILD')
                    if (trim(value_str) == 'yes' .or. &
                        trim(value_str) == 'true' .or. &
                        trim(value_str) == '1') then
                        selective_rebuild = .true.
                    else
                        selective_rebuild = .false.
                    end if
                
                case('INPUT_DATA')
                    input_data_file = trim(value_str)
                
                case('INPUT_TRAJECTORY')
                    input_trajectory = trim(value_str)
                
                case('OUTPUT_FILE')
                    output_data_file = trim(value_str)
            end select
        end do
        
        ! Verify that we found atom types
        if (.not. atom_types_found .or. n_types <= 0) then
            if (VERBOSE) write(*,*) "No valid ATOM_TYPES found in setup file"
            close(unit)
            return
        end if
        
        ! Second pass: read atom types and pair cutoffs
        rewind(unit)
        do
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            
            ! Skip comments and empty lines
            if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
            
            ! Parse line
            key = adjustl(line)
            value_str = ''
            if (index(line, '=') > 0) then
                key = adjustl(line(:index(line, '=')-1))
                value_str = adjustl(line(index(line, '=')+1:))
                
                ! Check for ignore comment
                comment_pos = index(value_str, '#')
                if (comment_pos > 0) then
                    comment = adjustl(value_str(comment_pos+1:))
                    value_str = adjustl(value_str(:comment_pos-1))
                    
                    ! Check if it's an ignore directive
                    if (index(comment, 'ignore') > 0) then
                        ! Extract the type number from TYPE_X
                        if (key(1:5) == 'TYPE_') then
                            read(key(6:), *, iostat=io_stat) i
                            if (io_stat == 0 .and. i > 0 .and. i <= n_types) then
                                type_include(i) = .false.
                                if (VERBOSE) write(*,*) "Ignoring atom type:", i
                            end if
                        end if
                    end if
                end if
            end if
            
            ! Process atom types
            if (key(1:5) == 'TYPE_') then
                read(key(6:), *, iostat=io_stat) i
                if (io_stat == 0 .and. i > 0 .and. i <= n_types) then
                    type_names(i) = trim(adjustl(value_str))
                end if
            end if
        end do
        
        ! Calculate number of unique pairs
        n_pairs = 0
        do i = 1, n_types
            if (.not. type_include(i)) cycle
            do j = i, n_types
                if (.not. type_include(j)) cycle
                n_pairs = n_pairs + 1
            end do
        end do
        
        ! Allocate array for atom types
        if (allocated(atom_info)) deallocate(atom_info)
        allocate(atom_info(n_types), stat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to allocate atom_info array", ERR_ALLOCATION)
            return
        end if
        
        ! Initialize atom info
        do i = 1, n_types
            atom_info(i)%type_id = i
            atom_info(i)%name = type_names(i)
            atom_info(i)%mass = type_masses(i)
            atom_info(i)%include = type_include(i)
        end do
        
        ! Allocate pairs array
        if (allocated(pairs)) deallocate(pairs)
        allocate(pairs(n_pairs), stat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to allocate pairs array", ERR_ALLOCATION)
            return
        end if
        
        ! Initialize pair indices first
        pair_idx = 1
        do i = 1, n_types
            if (.not. type_include(i)) cycle
            do j = i, n_types
                if (.not. type_include(j)) cycle
                
                ! Set default cutoff
                pairs(pair_idx)%type1_id = i
                pairs(pair_idx)%type2_id = j
                pairs(pair_idx)%cutoff = DEFAULT_CUTOFF
                pairs(pair_idx)%cutoff_sq = DEFAULT_CUTOFF**2
                
                pair_idx = pair_idx + 1
            end do
        end do
        
        ! Read cutoff values
        rewind(unit)
        do
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            
            ! Skip comments and empty lines
            if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
            
            ! Parse line
            if (index(line, '=') > 0) then
                key = adjustl(line(:index(line, '=')-1))
                value_str = adjustl(line(index(line, '=')+1:))
                
                ! Check for PAIR_CUTOFF entries
                if (key(1:12) == 'PAIR_CUTOFF_') then
                    ! Parse the cutoff key to extract type indices
                    call parse_cutoff_key(key, i, j, success_parse)
                    
                    if (success_parse .and. i > 0 .and. i <= n_types .and. j > 0 .and. j <= n_types) then
                        ! Find the pair index
                        do pair_idx = 1, n_pairs
                            ! Ensure i <= j for consistent lookup
                            if ((pairs(pair_idx)%type1_id == min(i,j)) .and. &
                                (pairs(pair_idx)%type2_id == max(i,j))) then
                                
                                read(value_str, *, iostat=io_stat) pairs(pair_idx)%cutoff
                                if (io_stat == 0) then
                                    pairs(pair_idx)%cutoff_sq = pairs(pair_idx)%cutoff**2
                                end if
                                exit
                            end if
                        end do
                    end if
                end if
            end if
        end do
        
        close(unit)
        
        ! Print summary if verbose
        if (VERBOSE) then
            write(*,'(A)') " Setup file successfully read:"
            write(*,'(A,I0)') "   Atom types: ", n_types
            write(*,'(A,I0)') "   Atom pairs: ", n_pairs
            write(*,'(A)') "   Frame selection:"
            if (start_frame > 0) then
                write(*,'(A,I0)') "     Start frame: ", start_frame
            else
                write(*,'(A)') "     Start frame: beginning"
            end if
            if (end_frame > 0) then
                write(*,'(A,I0)') "     End frame: ", end_frame
            else
                write(*,'(A)') "     End frame: end of trajectory"
            end if
            write(*,'(A,I0)') "   Cell update frequency: ", cell_update_freq
            write(*,'(A,L1)') "   Selective cell rebuilding: ", selective_rebuild
            write(*,'(A)') "   Atom type details:"
            do i = 1, n_types
                write(*,'(A,I0,A,A,A,L1)') "     Type ", i, ": ", trim(atom_info(i)%name), &
                                          ", Include: ", atom_info(i)%include
            end do
            write(*,'(A)') "   Pair cutoffs:"
            do i = 1, n_pairs
                write(*,'(A,I0,A,A,A,I0,A,A,A,F6.3)') "     Type ", pairs(i)%type1_id, &
                                       " (", trim(atom_info(pairs(i)%type1_id)%name), ") - Type ", pairs(i)%type2_id, &
                                       " (", trim(atom_info(pairs(i)%type2_id)%name), "): ", pairs(i)%cutoff
            end do
        end if
        
        ! Cleanup temporary arrays
        if (allocated(type_names)) deallocate(type_names)
        if (allocated(type_include)) deallocate(type_include)
        if (allocated(type_masses)) deallocate(type_masses)
        
        setup_file_found = .true.
        success = .true.
    end function read_setup_file

    subroutine parse_config_line(line, key, value)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: key
        character(len=*), intent(out) :: value
        integer :: eq_pos
        
        eq_pos = index(line, '=')
        if (eq_pos > 0) then
            key = adjustl(line(:eq_pos-1))
            value = adjustl(line(eq_pos+1:))
        else
            key = ''
            value = ''
        end if
    end subroutine parse_config_line

    subroutine store_config_value(key, value)
        character(len=*), intent(in) :: key
        character(len=*), intent(in) :: value
        integer :: i, empty_slot
        
        ! Find existing key or first empty slot
        empty_slot = -1
        do i = 1, size(config_data)
            if (config_data(i)%key == key) then
                config_data(i)%value = value
                return
            else if (empty_slot == -1 .and. len_trim(config_data(i)%key) == 0) then
                empty_slot = i
            end if
        end do
        
        ! Store in empty slot if found
        if (empty_slot > 0) then
            config_data(empty_slot)%key = key
            config_data(empty_slot)%value = value
        else
            call handle_error("Configuration storage full", ERR_INVALID_PARAM, fatal=.false.)
        end if
    end subroutine store_config_value

    function get_config_value(key, default) result(value)
        character(len=*), intent(in) :: key
        character(len=*), intent(in) :: default
        character(len=128) :: value
        integer :: i
        
        value = default
        if (.not. allocated(config_data)) return
        
        do i = 1, size(config_data)
            if (config_data(i)%key == key) then
                value = config_data(i)%value
                return
            end if
        end do
    end function get_config_value

    subroutine cleanup_config()
        if (allocated(config_data)) deallocate(config_data)
    end subroutine cleanup_config

end module config_mod
