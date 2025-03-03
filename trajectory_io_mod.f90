module trajectory_io_mod
    use types_mod
    use config_mod
    use error_mod
    implicit none
    private

    ! Public procedures
    public :: init_trajectory_reader
    public :: read_trajectory_frame
    public :: count_trajectory_frames
    public :: seek_to_frame
    public :: cleanup_trajectory_reader

    ! Module variables
    integer :: input_unit = -1
    logical :: input_file_open = .false.
    
    ! Column mapping type
    type :: column_map_type
        integer :: id = -1
        integer :: type = -1
        integer :: x = -1
        integer :: y = -1
        integer :: z = -1
        integer :: max_cols = -1
    end type
    
    type(column_map_type), save :: col_map
    
    ! Frame location cache for fast seeking
    integer, allocatable :: frame_offsets(:)
    logical :: have_frame_offsets = .false.
    
    ! Module-level buffer for values to avoid allocation in read_trajectory_frame
    real(dp), allocatable :: values_buffer(:)

contains
    subroutine init_trajectory_reader(filename)
        character(len=*), intent(in) :: filename
        integer :: io_stat
        
        ! Open trajectory file
        open(newunit=input_unit, file=filename, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to open trajectory file: " // trim(filename), ERR_FILE_IO)
        end if
        input_file_open = .true.
        
        ! Allocate module-level buffer
        if (allocated(values_buffer)) deallocate(values_buffer)
        allocate(values_buffer(1000))  ! Initial size, will grow if needed
    end subroutine init_trajectory_reader

    function count_trajectory_frames(filename) result(num_frames)
        character(len=*), intent(in) :: filename
        integer :: num_frames
        integer :: unit, io_stat, grep_stat, cmd_stat
        logical :: file_exists
        character(256) :: command, result_string
        character(len=32) :: temp_file
        character(len=32) :: grep_positions_file
        
        num_frames = 0
        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            write(*,*) "Warning: Could not find file: ", trim(filename)
            return
        end if
        
        if (VERBOSE) write(*,'(A)') ' Counting frames in trajectory file using grep...'
        
        temp_file = "frame_count_temp.txt"
        grep_positions_file = "frame_positions_temp.txt"
        
        ! Count frames and get byte offsets for each frame
        command = 'grep -b "ITEM: NUMBER OF ATOMS" ' // trim(filename) // ' > ' // trim(grep_positions_file)
        call execute_command_line(trim(command), exitstat=grep_stat, cmdstat=cmd_stat)
        
        ! Get just the count
        command = 'grep -c "ITEM: NUMBER OF ATOMS" ' // trim(filename) // ' > ' // trim(temp_file)
        call execute_command_line(trim(command), exitstat=grep_stat, cmdstat=cmd_stat)
        
        ! If command execution fails or grep isn't available, fall back immediately
        if (grep_stat /= 0 .or. cmd_stat /= 0) then
            if (VERBOSE) write(*,*) "Using direct frame counting method"
            num_frames = count_frames_direct(filename)
            
            ! Clean up temporary files if they were created
            inquire(file=temp_file, exist=file_exists)
            if (file_exists) then
                call execute_command_line('rm ' // trim(temp_file), wait=.true.)
            end if
            
            inquire(file=grep_positions_file, exist=file_exists)
            if (file_exists) then
                call execute_command_line('rm ' // trim(grep_positions_file), wait=.true.)
            end if
            
            return
        end if
        
        ! Read the frame count
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
        
        if (VERBOSE) write(*,'(A,I0,A)') ' Found ', num_frames, ' frames in trajectory file'
        if (have_frame_offsets .and. VERBOSE) then
            write(*,'(A)') ' Frame offset table created for fast seeking'
        end if
    end function count_trajectory_frames
    
    ! Parse the output of grep -b to get file offsets for each frame
    subroutine parse_frame_offsets(filename, num_frames)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: num_frames
        integer :: unit, io_stat, i, alloc_stat
        character(256) :: line
        integer(kind=8) :: offset
        
        ! Allocate array for frame offsets
        if (allocated(frame_offsets)) deallocate(frame_offsets)
        allocate(frame_offsets(num_frames), stat=alloc_stat)
        if (alloc_stat /= 0) then
            if (VERBOSE) write(*,*) "Warning: Could not allocate frame offsets array"
            have_frame_offsets = .false.
            return
        end if
        
        ! Open the grep output file
        open(newunit=unit, file=filename, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            have_frame_offsets = .false.
            return
        end if
        
        ! Read each line and extract the byte offset
        do i = 1, num_frames
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            
            ! Extract offset (format is "OFFSET:ITEM: NUMBER OF ATOMS")
            read(line(:index(line, ':')-1), *, iostat=io_stat) offset
            if (io_stat == 0) then
                frame_offsets(i) = offset
            else
                if (VERBOSE) write(*,*) "Warning: Error parsing frame offset at line", i
            end if
        end do
        
        close(unit)
        
        ! Check if we got all offsets
        if (i > num_frames) then
            have_frame_offsets = .true.
        else
            if (VERBOSE) write(*,*) "Warning: Only found", i-1, "frame offsets of", num_frames
            have_frame_offsets = (i-1 >= num_frames)
        end if
    end subroutine parse_frame_offsets

    function count_frames_direct(fname) result(count)
        character(len=*), intent(in) :: fname
        integer :: count
        integer :: cnt_unit, cnt_stat
        character(256) :: cnt_line
        
        count = 0
        open(newunit=cnt_unit, file=fname, status='old', action='read', iostat=cnt_stat)
        if (cnt_stat /= 0) return
        
        do
            read(cnt_unit, '(A)', iostat=cnt_stat) cnt_line
            if (cnt_stat /= 0) exit
            if (index(cnt_line, 'ITEM: NUMBER OF ATOMS') > 0) count = count + 1
        end do
        
        close(cnt_unit)
    end function count_frames_direct
    
    ! New function to seek to a specific frame
    subroutine seek_to_frame(frame_number)
        integer, intent(in) :: frame_number
        integer :: i, io_stat
        character(256) :: dummy_line
        
        if (.not. input_file_open) then
            call handle_error("File not open for seeking", ERR_FILE_IO)
            return
        end if
        
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
        
        ! Fall back to sequential reading
        rewind(input_unit)
        if (frame_number <= 0) return
        
        ! Skip frames until we reach the desired one
        do i = 1, frame_number - 1
            ! Find the start of a frame
            do
                read(input_unit, '(A)', iostat=io_stat) dummy_line
                if (io_stat /= 0) then
                    call handle_error("Reached end of file before frame " // trim(integer_to_string(frame_number)), &
                                     ERR_FILE_IO)
                    return
                end if
                
                if (index(dummy_line, 'ITEM: NUMBER OF ATOMS') > 0) exit
            end do
            
            ! Skip this frame's content
            do
                read(input_unit, '(A)', iostat=io_stat) dummy_line
                if (io_stat /= 0) then
                    call handle_error("Reached end of file while skipping frame", ERR_FILE_IO)
                    return
                end if
                
                if (index(dummy_line, 'ITEM: NUMBER OF ATOMS') > 0) then
                    ! We found the start of the next frame, step back
                    backspace(input_unit)
                    exit
                end if
            end do
        end do
    end subroutine seek_to_frame
    
    ! Helper function to convert integer to string
    function integer_to_string(int_value) result(str)
        integer, intent(in) :: int_value
        character(len=20) :: str
        
        write(str, '(I0)') int_value
    end function integer_to_string

    subroutine read_trajectory_frame(coords, atom_types, elements, box_length, frame_number, eof, atom_info)
        real(dp), intent(out) :: coords(:,:)
        integer, intent(out) :: atom_types(:)
        character(len=*), intent(out) :: elements(:)
        real(dp), intent(out) :: box_length(3)
        integer, intent(inout) :: frame_number
        logical, intent(out) :: eof
        type(atom_type_info), intent(in) :: atom_info(:)
        
        integer :: n_atoms_frame, atom_id, io_stat, i, timestep
        real(dp) :: xlo, xhi, ylo, yhi, zlo, zhi
        character(256) :: line
        
        eof = .false.
        
        ! Look for frame start
        do
            read(input_unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) then
                eof = .true.
                return
            end if
            
            if (index(line, 'ITEM: NUMBER OF ATOMS') > 0) then
                read(input_unit, *) n_atoms_frame
                if (n_atoms_frame /= size(coords,1)) then
                    write(*,*) "Warning: Atom count mismatch -", &
                              "Expected:", size(coords,1), &
                              "Found:", n_atoms_frame
                end if
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
        
        ! Ensure values_buffer is large enough
        if (.not. allocated(values_buffer) .or. size(values_buffer) < col_map%max_cols) then
            if (allocated(values_buffer)) deallocate(values_buffer)
            allocate(values_buffer(col_map%max_cols))
        end if
        
        ! Read atom data
        do i = 1, n_atoms_frame
            read(input_unit, *, iostat=io_stat) values_buffer(1:col_map%max_cols)
            if (io_stat /= 0) then
                write(*,*) "Error reading atom", i, "in frame", frame_number
                call handle_error("Error reading coordinates", ERR_FILE_IO)
            end if
            
            atom_id = nint(values_buffer(col_map%id))
            if (atom_id < 1 .or. atom_id > size(coords,1)) then
                write(*,*) "Invalid atom ID:", atom_id, "at line", i
                call handle_error("Invalid atom ID in trajectory", ERR_INCONSISTENT_DATA)
            end if
            
            coords(atom_id,1) = values_buffer(col_map%x)
            coords(atom_id,2) = values_buffer(col_map%y)
            coords(atom_id,3) = values_buffer(col_map%z)
            atom_types(atom_id) = nint(values_buffer(col_map%type))
            
            ! Set element name from atom_info instead of generic "Type-X"
            if (atom_types(atom_id) >= 1 .and. atom_types(atom_id) <= size(atom_info)) then
                elements(atom_id) = atom_info(atom_types(atom_id))%name
            else
                write(elements(atom_id),'(A,I0)') 'Type', atom_types(atom_id)
            end if
        end do
        
        frame_number = frame_number + 1
    end subroutine read_trajectory_frame

    subroutine setup_column_mapping(header_line)
        character(len=*), intent(in) :: header_line
        integer :: i, pos, prev_pos, n_cols
        character(len=32) :: column_name
        
        ! Reset column mapping
        col_map = column_map_type()
        
        ! Skip "ITEM: ATOMS" prefix
        prev_pos = index(header_line, "ITEM: ATOMS") + 11
        
        ! Count columns
        n_cols = 0
        pos = prev_pos
        do while (pos <= len_trim(header_line))
            if (header_line(pos:pos) /= ' ') then
                n_cols = n_cols + 1
                pos = pos + 1
                do while (pos <= len_trim(header_line) .and. header_line(pos:pos) /= ' ')
                    pos = pos + 1
                end do
            else
                pos = pos + 1
            end if
        end do
        
        ! Map columns
        pos = prev_pos
        do i = 1, n_cols
            ! Skip spaces
            do while (pos <= len_trim(header_line) .and. header_line(pos:pos) == ' ')
                pos = pos + 1
            end do
            
            ! Extract column name
            prev_pos = pos
            do while (pos <= len_trim(header_line) .and. header_line(pos:pos) /= ' ')
                pos = pos + 1
            end do
            column_name = adjustl(header_line(prev_pos:pos-1))
            
            ! Map column
            select case (trim(column_name))
                case ('id')
                    col_map%id = i
                case ('type')
                    col_map%type = i
                case ('x')
                    col_map%x = i
                case ('y')
                    col_map%y = i
                case ('z')
                    col_map%z = i
            end select
        end do
        
        ! Set maximum columns needed
        col_map%max_cols = max(col_map%id, col_map%type, col_map%x, col_map%y, col_map%z)
        
        ! Verify required columns
        if (col_map%id < 0 .or. col_map%type < 0 .or. &
            col_map%x < 0 .or. col_map%y < 0 .or. col_map%z < 0) then
            call handle_error("Missing required columns in trajectory file", ERR_FILE_IO)
        end if
        
        if (VERBOSE) then
            write(*,'(A)') ' Column mapping in trajectory:'
            write(*,'(A,I0)') '   ID column: ', col_map%id
            write(*,'(A,I0)') '   Type column: ', col_map%type
            write(*,'(A,I0)') '   X column: ', col_map%x
            write(*,'(A,I0)') '   Y column: ', col_map%y
            write(*,'(A,I0)') '   Z column: ', col_map%z
        end if
    end subroutine setup_column_mapping

    subroutine cleanup_trajectory_reader()
        if (input_file_open) then
            close(input_unit)
            input_file_open = .false.
        end if
        
        if (allocated(frame_offsets)) deallocate(frame_offsets)
        if (allocated(values_buffer)) deallocate(values_buffer)
        have_frame_offsets = .false.
    end subroutine cleanup_trajectory_reader

end module trajectory_io_mod
