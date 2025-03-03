module error_mod
    use iso_fortran_env, only: error_unit
    implicit none
    private

    ! Public procedures
    public :: handle_error
    public :: check_allocation
    public :: validate_parameters

    ! Error codes
    integer, parameter, public :: ERR_ALLOCATION = 1
    integer, parameter, public :: ERR_FILE_IO = 2
    integer, parameter, public :: ERR_INVALID_PARAM = 3
    integer, parameter, public :: ERR_INVALID_TYPE = 4
    integer, parameter, public :: ERR_INCONSISTENT_DATA = 5

contains

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

    subroutine validate_parameters(param_name, param_value, min_value, max_value)
        character(len=*), intent(in) :: param_name
        real, intent(in) :: param_value
        real, intent(in) :: min_value
        real, intent(in) :: max_value
        
        if (param_value < min_value .or. param_value > max_value) then
            write(error_unit,'(5A,2(A,G0.4))') "Parameter '", trim(param_name), &
                "' value ", "out of range [", trim(param_name), &
                "] = ", param_value, " not in [", min_value, ",", max_value, "]"
            call handle_error("Parameter validation failed", ERR_INVALID_PARAM)
        end if
    end subroutine validate_parameters

end module error_mod