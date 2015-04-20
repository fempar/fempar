! Copyright (C) 2014 Santiago Badia, Alberto F. Martín and Javier Principe
!
! This file is part of FEMPAR (Finite Element Multiphysics PARallel library)
!
! FEMPAR is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! FEMPAR is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with FEMPAR. If not, see <http://www.gnu.org/licenses/>.
!
! Additional permission under GNU GPL version 3 section 7
!
! If you modify this Program, or any covered work, by linking or combining it 
! with the Intel Math Kernel Library and/or the Watson Sparse Matrix Package 
! and/or the HSL Mathematical Software Library (or a modified version of them), 
! containing parts covered by the terms of their respective licenses, the
! licensors of this Program grant you additional permission to convey the 
! resulting work. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# include "debug.i90"
module stdio
  ! This module contains functions to manage io. It was 
  ! created by C. Labra (clabra@cimne.upc.edu) and slightly
  ! modified to fit our environment. The modifications are
  !
  ! 1) function io_open is modified making action optional and
  !    default to READWRITE
  !
  use, intrinsic :: iso_fortran_env, only: ERROR_UNIT, INPUT_UNIT, &
       &                                   OUTPUT_UNIT, IOSTAT_END, IOSTAT_EOR
  implicit none

  ! private module
  private

  integer, parameter :: IOSTAT_OK = 0  ! extra parameter for iso_fortran_env

  integer, parameter :: BUFFER_STRING_SIZE  = 30
  integer, parameter :: BUFFER_SIZE         = 1024
  integer, parameter :: CHAR_SIZE           = 56
  integer, parameter :: NAME_SIZE           = 128

  character(*), parameter :: REAL4_FORMAT = '(F8.4)'
  character(*), parameter :: REAL8_FORMAT = '(E16.6)'
  character(*), parameter :: INT_FORMAT = '(I16)'

  integer, parameter :: FILE_SIZE     = 256
  integer, parameter :: FILE_NULL     = -1

  integer,      parameter :: FILENAME_SIZE = 512
  character(*), parameter :: FILENAME_NULL = ''

  integer, parameter :: MIN_UNIT = 11
  integer, parameter :: MAX_UNIT = 9999

  character(*), parameter ::                                 &
       COMMENT_SYMBOL   = '#$!',                         &
       STRING_SEPARATOR = ' =:'//achar(9),               &
       CHAR_LOWERCASE   = 'abcdefghijklmnopqrstuvwxyz',  &
       CHAR_UPPERCASE   = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ',  &
       CHAR_NUMERIC     = '0123456789-+.'!,               &
!       CHAR_NEXTLINE    = '\'
       
       character(*), parameter :: VALID_CHARS =  &
       CHAR_LOWERCASE//CHAR_UPPERCASE//CHAR_NUMERIC//COMMENT_SYMBOL !//CHAR_NEXTLINE


  character(CHAR_SIZE), parameter :: CHAR_NULL = ''
  character(NAME_SIZE), parameter :: NAME_NULL = ''

  ! standart units
  integer, save :: stdin  = INPUT_UNIT  ! iso_fortran_env
  integer, save :: stdout = OUTPUT_UNIT ! iso_fortran_env
  integer, save :: stderr = ERROR_UNIT  ! iso_fortran_env
  character(FILE_SIZE), save :: inputfile  = ''
  character(FILE_SIZE), save :: outputfile = ''
  character(FILE_SIZE), save :: errorfile  = ''
  logical, save :: stdin_file  = .false.
  logical, save :: stdout_file = .false.
  logical, save :: stderr_file = .false.

  character(*), parameter :: newline = achar(10)

  type buffer_t
     integer :: length
     character(CHAR_SIZE), dimension(BUFFER_STRING_SIZE) :: string
  end type buffer_t

#if defined(_DEBUG)
  logical, parameter :: stdio_debug = .true.
#else
  logical, parameter :: stdio_debug = .false.
#endif


  interface assignment (=)
     module procedure char_to_real4
     module procedure char_to_real8
     module procedure char_to_integer
     module procedure char_to_logical
  end interface assignment (=)

  interface ch
     module procedure real4_to_char
     module procedure real4_to_char_fmt
     module procedure real8_to_char
     module procedure real8_to_char_fmt
     module procedure integer1_to_char
     module procedure integer_to_char
     module procedure integer_to_char_fmt
     module procedure logical_to_char
  end interface ch

  interface operator(//)
     module procedure cat_char_int
     module procedure cat_char_int1
     module procedure cat_char_real4
     module procedure cat_char_real8
     module procedure cat_char_logical
     module procedure cat_int_char
     module procedure cat_real4_char
     module procedure cat_real8_char
     module procedure cat_logical_char
  end interface operator(//)

  interface io_get
     module procedure io_get_int
     module procedure io_get_int_array
     module procedure io_get_r4
     module procedure io_get_r4_array
     module procedure io_get_r8
     module procedure io_get_r8_array
     module procedure io_get_line
     module procedure io_get_buffer
  end interface io_get

  interface io_get_screen
     module procedure io_get_screen_char
     module procedure io_get_screen_int
     module procedure io_get_screen_real4
     module procedure io_get_screen_real8
  end interface io_get_screen

  interface lowercase
     module procedure lowcase
  end interface lowercase

  interface print
     module procedure print_screen
     module procedure print_file
  end interface print

  ! public types
  public :: buffer_t

  ! public paramters
  public :: BUFFER_STRING_SIZE, BUFFER_SIZE, CHAR_SIZE, FILE_SIZE, &
       CHAR_NULL, NAME_SIZE, NAME_NULL, FILE_NULL, FILENAME_SIZE, FILENAME_NULL, &
       NEWLINE

  ! public routines
  public :: assignment(=), io_get, io_open, is_number, io_file_opened,  &
       screen, io_get_screen, io_file_exist, ch, io_end, io_eor,   &
       io_get_extension, io_remove_extension, io_mkdir, io_remove, &
       io_pwd, io_rewind, io_close, io_close_all, io_get_unit,     &
       lowercase, operator(//), numbered_filename_compose, numbered_filename_compose_deferred_length, &
       count_digits

  ! public units
  public :: stdin, stdout, stderr, io_open_stdin, io_open_stdout,  &
       io_open_stderr, stdin_file, stdout_file, stderr_file

contains

  !********************************************************************!

  subroutine char_to_real4( value, string )

    real(4),      intent(out) :: value
    character(*), intent(in)  :: string
    integer :: io_status

    read(string,fmt=*,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'STRING_TO_REAL4','error reading string :'//string)

  end subroutine char_to_real4

  !********************************************************************!

  subroutine char_to_real8( value, string )

    real(8),      intent(out) :: value
    character(*), intent(in)  :: string
    integer :: io_status

    read(string,fmt=*,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'STRING_TO_REAL8','error reading string :'//string)

  end subroutine char_to_real8

  !********************************************************************!

  subroutine char_to_integer( value, string )

    integer,      intent(out) :: value
    character(*), intent(in)  :: string
    integer :: io_status

    read(string,fmt=*,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'STRING_TO_INTEGER','error reading string :'//string)

  end subroutine char_to_integer

  !********************************************************************!

  subroutine char_to_logical( value, string )

    logical,      intent(out) :: value
    character(*), intent(in)  :: string
    integer :: io_status

    read(string,fmt=*,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'STRING_TO_INTEGER','error reading string :'//string)

  end subroutine char_to_logical

  !********************************************************************!

  function real4_to_char( value ) result( string )

    real(kind=4), intent(in) :: value
    character(CHAR_SIZE) :: string
    integer :: io_status

    write(string,fmt=REAL4_FORMAT,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'REAL4_TO_STRING',' fatal error.')
    string = adjustl(string)

  end function real4_to_char

  !********************************************************************!

  function real4_to_char_fmt( value, io_format ) result( string )

    real(kind=4), intent(in) :: value
    character(*), intent(in) :: io_format
    character(CHAR_SIZE) :: string
    integer :: io_status

    write(string,fmt=io_format,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'REAL4_TO_STRING',' fatal error.')
    string = adjustl(string)

  end function real4_to_char_fmt

  !********************************************************************!

  function real8_to_char( value ) result( string )

    real(kind=8), intent(in) :: value
    character(CHAR_SIZE) :: string
    integer :: io_status

    write(string,fmt=REAL8_FORMAT,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'REAL8_TO_STRING',' fatal error.')
    string = adjustl(string)

  end function real8_to_char

  !********************************************************************!

  function real8_to_char_fmt( value, io_format ) result( string )

    real(kind=8), intent(in) :: value
    character(*), intent(in) :: io_format
    character(CHAR_SIZE) :: string
    integer :: io_status

    write(string,fmt=io_format,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'REAL8_TO_STRING',' fatal error.')
    string = adjustl(string)

  end function real8_to_char_fmt

  !********************************************************************!

  function integer_to_char( value ) result( string )

    integer,  intent(in) :: value
    character(CHAR_SIZE) :: string
    integer :: io_status

    write(string,fmt=INT_FORMAT,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'INTEGER_TO_STRING',' fatal error.')
    string = adjustl(string)

  end function integer_to_char

  !********************************************************************!

  function integer1_to_char( value ) result( string )

    integer(1), intent(in) :: value
    character(CHAR_SIZE)   :: string
    integer :: io_status

    write(string,fmt=INT_FORMAT,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'INTEGER_TO_STRING',' fatal error.')
    string = adjustl(string)

  end function integer1_to_char

  !********************************************************************!

  function integer_to_char_fmt( value, io_format ) result( string )

    integer,      intent(in) :: value
    character(*), intent(in) :: io_format
    character(CHAR_SIZE) :: string
    integer :: io_status

    write(string,fmt=io_format,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'INTEGER_TO_STRING',' fatal error.')
    string = adjustl(string)

  end function integer_to_char_fmt

  !********************************************************************!

  function logical_to_char( value ) result( string )

    logical, intent(in) :: value
    character(CHAR_SIZE) :: string
    integer :: io_status

    if( value ) then
       write(string,fmt='(A)',iostat=io_status) 'TRUE'
    else
       write(string,fmt='(A)',iostat=io_status) 'FALSE'
    end if
    if(stdio_debug) &
         call iostat_error(io_status,'LOGICAL_TO_STRING',' fatal error.')

  end function logical_to_char

  !********************************************************************!

  function cat_char_int1(string,number) result(new_string)
    character(*),  intent(in) :: string
    integer(1),    intent(in) :: number
    character(buffer_size) :: new_string
    new_string = string//trim(integer1_to_char(number))
  end function cat_char_int1

  !********************************************************************!

  function cat_char_int(string,number) result(new_string)
    character(*),  intent(in) :: string
    integer,       intent(in) :: number
    character(buffer_size) :: new_string
    new_string = string//trim(integer_to_char(number))
  end function cat_char_int

  !********************************************************************!

  function cat_int_char(number,string) result(new_string)
    character(*),  intent(in) :: string
    integer,       intent(in) :: number
    character(buffer_size) :: new_string
    new_string = trim(integer_to_char(number))//string
  end function cat_int_char

  !********************************************************************!

  function cat_char_real4(string,number) result(new_string)
    character(*),  intent(in) :: string
    real(4),       intent(in) :: number
    character(buffer_size) :: new_string
    new_string = string//trim(real4_to_char(number))
  end function cat_char_real4

  !********************************************************************!

  function cat_real4_char(number,string) result(new_string)
    character(*),  intent(in) :: string
    real(4),       intent(in) :: number
    character(buffer_size) :: new_string
    new_string = trim(real4_to_char(number))//string
  end function cat_real4_char

  !********************************************************************!

  function cat_char_real8(string,number) result(new_string)
    character(*),  intent(in) :: string
    real(8),       intent(in) :: number
    character(buffer_size) :: new_string
    new_string = string//trim(real8_to_char(number))
  end function cat_char_real8

  !********************************************************************!

  function cat_real8_char(number,string) result(new_string)
    character(*),  intent(in) :: string
    real(8),       intent(in) :: number
    character(buffer_size) :: new_string
    new_string = trim(real8_to_char(number))//string
  end function cat_real8_char

  !********************************************************************!

  function cat_char_logical(string,number) result(new_string)
    character(*),  intent(in) :: string
    logical,       intent(in) :: number
    character(buffer_size) :: new_string
    new_string = string//trim(logical_to_char(number))
  end function cat_char_logical

  !********************************************************************!

  function cat_logical_char(number,string) result(new_string)
    character(*),  intent(in) :: string
    logical,       intent(in) :: number
    character(buffer_size) :: new_string
    new_string = trim(logical_to_char(number))//string
  end function cat_logical_char

  !********************************************************************!

  subroutine iostat_error(iostat,function,message)

    integer,      intent(in) :: iostat
    character(*), optional, intent(in) :: function
    character(*), optional, intent(in) :: message

    if( iostat /= IOSTAT_OK ) then
       write(stderr,'(A)') 'IOSTAT_ERROR: Fatal!'
       if(present(function)) write(stderr,'(A)') ' Function : '//trim(function)
       if(present(message))  write(stderr,'(A)') '  comment : '//trim(message)
       stop
    end if

  end subroutine iostat_error

  !********************************************************************!

  subroutine iostat_warning(iostat,message)

    integer,      intent(in) :: iostat
    character(*), intent(in) :: message

    if( iostat /= IOSTAT_OK ) then
       write(stderr,'(A)') message
       write(stderr,'(A)') ' Warning!.'
    end if

  end subroutine iostat_warning

  !********************************************************************!

  integer function io_open( file, action, status, form, position, recl )

    character(len=*), intent(in)           :: file
    character(len=*), intent(in), optional :: action, status, form, position
    integer,          intent(in), optional :: recl

    character(CHAR_SIZE)   :: action_, status_, form_, position_
    integer                :: iostat, io

    action_ = 'readwrite'
    if(present(action  )) action_   = action
    status_ = 'unknown'
    if(present(status  )) status_   = status
    form_   = 'formatted'
    if(present(form    )) form_     = form
    position_ = 'asis'
    if(present(position)) position_ = position

    io = io_get_unit()

    if(present(recl)) then
       open(unit=io, file=trim(file), status=trim(status_), form=trim(form_), &
            recl=recl, action=trim(action_), position=trim(position_), iostat=iostat)
    else
       open(unit=io, file=trim(file), status=trim(status_), form=trim(form_), &
            action=trim(action_), position=trim(position_), iostat=iostat)
    endif
    call iostat_error(iostat,'IO_OPEN','could not open file: '//file)
    io_open = io

  end function io_open

  !********************************************************************!

  subroutine io_rewind( io )
    integer, intent(in) :: io
    integer :: io_status

    rewind( io, iostat=io_status )
    if(stdio_debug) &
         call iostat_error(io_status,'IO_REWIND',' unit = '//io)

  end subroutine io_rewind

  !********************************************************************!

  subroutine io_close( io, status )
    integer,                intent(in) :: io
    character(*), optional, intent(in) :: status
    integer :: io_status

    if( present(status) ) then
       close(unit=io,iostat=io_status, status=status)
    else
       close(unit=io,iostat=io_status)
    end if
    if(stdio_debug) &
         call iostat_error(io_status,'IO_CLOSE','unit = '//io)

  end subroutine io_close

  !********************************************************************!

  subroutine io_close_all()

    integer :: io, io_status
    logical :: is_opened

    do io = MIN_UNIT, MAX_UNIT
       inquire (unit=io, opened=is_opened, iostat=io_status )
       if(stdio_debug) &
            call iostat_error(io_status,'IO_CLOSE_ALL',' unit = '//io)
       if(is_opened) call io_close(io)
    end do

  end subroutine io_close_all

  !********************************************************************!

  subroutine io_get_line( io, text )

    integer,          intent(in)  :: io
    character(*), intent(out) :: text

    integer :: stat, i

    do
       read( io, fmt='(a)', iostat=stat ) text
       if(stat/=IOSTAT_OK) return
       i = max(scan(text,VALID_CHARS),1)
       if( verify(text(i:i),COMMENT_SYMBOL)==1 .and. len_trim(text) > 0 ) exit
    end do

    text = text(i:)
    call lowcase( text )

  end subroutine io_get_line

  !********************************************************************!

  subroutine io_get_buffer( io, buffer )

    integer, intent(in) :: io
    type(buffer_t), intent(out) :: buffer

    integer :: i, nr_words, n, len, next
    character(BUFFER_SIZE) :: text


    buffer%length = 0
    nr_words = 0

    i = 1
    do
       call io_get_line( io, text(i:) )
       len = len_trim(text)
!!$       if( text(len:len) == CHAR_NEXTLINE ) then
!!$          i = len
!!$       else
!!$          exit
!!$       end if
    end do

    !i = 1
    i = verify( text(1:len), STRING_SEPARATOR )
    word_loop: do
       n = scan( text(i:len), STRING_SEPARATOR )
       if ( n == 0 ) then
          nr_words = nr_words + 1
          buffer%string(nr_words) = text(i:len)
          exit word_loop
       end if

       nr_words = nr_words + 1
       buffer%string(nr_words) = text(i:i+n-2)

       next = verify( text(i+n-1:len), STRING_SEPARATOR )
       if ( next == 0 ) exit

       i = next + i + n - 2
    end do word_loop

    buffer%length = nr_words

  end subroutine io_get_buffer

  !********************************************************************!

  subroutine io_get_int( io, value )

    integer, intent(in)  :: io
    integer, intent(out) :: value

    integer :: io_status

    read(unit=io,fmt=*,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'IO_GET_INTEGER','fatal error.')

  end subroutine io_get_int

  !********************************************************************!

  subroutine io_get_r4( io, value )
    integer,      intent(in)  :: io
    real(kind=4), intent(out) :: value
    integer :: io_status

    read(unit=io,fmt=*,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'IO_GET_REAL4','fatal error.')

  end subroutine io_get_r4

  !********************************************************************!

  subroutine io_get_r8( io, value )
    integer,      intent(in)  :: io
    real(kind=8), intent(out) :: value
    integer :: io_status

    read(unit=io,fmt=*,iostat=io_status) value
    if(stdio_debug) &
         call iostat_error(io_status,'IO_GET_REAL8','fatal error.')

  end subroutine io_get_r8

  !********************************************************************!

  subroutine io_get_int_array( io, array )

    integer, intent(in)  :: io
    integer, intent(out) :: array(:)
    integer :: io_status

    read(unit=io,fmt=*,iostat=io_status) array
    if(stdio_debug) &
         call iostat_error(io_status,'IO_GET_INTEGER_ARRAY','fatal error.')

  end subroutine io_get_int_array

  !********************************************************************!

  subroutine io_get_r4_array( io, array )

    integer, intent(in)  :: io
    real(4), intent(out) :: array(:)
    integer :: io_status

    read(unit=io,fmt=*,iostat=io_status) array
    if(stdio_debug) &
         call iostat_error(io_status,'IO_GET_REAL4_ARRAY','fatal error.')

  end subroutine io_get_r4_array

  !********************************************************************!

  subroutine io_get_r8_array( io, array )

    integer, intent(in)  :: io
    real(8), intent(out) :: array(:)
    integer :: io_status

    read(unit=io,fmt=*,iostat=io_status) array
    if(stdio_debug) &
         call iostat_error(io_status,'IO_GET_REAL8_ARRAY','fatal error.')

  end subroutine io_get_r8_array

  !********************************************************************!

  subroutine io_open_stdin( file )

    character(*), intent(in) :: file

    inputfile = trim(file)
    stdin = io_open(file,action='read')
    stdin_file = .true.

  end subroutine io_open_stdin

  !********************************************************************!

  subroutine io_open_stdout( file )

    character(*), intent(in) :: file

    outputfile = trim(file)
    stdout = io_open(file,action='write')
    stdout_file = .true.

  end subroutine io_open_stdout

  !********************************************************************!

  subroutine io_open_stderr( file )

    character(len=*), intent(in) :: file

    outputfile = trim(file)
    stderr     = io_open(file,action='write')
    stderr_file = .true.

  end subroutine io_open_stderr

  !********************************************************************!

  subroutine screen( message, advance )

    character(len=*),           intent(in) :: message
    logical,          optional, intent(in) :: advance

    if ( present(advance) ) then
       if( advance ) then
          write(stdout,'(A)') trim(message)
       else
          write(stdout,'(A)', advance='no') trim(message)
       end if
    else
       write(stdout,'(A)') trim(message)
    end if

  end subroutine screen

  !********************************************************************!

  subroutine print_screen( message, advance )

    character(len=*),           intent(in) :: message
    logical,          optional, intent(in) :: advance

    if ( present(advance) ) then
       if( advance ) then
          write(stdout,'(A)') trim(message)
       else
          write(stdout,'(A)', advance='no') trim(message)
       end if
    else
       write(stdout,'(A)') trim(message)
    end if

  end subroutine print_screen

  !********************************************************************!

  subroutine print_file( ounit, message, advance )

    integer,                    intent(in) :: ounit
    character(len=*),           intent(in) :: message
    logical,          optional, intent(in) :: advance

    if ( present(advance) ) then
       if( advance ) then
          write(ounit,'(A)') trim(message)
       else
          write(ounit,'(A)', advance='no') trim(message)
       end if
    else
       write(ounit,'(A)') trim(message)
    end if

  end subroutine print_file

  !********************************************************************!

  subroutine io_get_screen_char( question, text, len )

    character(*), intent(in)  :: question
    character(*), intent(out) :: text
    integer,      intent(out) :: len

    write(stdout,'(A)') question
    read (stdin, '(A)') text
    len = len_trim( text )

  end subroutine io_get_screen_char

  !********************************************************************!

  subroutine io_get_screen_int( question, value )

    character(*), intent(in)  :: question
    integer,      intent(out) :: value

    write(stdout,'(A)') question
    read (stdin,*) value

  end subroutine io_get_screen_int

  !********************************************************************!

  subroutine io_get_screen_real4( question, value )

    character(len=*), intent(in)  :: question
    real(kind=4),     intent(out) :: value

    write(stdout,'(A)') question
    read (stdin,*) value

  end subroutine io_get_screen_real4

  !********************************************************************!

  subroutine io_get_screen_real8( question, value )

    character(len=*), intent(in)  :: question
    real(kind=8),     intent(out) :: value

    write(stdout,'(A)') question
    read (stdin,*) value

  end subroutine io_get_screen_real8

  !********************************************************************!

  logical function io_file_exist( file )
    character(*), intent(in) :: file

    inquire(file=trim(file),exist=io_file_exist)

  end function io_file_exist

  !********************************************************************!

  logical function io_file_opened( file, unit )
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    if( present(unit) ) then
       inquire( unit=unit , opened=io_file_opened )
    else if ( present(file) ) then
       inquire( file=file , opened=io_file_opened )
    else
       write(stderr,'(A)') 'io_opened error: empty dummy argument.'
    end if

  end function io_file_opened

  !********************************************************************!

  subroutine io_get_extension( filename, ext )
    character(*), intent(in)  :: filename
    character(*), intent(out) :: ext
    integer :: i, j

    i = index(filename, '.', back = .true.)
    j = index(filename(i+1:), '/')
    if( i==0 .or. j==0 ) then
       ext = CHAR_NULL
    else
       ext = filename(i+1:)
    end if

  end subroutine io_get_extension

  !********************************************************************!

  subroutine io_remove_extension( filename )
    character(*), intent(inout) :: filename
    integer :: i

    i = index(filename, '.', back = .true.)
    if( i/=0 ) filename(i:) = CHAR_NULL

  end subroutine io_remove_extension

  !********************************************************************!

  elemental logical function io_end( io_status )
    integer, intent(in) :: io_status
    io_end = ( io_status == IOSTAT_END )
  end function io_end

  !********************************************************************!

  elemental logical function io_eor( io_status )
    integer, intent(in) :: io_status
    io_eor = ( io_status == IOSTAT_EOR )
  end function io_eor

  !********************************************************************!

  logical function is_number( string )

    character(*), intent(in) :: string

    real(kind=8) :: value
    integer :: io_status

    read(string,*,iostat=io_status) value
    is_number = ( io_status == IOSTAT_OK )

  end function is_number

  !********************************************************************!

  integer function io_get_unit()

    integer :: io_unit, io_status
    logical :: is_opened

    io_get_unit = FILE_NULL
    do io_unit = MIN_UNIT, MAX_UNIT
       inquire ( unit = io_unit, opened = is_opened, iostat = io_status )
       if ( io_status == IOSTAT_OK .and. .not.is_opened ) then
          io_get_unit = io_unit
          exit
       end if
    end do
    if(io_get_unit == FILE_NULL) then
       write(stderr, '(A)') 'io_open error: Too many files opened'
       stop
    end if

  end function io_get_unit

  !********************************************************************!

  subroutine io_mkdir( dir )

    character(*), intent(in) :: dir

    write(stdout,'(A)') ' io_mkdir warning: subroutine not implemented'
    write(stdout,'(A)') '        directory: '//dir

  end subroutine io_mkdir

  !********************************************************************!

  subroutine io_remove( filename )

    character(*), intent(in) :: filename

    logical :: exist, is_opened
    integer :: io_unit, io_status

    !write(stdout,'(A)') ' io_remove warning: subroutine not implemented'

    inquire(file=filename,number=io_unit,exist=exist,opened=is_opened,iostat=io_status)

    if(stdio_debug) &
         call iostat_error(io_status,'IO_REMOVE',' inquire error on file: '//filename)

    if (exist) then
       if (.not.is_opened) io_unit = io_open(filename,action='read')
       call io_close(io_unit,status='DELETE')
    end if

  end subroutine io_remove

  !********************************************************************!

  subroutine io_pwd( path )
    character(*), intent(out) :: path

    path = CHAR_NULL
    write(stdout,'(A)') ' io_pwd warning: subroutine not implemented'

  end subroutine io_pwd

  !********************************************************************!

  elemental subroutine lowcase( str )

    character(*), intent(inout) :: str
    integer :: i,j

    do i=1,len_trim(str)
       j = index( CHAR_UPPERCASE, str(i:i) )
       if( j /= 0 ) str(i:i) = CHAR_LOWERCASE(j:j)
    end do

  end subroutine lowcase

  !********************************************************************!
  subroutine numbered_filename_compose( i, n, file ) ! AFM: to be removed in the future
    implicit none 
    ! Parameters
    integer         , intent(in)    :: i, n 
    character(len=*), intent(inout) :: file

    integer         :: j, ndigs_i, ndigs_n
    character(256)  :: zeros
    character(256)  :: chari


    ndigs_n = count_digits( n )
    zeros = ' '
    ndigs_i = count_digits( i )

    do j= 1, ndigs_n - ndigs_i
       zeros (j:j) = '0'
    end do
    chari = ch(i)
   
    file = trim(file) // '.' // trim(zeros) // trim(chari)

  end subroutine numbered_filename_compose

  !********************************************************************!
  subroutine numbered_filename_compose_deferred_length( i, n, file )
    implicit none 
    ! Parameters
    integer                      , intent(in)    :: i, n 
    character(len=:), allocatable, intent(inout) :: file

    integer          :: j, ndigs_i, ndigs_n, istat
    character(len=:), allocatable :: zeros
    character(len=:), allocatable :: chari


    ndigs_n = count_digits( n )
    ndigs_i = count_digits( i )

    allocate(character(ndigs_n - ndigs_i) :: zeros, stat=istat)
    check(istat==0)

    do j= 1, ndigs_n - ndigs_i
       zeros (j:j) = '0'
    end do
    chari = ch(i)
   
    file = trim(file) // '.' // trim(zeros) // trim(chari)

    deallocate(zeros, stat=istat)
    check(istat==0)

    deallocate(chari, stat=istat)
    check(istat==0)

  end subroutine numbered_filename_compose_deferred_length

  !********************************************************************!

  function count_digits ( i )
    implicit none
    ! Parameters
    integer, intent(in) :: i 
    integer             :: count_digits
    ! Locals   
    integer             :: x 
    x = i 
    if (x < 0) x = -x;
    count_digits = 1;
    x = x/10;
    do while( x > 0)
       count_digits = count_digits + 1
       x = x/10;
    end do
  end function count_digits

end module stdio
