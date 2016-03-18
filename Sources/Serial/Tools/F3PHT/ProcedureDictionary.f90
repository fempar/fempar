!-----------------------------------------------------------------
! FPFPHT (Fortran Polymorphic Procedure Pointer Hash Table)
! Copyright (c) 2015 Santiago Badia, Alberto F. Martín, 
! Javier Principe and Víctor Sande.
! All rights reserved.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3.0 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library.
!-----------------------------------------------------------------

module ProcedureDictionary

USE iso_fortran_env, only: OUTPUT_UNIT
USE ProcedureHashTable

implicit none
private
save

    type :: ProcedureDictionary_t
    private
        type(ProcedureHashTable_t) :: Dictionary
    contains
    private
        procedure, public :: Set            => ProcedureDictionary_Set
        procedure, public :: Get            => ProcedureDictionary_Get
        procedure, public :: Del            => ProcedureDictionary_RemoveEntry
        procedure, public :: Init           => ProcedureDictionary_Init
        procedure, public :: isInitialized  => ProcedureDictionary_isInitialized
        procedure, public :: isPresent      => ProcedureDictionary_isPresent
        procedure, public :: Free           => ProcedureDictionary_Free
        procedure, public :: Print          => ProcedureDictionary_Print
        procedure, public :: Length         => ProcedureDictionary_Length
        final             ::                   ProcedureDictionary_Finalize
    end type ProcedureDictionary_t

public :: ProcedureDictionary_t

contains


    subroutine ProcedureDictionary_Init(this,Size)
    !-----------------------------------------------------------------
    !< Initialize the dictionary
    !-----------------------------------------------------------------
        class(ProcedureDictionary_t),      intent(INOUT) :: this      !< Procedure Dictionary
        integer, optional,                 intent(IN)    :: Size      !< Dictionary Size
    !-----------------------------------------------------------------
        call this%Free()
        if (present(Size)) then
            call this%Dictionary%Init(Size = Size)
        else
            call this%Dictionary%Init()
        endif
    end subroutine ProcedureDictionary_Init


    function ProcedureDictionary_isInitialized(this) result(Initialized)
    !-----------------------------------------------------------------
    !< Initialize the dictionary
    !-----------------------------------------------------------------
        class(ProcedureDictionary_t), intent(IN) :: this              !< Procedure Dictionary
        logical                                  :: Initialized       !< Boolean flag to check if the dictionary is initialized               
    !-----------------------------------------------------------------
        Initialized = this%Dictionary%isInitialized()
    end function ProcedureDictionary_isInitialized


    subroutine ProcedureDictionary_Free(this)
    !-----------------------------------------------------------------
    !< Free the dictionary
    !-----------------------------------------------------------------
        class(ProcedureDictionary_t),      intent(INOUT) :: this      !< Procedure Dictionary
    !-----------------------------------------------------------------
        call this%Dictionary%Free()
    end subroutine ProcedureDictionary_Free


    subroutine ProcedureDictionary_Finalize(this)
    !-----------------------------------------------------------------
    !< Destructor procedure
    !-----------------------------------------------------------------
        type(ProcedureDictionary_t),       intent(INOUT) :: this      !< Procedure Hash Table List
    !-----------------------------------------------------------------
        call this%Free()
    end subroutine ProcedureDictionary_Finalize


    subroutine ProcedureDictionary_Set(this,Key,Proc)
    !-----------------------------------------------------------------
    !< Set a Key/Value pair into the Dictionary
    !-----------------------------------------------------------------
        class(ProcedureDictionary_t),      intent(INOUT) :: this      !< Procedure Hash Table List
        character(len=*),                  intent(IN)    :: Key       !< String Key
        procedure(),           pointer                   :: Proc      !< Polymorphic procedure pointer
    !-----------------------------------------------------------------
        call this%Dictionary%Set(Key=Key,Proc=Proc)
    end subroutine ProcedureDictionary_Set


    subroutine ProcedureDictionary_Get(this,Key,Proc)
    !-----------------------------------------------------------------
    !< Return a scalar Value given the Key
    !-----------------------------------------------------------------
        class(ProcedureDictionary_t),      intent(IN)    :: this      !< Procedure Hash Table List
        character(len=*),                  intent(IN)    :: Key       !< String Key
        procedure(),      pointer                        :: Proc      !< Returned procedure pointer
    !-----------------------------------------------------------------
        call this%Dictionary%Get(Key=Key,Proc=Proc)
    end subroutine ProcedureDictionary_Get


    function ProcedureDictionary_isPresent(this,Key) result(isPresent)
    !-----------------------------------------------------------------
    !< Check if a Key is present at the DataBase
    !-----------------------------------------------------------------
        class(ProcedureDictionary_t),      intent(IN) :: this         !< Procedure Hash Table List
        character(len=*),                  intent(IN) :: Key          !< String Key
        logical                                       :: isPresent    !< Boolean flag to check if a Key is present
    !-----------------------------------------------------------------
        isPresent = this%Dictionary%IsPresent(Key=Key)
    end function ProcedureDictionary_isPresent


    subroutine ProcedureDictionary_RemoveEntry(this, Key)
    !-----------------------------------------------------------------
    !< Remove an Entry given a Key
    !-----------------------------------------------------------------
        class(ProcedureDictionary_t),      intent(INOUT) :: this      !< Procedure Hash Table List
        character(len=*),                  intent(IN)    :: Key       !< String Key
    !-----------------------------------------------------------------
        call this%Dictionary%Del(Key=Key)
    end subroutine ProcedureDictionary_RemoveEntry


    function ProcedureDictionary_Length(this) result(Length)
    !-----------------------------------------------------------------
    !< Return the number of ProcedureDictionaryEntries contained in the DataBase
    !-----------------------------------------------------------------
        class(ProcedureDictionary_t),      intent(IN) :: this         !< Procedure DictionaryEntry Containter type
        integer                                       :: Length       !< Number of parameters in database
        integer                                       :: DBIterator   !< Database Iterator index 
    !-----------------------------------------------------------------
        Length = this%Dictionary%Length()
    end function ProcedureDictionary_Length


    subroutine ProcedureDictionary_Print(this, unit, prefix, iostat, iomsg)
    !-----------------------------------------------------------------
    !< Print the content of the DataBase
    !-----------------------------------------------------------------
        class(ProcedureDictionary_t),         intent(IN)  :: this    !< Linked List
        integer,      optional,               intent(IN)  :: unit    !< Logic unit.
        character(*), optional,               intent(IN)  :: prefix  !< Prefixing string.
        integer,      optional,               intent(OUT) :: iostat  !< IO error.
        character(*), optional,               intent(OUT) :: iomsg   !< IO error message.
        character(len=:), allocatable                     :: prefd   !< Prefixing string.
        integer                                           :: unitd   !< Logic unit.
        integer                                           :: iostatd !< IO error.
        character(500)                                    :: iomsgd  !< Temporary variable for IO error message.
    !-----------------------------------------------------------------
        prefd = '' ; if (present(prefix)) prefd = prefix
        unitd = OUTPUT_UNIT; if(present(unit)) unitd = unit
        write(*,fmt='(A)') prefd//' Procedure Dictionary Content:'
        write(*,fmt='(A)') prefd//' -----------------------------'
        call this%Dictionary%Print(unit=unitd, prefix=prefd, iostat=iostatd, iomsg=iomsgd)
        if (present(iostat)) iostat = iostatd
        if (present(iomsg))  iomsg  = iomsgd
    end subroutine ProcedureDictionary_Print


end module ProcedureDictionary
