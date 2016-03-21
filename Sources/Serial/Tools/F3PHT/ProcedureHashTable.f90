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

!-----------------------------------------------------------------
! ProcedureHashTable is a datatype containing a Database
! array of ProcedureListRoot made to store diferent Entries
! depending on the hash of its Key.
!
! This work takes as a starting point the previou work of
! Stefano Zaghi (@szaghi, https://github.com/szaghi).
!
! You can find the original source at:
! https://github.com/szaghi/OFF/blob/95691ca15e6d68128ba016e40df74e42123f1c54/src/Data_Type_Hash_Table.f90
!-----------------------------------------------------------------

module ProcedureHashTable

USE ProcedureListNode
USE ProcedureListRoot

implicit None
private

    integer, parameter:: DefaultDataBaseSize = 100

    type :: ProcedureHashTable_t
    private
        type(ProcedureListRoot_t), allocatable :: DataBase(:)
        integer                            :: Size = 0
    contains
    private
        procedure, non_overridable         :: Hash          => ProcedureHashTable_Hash
        procedure, non_overridable, public :: Init          => ProcedureHashTable_Init
        procedure, non_overridable, public :: isInitialized => ProcedureHashTable_isInitialized
        procedure, non_overridable, public :: Set           => ProcedureHashTable_Set
        procedure, non_overridable, public :: Get           => ProcedureHashTable_Get
        procedure, non_overridable, public :: Del           => ProcedureHashTable_Delete
        procedure, non_overridable, public :: IsPresent     => ProcedureHashTable_IsPresent
        procedure, non_overridable, public :: Length        => ProcedureHashTable_Length
        procedure, non_overridable, public :: Print         => ProcedureHashTable_Print
        procedure, non_overridable, public :: Free          => ProcedureHashTable_Free
        final                              ::                  ProcedureHashTable_Finalize
    end type

public :: ProcedureHashTable_t

contains


    function ProcedureHashTable_Hash(this,Key) result(Hash)
    !-----------------------------------------------------------------
    !< String hash function
    !-----------------------------------------------------------------
        class(ProcedureHashTable_t), intent(IN) :: this               !< Procedure Hash Table
        character(len=*),            intent(IN) :: Key                !< String Key
        integer                                 :: Hash               !< Hash code
        character, dimension(len(Key))          :: CharArray          !< Character array containing the Key
        integer                                 :: CharIterator       !< Char iterator index
    !-----------------------------------------------------------------
        forall (CharIterator=1:LEN(Key))
            CharArray(CharIterator) = Key(CharIterator:CharIterator)
        end forall
        Hash = MOD(SUM(ICHAR(CharArray)), this%Size)
    end function ProcedureHashTable_Hash


    subroutine ProcedureHashTable_Init(this,Size)
    !-----------------------------------------------------------------
    !< Allocate the database with a given Szie of DefaultDataBaseSize
    !-----------------------------------------------------------------
        class(ProcedureHashTable_t), intent(INOUT) :: this            !< Procedure Hash Table
        integer, optional,           intent(IN)    :: Size            !< DataBase Size
    !-----------------------------------------------------------------
        call this%Free()
        if (present(Size)) then
            this%Size = Size
        else
            this%Size = DefaultDataBaseSize
        endif
        allocate(this%DataBase(0:this%Size-1))
    end subroutine ProcedureHashTable_Init


    function ProcedureHashTable_isInitialized(this) result(isInitialized)
    !-----------------------------------------------------------------
    !< Check if a the hash table is allocated
    !-----------------------------------------------------------------
        class(ProcedureHashTable_t), intent(IN) :: this               !< Procedure Hash Table
        logical                                 :: isInitialized      !< Boolean flag to check if the database is initialized
    !-----------------------------------------------------------------
        isInitialized = allocated(this%DataBase)
    end function ProcedureHashTable_isInitialized


    function ProcedureHashTable_isPresent(this,Key) result(isPresent)
    !-----------------------------------------------------------------
    !< Check if a Key is present in the DataBase
    !-----------------------------------------------------------------
        class(ProcedureHashTable_t), intent(IN) :: this               !< Procedure Hash Table
        character(len=*),            intent(IN) :: Key                !< String Key
        logical                                 :: isPresent          !< Boolean flag to check if a Key is present
    !-----------------------------------------------------------------
        isPresent = this%DataBase(this%Hash(Key=Key))%isPresent(Key=Key)
    end function ProcedureHashTable_isPresent


    subroutine ProcedureHashTable_Set(this,Key,Proc)
    !-----------------------------------------------------------------
    !< Set a Key/Procedure pair into the DataBase
    !-----------------------------------------------------------------
        class(ProcedureHashTable_t), intent(INOUT) :: this            !< Procedure Hash Table
        character(len=*),            intent(IN)    :: Key             !< String Key
        procedure(),      pointer                  :: Proc            !< Procedure pointer
    !-----------------------------------------------------------------
        call this%DataBase(this%Hash(Key=Key))%AddEntry(Key=Key,Proc=Proc)
    end subroutine ProcedureHashTable_Set


    subroutine ProcedureHashTable_Get(this,Key,Proc)
    !-----------------------------------------------------------------
    !< Return a Procedure given the Key
    !-----------------------------------------------------------------
        class(ProcedureHashTable_t), intent(IN) :: this               !< Procedure Hash Table
        character(len=*),            intent(IN) :: Key                !< String Key
        procedure(),                pointer     :: Proc               !< Returned procedure pointer
        class(ProcedureListNode_t), pointer     :: Entry              !< Pointer to a Parameter List
    !-----------------------------------------------------------------
        Entry => this%DataBase(this%Hash(Key=Key))%GetEntry(Key=Key)
        if(associated(Entry)) call Entry%GetProcedure(Proc=Proc)
    end subroutine ProcedureHashTable_Get


    subroutine ProcedureHashTable_Delete(this, Key)
    !-----------------------------------------------------------------
    !< Remove an Entry given a Key
    !-----------------------------------------------------------------
        class(ProcedureHashTable_t), intent(INOUT) :: this            !< Procedure Hash Table
        character(len=*),            intent(IN)    :: Key             !< String Key
    !-----------------------------------------------------------------
        call this%DataBase(this%Hash(Key=Key))%RemoveEntry(Key=Key)
    end subroutine ProcedureHashTable_Delete


    function ProcedureHashTable_Length(this) result(Length)
    !-----------------------------------------------------------------
    !< Return the number of Nodes contained in the DataBase
    !-----------------------------------------------------------------
        class(ProcedureHashTable_t), intent(IN) :: this               !< Procedure Hash Table
        integer                                 :: Length             !< Number of parameters in database
        integer                                 :: DBIterator         !< Database Iterator index 
    !-----------------------------------------------------------------
        Length = 0
        if (allocated(this%DataBase)) THEN
            do DBIterator=lbound(this%DataBase,dim=1),ubound(this%DataBase,dim=1)
                    Length = Length + this%DataBase(DBIterator)%Length()
            enddo
        endif
    end function ProcedureHashTable_Length


    subroutine ProcedureHashTable_Free(this)
    !-----------------------------------------------------------------
    !< Free the DataBase and its nodes
    !-----------------------------------------------------------------
        class(ProcedureHashTable_t), intent(INOUT) :: this            !< Procedure Hash Table
        integer                                    :: DBIterator      !< Database Iterator index 
    !-----------------------------------------------------------------
        if (allocated(this%DataBase)) THEN
            do DBIterator=lbound(this%DataBase,dim=1),ubound(this%DataBase,dim=1)
                call this%DataBase(DBIterator)%Free()
            enddo
            deallocate(this%DataBase)
        endif
        this%Size = 0    
    end subroutine ProcedureHashTable_Free


    subroutine ProcedureHashTable_Finalize(this)
    !-----------------------------------------------------------------
    !< Destructor procedure
    !-----------------------------------------------------------------
        type(ProcedureHashTable_t), intent(INOUT) :: this             !< Procedure Hash Table
    !-----------------------------------------------------------------
        call this%Free()
    end subroutine ProcedureHashTable_Finalize


    subroutine ProcedureHashTable_Print(this, unit, prefix, iostat, iomsg)
    !-----------------------------------------------------------------
    !< Print the content of the DataBase
    !-----------------------------------------------------------------
        class(ProcedureHashTable_t),           intent(IN)  :: this    !< Linked List
        integer,                               intent(IN)  :: unit    !< Logic unit.
        character(*), optional,                intent(IN)  :: prefix  !< Prefixing string.
        integer,      optional,                intent(OUT) :: iostat  !< IO error.
        character(*), optional,                intent(OUT) :: iomsg   !< IO error message.
        character(len=:), allocatable                      :: prefd   !< Prefixing string.
        integer                                            :: iostatd !< IO error.
        character(500)                                     :: iomsgd  !< Temporary variable for IO error message.
        integer                                            :: DBIter  !< Database iterator
        character(64)                                      :: iterchar!< aux string variable
    !-----------------------------------------------------------------
        prefd = '' ; if (present(prefix)) prefd = prefix
        if (allocated(this%DataBase)) then
            do DBIter=lbound(this%DataBase,dim=1), ubound(this%DataBase,dim=1)
                write(iterchar,*) DBIter
                call this%DataBase(DBIter)%Print(unit=unit,                         &
                    prefix=prefd//'  ['//trim(adjustl(iterchar))//'] ', &
                    iostat=iostatd,iomsg=iomsgd)
            enddo
        endif
        if (present(iostat)) iostat = iostatd
        if (present(iomsg))  iomsg  = iomsgd
    end subroutine ProcedureHashTable_Print


end module ProcedureHashTable
