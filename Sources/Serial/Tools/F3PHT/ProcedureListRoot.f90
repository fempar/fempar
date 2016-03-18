!-----------------------------------------------------------------
! F3PHT (Fortran Polymorphic Procedure Pointer Hash Table)
! Copyright (c) 2016 Santiago Badia, Alberto F. Martín, 
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

module ProcedureListRoot

USE ProcedureListNode

implicit none
private

    type :: ProcedureListRoot_t
    private
        class(ProcedureListNode_t), pointer :: Root => null()
    contains
    private
        procedure, non_overridable        :: Init             => ProcedureListRoot_Init
        procedure, non_overridable         :: HasRoot          => ProcedureListRoot_HasRoot
        procedure, non_overridable         :: SetRoot          => ProcedureListRoot_SetRoot
        procedure, non_overridable         :: GetRoot          => ProcedureListRoot_GetRoot
        procedure, non_overridable         :: NullifyRoot      => ProcedureListRoot_NullifyRoot
        procedure, non_overridable         :: DeallocateRoot   => ProcedureListRoot_DeallocateRoot
        procedure, non_overridable, public :: GetEntry         => ProcedureListRoot_GetEntry
        procedure, non_overridable, public :: GetPreviousEntry => ProcedureListRoot_GetPreviousEntry
        procedure, non_overridable, public :: Print            => ProcedureListRoot_Print
        procedure, non_overridable, public :: isPresent        => ProcedureListRoot_isPresent
        procedure, non_overridable, public :: Length           => ProcedureListRoot_Length
        procedure, non_overridable, public :: RemoveEntry      => ProcedureListRoot_RemoveEntry
        procedure, non_overridable, public :: AddEntry         => ProcedureListRoot_AddEntry
        procedure, non_overridable, public :: Free             => ProcedureListRoot_Free
        final                              ::                     ProcedureListRoot_Finalize 
    end type


public :: ProcedureListRoot_t

contains


    subroutine ProcedureListRoot_SetRoot(this, Root)
    !-----------------------------------------------------------------
    !< Set the Root of the list
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t),          intent(INOUT) :: this    !< Procedure List Root
        class(ProcedureListNode_t), pointer, intent(IN)    :: Root    !< Parameter Entry correspoing to the head of the list
    !-----------------------------------------------------------------
        this%Root => Root
    end subroutine ProcedureListRoot_SetRoot


    function ProcedureListRoot_GetRoot(this) result(Root)
    !-----------------------------------------------------------------
    !< Return a pointer to the Root of the list
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t),         intent(IN) :: this        !< Procedure List Root
        class(ProcedureListNode_t), pointer            :: Root        !< Parameter Entry correspoing to the head of the list
    !-----------------------------------------------------------------
        Root => this%Root
    end function ProcedureListRoot_GetRoot


    function ProcedureListRoot_HasRoot(this) result(HasRoot)
    !-----------------------------------------------------------------
    !< Return a pointer to the Root of the list
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t),         intent(IN) :: this        !< Procedure List Root
        logical                                        :: hasRoot     !< Check if Root is associated
    !-----------------------------------------------------------------
        hasRoot = associated(this%GetRoot())
    end function ProcedureListRoot_HasRoot


    subroutine ProcedureListRoot_NullifyRoot(this)
    !-----------------------------------------------------------------
    !< Set the Root of the list
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t),   intent(INOUT) :: this           !< Procedure List Root
    !-----------------------------------------------------------------
        nullify(this%Root)
    end subroutine ProcedureListRoot_NullifyRoot


    subroutine ProcedureListRoot_DeallocateRoot(this)
    !-----------------------------------------------------------------
    !< Set the Root of the list
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t),   intent(INOUT) :: this           !< Procedure List Root
    !-----------------------------------------------------------------
        if(this%HasRoot()) then
            call this%Root%Free()
            deallocate(this%Root)
        endif
    end subroutine ProcedureListRoot_DeallocateRoot


    subroutine ProcedureListRoot_Init(this, Key, Proc)
    !-----------------------------------------------------------------
    !< Initialize the Root of the list
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t), intent(INOUT) :: this             !< Procedure List Root
        character(len=*),           intent(IN)    :: Key              !< Key (unique) of the current node.
        procedure(), pointer                      :: Proc             !< Procedure pointer
    !-----------------------------------------------------------------
        if(.not. this%HasRoot()) allocate(ProcedureListNode_t::this%Root)
        call this%Root%SetKey(Key=Key)
        call this%Root%SetProcedure(Proc=Proc)
    end subroutine ProcedureListRoot_Init


    function ProcedureListRoot_IsPresent(this, Key) result(isPresent)
    !-----------------------------------------------------------------
    !< Check if a Key is present in the List
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t), intent(IN)  :: this               !< Procedure List Root
        character(len=*),           intent(IN)  :: Key                !< String Key
        logical                                 :: isPresent          !< Boolean flag to check if a Key is present
    !-----------------------------------------------------------------
        isPresent = associated(this%GetEntry(Key))
    end function ProcedureListRoot_IsPresent


    subroutine ProcedureListRoot_AddEntry(this,Key, Proc)
    !-----------------------------------------------------------------
    !< Add a new Node if key does not Exist
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t),         intent(INOUT) :: this         !< Procedure List Root
        character(len=*),                   intent(IN)    :: Key          !< Key (unique) of the current node.
        procedure(),                pointer               :: Proc         !< Procedure pointer
        class(ProcedureListNode_t), pointer               :: NextEntry    !< Parameter Entry
        class(ProcedureListNode_t), pointer               :: NewEntry     !< New Parameter Entry
        character(len=:), allocatable                     :: NextEntryKey !< Key of the NextEntry
    !-----------------------------------------------------------------
        if(.not. this%HasRoot()) then
            call this%Init(Key=Key, Proc=Proc)
        else
            NextEntry => this%GetRoot()
            do while(associated(NextEntry))
                NExtEntryKey = NextEntry%GetKey()
                if (NextEntryKey/=Key) then
                    if (.not. NextEntry%hasNext()) then 
                        ! I reached the end of the list
                        allocate(ProcedureListNode_t::NewEntry)
                        call NewEntry%SetKey(Key=Key)
                        call NewEntry%SetProcedure(Proc=Proc)
                        call NextEntry%SetNext(NExt=NewEntry)
                        exit
                    else
                        NextEntry => NextEntry%GetNext()
                    endif
                else
                    call NextEntry%SetProcedure(Proc=Proc)
                    exit
                endif
            enddo
            if(allocated(NextEntryKey)) deallocate(NextEntryKey)
        endif
    end subroutine ProcedureListRoot_AddEntry


    subroutine ProcedureListRoot_RemoveEntry(this, Key)
    !-----------------------------------------------------------------
    !< Remove an Entry given a Key
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t),           intent(INOUT) :: this            !< Procedure List Root
        character(len=*),                     intent(IN)    :: Key             !< String Key
        character(len=:), allocatable                       :: CurrentEntryKey !< Current Entry Key
        class(ProcedureListNode_t),     pointer             :: PreviousEntry   !< The Previous Entry of a given key
        class(ProcedureListNode_t),     pointer             :: CurrentEntry    !< Entry of a given key
        class(ProcedureListNode_t),     pointer             :: NextEntry       !< The Next Entry of a given key
    !-----------------------------------------------------------------
        if(this%HasRoot()) then
            CurrentEntry => this%GetRoot()
            CurrentEntryKey = CurrentEntry%GetKey()
            if(CurrentEntryKey == Key) then
                NextEntry => CurrentEntry%GetNext()
                call CurrentEntry%DeallocateKey()    
                call CurrentEntry%NullifyProcedure()
                call CurrentEntry%NullifyNext()
                deallocate(CurrentEntry)
                call this%NullifyRoot()
                if(allocated(CurrentEntryKey)) deallocate(CurrentEntryKey)
            else
                PreviousEntry     => this%GetPreviousEntry(Key=Key)
                if(associated(PreviousEntry)) then
                    CurrentEntry  => PreviousEntry%GetNext()
                    NextEntry     => CurrentEntry%GetNext()
                    call CurrentEntry%DeallocateKey()    
                    call CurrentEntry%NullifyProcedure()
                    call CurrentEntry%NullifyNext()
                    deallocate(CurrentEntry)
                    call PreviousEntry%NullifyNext()
                    if(associated(NextEntry)) call PreviousEntry%SetNext(Next=NextEntry)
                endif   
            endif
            if(associated(NextEntry)) call this%SetRoot(Root = NextEntry)
        endif
    end subroutine ProcedureListRoot_RemoveEntry



    function ProcedureListRoot_GetEntry(this,Key) result(Entry)
    !-----------------------------------------------------------------
    !< Return a pointer to a ParameterEntry given a Key
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t),            intent(IN) :: this     !< Procedure List Root
        character(len=*),                      intent(IN) :: Key      !< String Key
        class(ProcedureListNode_t),     pointer           :: Entry    !< Parameter Entry
        character(len=:), allocatable                     :: EntryKey !< Entry Key
    !-----------------------------------------------------------------
        Entry => this%GetRoot()
        do while(associated(Entry))
            EntryKey = Entry%GetKey()
            if (EntryKey==Key) exit
            Entry => Entry%GetNext()
        enddo
        if(allocated(EntryKey)) deallocate(EntryKey)
    end function ProcedureListRoot_GetEntry


    function ProcedureListRoot_GetPreviousEntry(this,Key) result(PreviousEntry)
    !-----------------------------------------------------------------
    !< Return a pointer to the provious node of a Parameter List given a Key
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t),        intent(IN) :: this         !< Parameter List
        character(len=*),                  intent(IN) :: Key          !< String Key
        class(ProcedureListNode_t), pointer           :: PreviousEntry!< Parameter List Entry
        class(ProcedureListNode_t), pointer           :: NextEntry    !< Parameter List Next Entry
        character(len=:), allocatable                 :: NExtEntryKey !< NextEntry Key
    !-----------------------------------------------------------------
        PreviousEntry => this%GetRoot()
        do while(associated(PreviousEntry))
            if (PreviousEntry%HasNext()) then
                NextEntry => PreviousEntry%GetNext()
                NextEntryKey = NextEntry%GetKey()
                if (NextEntryKey==Key) then
                    exit
                else
                    PreviousEntry => NextEntry
                endif
            else
                nullify(PreviousEntry)
                exit
            endif
        enddo    
        if(allocated(NextEntryKey)) deallocate(NextEntryKey)
    end function ProcedureListRoot_GetPreviousEntry


    function ProcedureListRoot_Length(this) result(Length)
    !-----------------------------------------------------------------
    !< Return the length of the list
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t),       intent(IN) :: this          !< Procedure List Root
        integer                                      :: Length        !< Length of the list
        type(ProcedureListNode_t), pointer           :: NextEntry     !< Next Parameter Entry
    !-----------------------------------------------------------------
        Length = 0
        NextEntry => this%GetRoot()
        do while (associated(NextEntry))
            Length = Length + 1
            NextEntry => NextEntry%GetNext()
        enddo
        nullify(NextEntry)
    end function ProcedureListRoot_Length



    subroutine ProcedureListRoot_Free(this)
    !-----------------------------------------------------------------
    !< Free the list
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t),         intent(INOUT) :: this     !< Procedure List Root
        class(ProcedureListNode_t), pointer               :: Current  !< Current Parameter List Node
        class(ProcedureListNode_t), pointer               :: Next     !< Next Parameter List Node
    !-----------------------------------------------------------------
        do while(this%HasRoot()) 
            Next => this%Root%GetNext()
            call this%Root%Free()
            call this%DeallocateRoot()
            call this%SetRoot(Root=Next)
        enddo
    end subroutine ProcedureListRoot_Free



    subroutine ProcedureListRoot_Print(this, unit, prefix, iostat, iomsg)
    !-----------------------------------------------------------------
    !< Print the keys/value pair contained in the parameter list
    !-----------------------------------------------------------------
        class(ProcedureListRoot_t),           intent(IN)  :: this     !< Procedure List Root
        integer,                              intent(IN)  :: unit     !< Logic unit.
        character(*),               optional, intent(IN)  :: prefix   !< Prefixing string.
        integer,                    optional, intent(OUT) :: iostat   !< IO error.
        character(*),               optional, intent(OUT) :: iomsg    !< IO error message.
        character(len=:), allocatable                     :: prefd    !< Prefixing string.
        integer                                           :: iostatd  !< IO error.
        character(500)                                    :: iomsgd   !< Temporary variable for IO error message.
        class(ProcedureListNode_t), pointer               :: NextEntry!< Pointer for scanning the list.
        integer                                           :: counter  !< List node counter
        character(64)                                     :: ch       !< List node counter in text
    !-----------------------------------------------------------------
        iostatd = 0 ; iomsgd = ''; prefd = '';if (present(prefix)) prefd = prefix
        if(this%HasRoot()) then
            counter = 0
            NextEntry => this%GetRoot()
            do while(associated(NextEntry))
                write(ch,*) counter
                call NextEntry%Print(unit=unit, prefix=prefd//'['//trim(adjustl(ch))//']', iostat=iostatd, iomsg=iomsgd )
                NextEntry => NextEntry%GetNext()
                counter = counter+1
            enddo
        endif
        if (present(iostat)) iostat = iostatd
        if (present(iomsg))  iomsg  = iomsgd
    end subroutine ProcedureListRoot_Print


    subroutine ProcedureListRoot_Finalize(this)
    !-----------------------------------------------------------------
    !< Finalize procedure
    !-----------------------------------------------------------------
        type(ProcedureListRoot_t), intent(INOUT):: this     !< Procedure List Root
    !-----------------------------------------------------------------
        call this%Free()
    end subroutine ProcedureListRoot_Finalize


end module ProcedureListRoot
