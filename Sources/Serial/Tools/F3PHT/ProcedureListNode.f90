!-----------------------------------------------------------------
! F3PHT (Fortran Polymorphic Procedure Pointer Hash Table)
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

module ProcedureListNode

implicit none
private

    type :: ProcedureListNode_t
    private
        character(len=:), allocatable               :: Key
        procedure(),                pointer, nopass :: Proc => null()
        class(ProcedureListNode_t), pointer         :: Next => null()
    contains
    private
        procedure, non_overridable, public :: Print            => ProcedureListNode_Print
        procedure, non_overridable, public :: HasNext          => ProcedureListNode_HasNext
        procedure, non_overridable, public :: SetNext          => ProcedureListNode_SetNext
        procedure, non_overridable, public :: GetNext          => ProcedureListNode_GetNext
        procedure, non_overridable, public :: NullifyNext      => ProcedureListNode_NullifyNext
        procedure, non_overridable, public :: HasKey           => ProcedureListNode_HasKey
        procedure, non_overridable, public :: SetKey           => ProcedureListNode_SetKey
        procedure, non_overridable, public :: GetKey           => ProcedureListNode_GetKey
        procedure, non_overridable, public :: DeallocateKey    => ProcedureListNode_DeallocateKey
        procedure, non_overridable, public :: HasProcedure     => ProcedureListNode_HasProcedure
        procedure, non_overridable, public :: SetProcedure     => ProcedureListNode_SetProcedure
        procedure, non_overridable, public :: GetProcedure     => ProcedureListNode_GetProcedure
        procedure, non_overridable, public :: NullifyProcedure => ProcedureListNode_NullifyProcedure
        procedure, non_overridable, public :: Free             => ProcedureListNode_Free
        final                              ::                     ProcedureListNode_Finalize 
    end type ProcedureListNode_t

public :: ProcedureListNode_t

contains


    function ProcedureListNode_HasNext(this) result(hasNext)
    !-----------------------------------------------------------------
    !< Check if Next is associated for the current Node
    !-----------------------------------------------------------------
        class(ProcedureListNode_t), intent(IN) :: this                !< Procedure List Node
        logical                                :: hasNext             !< Check if Next is associated
    !-----------------------------------------------------------------
        hasNext = associated(this%Next)
    end function ProcedureListNode_HasNext


    subroutine ProcedureListNode_SetNext(this, Next)
    !-----------------------------------------------------------------
    !< Set the pointer to the Next node
    !-----------------------------------------------------------------
        class(ProcedureListNode_t),          intent(INOUT) :: this    !< Procedure List Node
        class(ProcedureListNode_t), pointer, intent(IN)    :: Next    !< Pointer to Next 
    !-----------------------------------------------------------------
        this%Next => Next
    end subroutine ProcedureListNode_SetNext


    function ProcedureListNode_GetNext(this) result(Next)
    !-----------------------------------------------------------------
    !< Return a pointer to the Next node
    !-----------------------------------------------------------------
        class(ProcedureListNode_t), intent(IN) :: this                !< Procedure List Node
        class(ProcedureListNode_t), pointer    :: Next                !< Pointer to Next
    !-----------------------------------------------------------------
        nullify(Next)
        if(this%HasNext()) Next => this%Next
    end function ProcedureListNode_GetNext


    subroutine ProcedureListNode_NullifyNext(this)
    !-----------------------------------------------------------------
    !< Nullify Next
    !-----------------------------------------------------------------
        class(ProcedureListNode_t), intent(INOUT) :: this             !< Procedure List Node
    !-----------------------------------------------------------------
        nullify(this%Next)
    end subroutine ProcedureListNode_NullifyNext


    function ProcedureListNode_HasKey(this) result(hasKey)
    !-----------------------------------------------------------------
    !< Check if Key is allocated for the current Node
    !-----------------------------------------------------------------
        class(ProcedureListNode_t),     intent(IN) :: this            !< Procedure List Node
        logical                                    :: hasKey          !< Check if Key is associated
    !-----------------------------------------------------------------
        hasKey = allocated(this%Key)
    end function ProcedureListNode_HasKey


    subroutine ProcedureListNode_SetKey(this, Key) 
    !-----------------------------------------------------------------
    !< Check if Next is associated for the current Node
    !-----------------------------------------------------------------
        class(ProcedureListNode_t),      intent(INOUT) :: this        !< Procedure List Node
        character(len=*),                intent(IN)    :: Key         !< Key
    !-----------------------------------------------------------------
        this%Key = Key
    end subroutine ProcedureListNode_SetKey


    function ProcedureListNode_GetKey(this) result(Key)
    !-----------------------------------------------------------------
    !< Check if Next is associated for the current Node
    !-----------------------------------------------------------------
        class(ProcedureListNode_t),     intent(IN) :: this           !< Procedure List Node
        character(len=:), allocatable              :: Key            !< Key
    !-----------------------------------------------------------------
        Key = this%Key
    end function ProcedureListNode_GetKey


    subroutine ProcedureListNode_DeallocateKey(this)
    !-----------------------------------------------------------------
    !< Deallocate Key if allocated
    !-----------------------------------------------------------------
        class(ProcedureListNode_t), intent(INOUT) :: this             !< Procedure List Node
    !-----------------------------------------------------------------
        if(this%HasKey()) deallocate(this%Key)
    end subroutine ProcedureListNode_DeallocateKey


    subroutine ProcedureListNode_Free(this)
    !-----------------------------------------------------------------
    !< Free the Entry
    !-----------------------------------------------------------------
        class(ProcedureListNode_t), intent(INOUT) :: this             !< Procedure List Node
    !-----------------------------------------------------------------
        call this%DeallocateKey()
        call this%NullifyProcedure()
        call this%NullifyNext()
    end subroutine ProcedureListNode_Free


    function ProcedureListNode_HasProcedure(this) result(hasProcedure)
    !-----------------------------------------------------------------
    !< Check if Value is allocated for the current Node
    !-----------------------------------------------------------------
        class(ProcedureListNode_t), intent(IN) :: this                !< Procedure List Node
        logical                                :: hasProcedure        !< Check if Procedure is associated
    !-----------------------------------------------------------------
        hasProcedure = associated(this%Proc)
    end function ProcedureListNode_HasProcedure


    subroutine ProcedureListNode_SetProcedure(this, Proc)
    !-----------------------------------------------------------------
    !< Set the procedure pointer
    !-----------------------------------------------------------------
        class(ProcedureListNode_t), intent(INOUT)  :: this            !< Procedure List Node
        procedure(),       pointer                 :: Proc            !< Procedure pointer
    !-----------------------------------------------------------------
        this%Proc => Proc
    end subroutine ProcedureListNode_SetProcedure


    subroutine ProcedureListNode_GetProcedure(this, Proc)
    !-----------------------------------------------------------------
    !< Return a procedure pointer
    !-----------------------------------------------------------------
        class(ProcedureListNode_t), intent(IN)  :: this               !< Procedure List Node
        procedure(),       pointer              :: Proc               !< Procedure pointer
    !-----------------------------------------------------------------
        Proc => this%Proc
    end subroutine ProcedureListNode_GetProcedure


    subroutine ProcedureListNode_NullifyProcedure(this)
    !-----------------------------------------------------------------
    !< Deallocate Key if allocated
    !-----------------------------------------------------------------
        class(ProcedureListNode_t), intent(INOUT) :: this             !< Procedure List Node
    !-----------------------------------------------------------------
        nullify(this%Proc)
    end subroutine ProcedureListNode_NullifyProcedure


    subroutine ProcedureListNode_Finalize(this)
    !-----------------------------------------------------------------
    !< Finalize procedure
    !-----------------------------------------------------------------
        type(ProcedureListNode_t), intent(INOUT):: this               !< Procedure List Node
    !-----------------------------------------------------------------
        call this%Free()
    end subroutine ProcedureListNode_Finalize


    subroutine ProcedureListNode_Print(this, unit, prefix, iostat, iomsg)
    !-----------------------------------------------------------------
    !< Print the keys/value pair contained in the Procedure Hash Table List
    !-----------------------------------------------------------------
        class(ProcedureListNode_t),       intent(IN)  :: this         !< Procedure List Node
        integer,                          intent(IN)  :: unit         !< Logic unit.
        character(*), optional,           intent(IN)  :: prefix       !< Prefixing string.
        integer,      optional,           intent(OUT) :: iostat       !< IO error.
        character(*), optional,           intent(OUT) :: iomsg        !< IO error message.
        character(len=:),       allocatable           :: prefd        !< Prefixing string.
        integer                                       :: iostatd      !< IO error.
        character(500)                                :: iomsgd       !< Temporary variable for IO error message.
    !-----------------------------------------------------------------
        iostatd = 0 ; iomsgd = ''; prefd = '';if (present(prefix)) prefd = prefix
        if(this%HasKey()) then
            write(unit=unit,fmt='(A,L1)',iostat=iostatd,iomsg=iomsgd)prefd//' Key = "'//this%GetKey()//'", Procedure = ',this%HasProcedure()
        endif
        if (present(iostat)) iostat = iostatd
        if (present(iomsg))  iomsg  = iomsgd
    end subroutine ProcedureListNode_Print


end module ProcedureListNode
