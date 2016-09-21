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
!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
module sparse_matrix_utils_names

  USE types_names

implicit none

    !---------------------------------------------------------------------
    !< AUX PROCEDURES INTERFACES
    !---------------------------------------------------------------------

    interface
        subroutine duplicates_operation(input, output)
            import rp
            real(rp), intent(in)    :: input
            real(rp), intent(inout) :: output
        end subroutine duplicates_operation
    end interface

contains

!---------------------------------------------------------------------
!< AUX PROCEDURES
!---------------------------------------------------------------------

    subroutine assign_value(input, output)
    !-----------------------------------------------------------------
    ! Assign an input value to the autput
    !-----------------------------------------------------------------
        real(rp), intent(in)    :: input
        real(rp), intent(inout) :: output
    !-------------------------------------------------------------
        output = input
    end subroutine assign_value


    subroutine sum_value(input, output)
    !-----------------------------------------------------------------
    ! Sum an input value to the autput
    !-----------------------------------------------------------------
        real(rp), intent(in)    :: input
        real(rp), intent(inout) :: output
    !-----------------------------------------------------------------
        output = output + input
    end subroutine sum_value


    function  binary_search(key,size,vector) result(ipos)
    !-----------------------------------------------------------------
    ! Perform binary search in a vector
    !-----------------------------------------------------------------
        integer(ip), intent(in) :: key
        integer(ip), intent(in) :: size
        integer(ip), intent(in) :: vector(size)
        integer(ip)             :: ipos
        integer(ip)             :: lowerbound
        integer(ip)             :: upperbound
        integer(ip)             :: midpoint
    !-----------------------------------------------------------------
        lowerbound = 1 
        upperbound = size
        ipos = -1 

        do while (lowerbound.le.upperbound) 
            midpoint = (lowerbound+upperbound)/2
            if (key.eq.vector(midpoint))  then
                ipos = midpoint
                lowerbound  = upperbound + 1
            else if (key < vector(midpoint))  then
                upperbound = midpoint-1
            else 
                lowerbound = midpoint + 1
            end if
        enddo
        return
    end function binary_search


    subroutine mergesort_link_list(n,k,l,iret)
    !-------------------------------------------------------------
    !   This subroutine sorts an integer array into ascending order.
    !
    ! Arguments:
    !   n    -  integer           Input: size of the array 
    !   k    -  integer(*)        input: array of keys to be sorted
    !   l    -  integer(0:n+1)   output: link list 
    !   iret -  integer          output: 0 Normal termination
    !                                    1 the array was already sorted 
    !
    ! REFERENCES  = (1) D. E. Knuth
    !                   The Art of Computer Programming,
    !                     vol.3: Sorting and Searching
    !                   Addison-Wesley, 1973
    !-------------------------------------------------------------
        use types_names
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: k(n)
        integer(ip), intent(out) :: l(0:n+1)
        integer(ip), intent(out) :: iret
        integer(ip)              :: p,q,s,t
    !-------------------------------------------------------------
        iret = 0
        !  first step: we are preparing ordered sublists, exploiting
        !  what order was already in the input data; negative links
        !  mark the end of the sublists
        l(0) = 1
        t = n + 1
        do  p = 1,n - 1
        if (k(p) <= k(p+1)) then
            l(p) = p + 1
            else
                l(t) = - (p+1)
                t = p
            end if
        end do
        l(t) = 0
        l(n) = 0
        ! see if the input was already sorted
        if (l(n+1) == 0) then
            iret = 1
            return 
        else
            l(n+1) = abs(l(n+1))
        end if

        mergepass: do 
            ! otherwise, begin a pass through the list.
            ! throughout all the subroutine we have:
            !  p, q: pointing to the sublists being merged
            !  s: pointing to the most recently processed record
            !  t: pointing to the end of previously completed sublist
            s = 0
            t = n + 1
            p = l(s)
            q = l(t)
            if (q == 0) exit mergepass

            outer: do 

                if (k(p) > k(q)) then 
                    l(s) = sign(q,l(s))
                    s = q
                    q = l(q)
                    if (q > 0) then 
                        do 
                            if (k(p) <= k(q)) cycle outer
                            s = q
                            q = l(q)
                            if (q <= 0) exit
                        end do
                    end if
                    l(s) = p
                    s = t
                    do 
                        t = p
                        p = l(p)
                        if (p <= 0) exit
                    end do

                else 
                    l(s) = sign(p,l(s))
                    s = p
                    p = l(p)
                    if (p>0) then 
                        do 
                            if (k(p) > k(q)) cycle outer 
                            s = p
                            p = l(p)
                            if (p <= 0) exit
                        end do
                    end if
                    !  otherwise, one sublist ended, and we append to it the rest
                    !  of the other one.
                    l(s) = q
                    s = t
                    do 
                        t = q
                        q = l(q)
                        if (q <= 0) exit
                    end do
                end if

                p = -p
                q = -q
                if (q == 0) then
                    l(s) = sign(p,l(s))
                    l(t) = 0
                    exit outer 
                end if
            end do outer
        end do mergepass

    end subroutine mergesort_link_list


    subroutine reorder_ip_from_link_list(n,i1,iaux)
    !-------------------------------------------------------------
    !  Reorder (an) input vector(s) based on a list sort output.
    !  Based on: D. E. Knuth: The Art of Computer Programming
    !            vol. 3: Sorting and Searching, Addison Wesley, 1973
    !            ex. 5.2.12
    !-------------------------------------------------------------
        use types_names
        integer(ip), intent(in)    :: n
        integer(ip), intent(inout) :: i1(*)
        integer(ip), intent(inout) :: iaux(0:*) 
        integer(ip) :: lswap, lp, k, isw1
    !-------------------------------------------------------------
        lp = iaux(0)
        k  = 1
        do 
            if ((lp == 0).or.(k>n)) exit
            do 
                if (lp >= k) exit
                lp = iaux(lp)
            end do
            isw1     = i1(lp)
            i1(lp)   = i1(k)
            i1(k)    = isw1
            lswap    = iaux(lp)
            iaux(lp) = iaux(k)
            iaux(k)  = lp
            lp = lswap 
            k  = k + 1
        enddo
        return
    end subroutine reorder_ip_from_link_list

    
    subroutine reorder_ip_rp_from_link_list(n,x,i1,iaux)
    !-------------------------------------------------------------
    !  Reorder (an) input vector(s) based on a list sort output.
    !  Based on: D. E. Knuth: The Art of Computer Programming
    !            vol. 3: Sorting and Searching, Addison Wesley, 1973
    !            ex. 5.2.12
    !-------------------------------------------------------------
        use types_names
        integer(ip), intent(in)    :: n
        real(rp),    intent(inout) :: x(*)
        integer(ip), intent(inout) :: i1(*)
        integer(ip), intent(inout) :: iaux(0:*) 
        integer(ip) :: lswap, lp, k, isw1
        real(rp)    :: swap
    !-------------------------------------------------------------
        lp = iaux(0)
        k  = 1
        do 
            if ((lp == 0).or.(k>n)) exit
            do 
                if (lp >= k) exit
                lp = iaux(lp)
            end do
            swap     = x(lp)
            x(lp)    = x(k)
            x(k)     = swap
            isw1     = i1(lp)
            i1(lp)   = i1(k)
            i1(k)    = isw1
            lswap    = iaux(lp)
            iaux(lp) = iaux(k)
            iaux(k)  = lp
            lp = lswap 
            k  = k + 1
        enddo
        return
    end subroutine reorder_ip_rp_from_link_list


    subroutine reorder_numeric_coo_from_link_list(n,x,i1,i2,iaux)
    !-------------------------------------------------------------
    !  Reorder (an) input vector(s) based on a list sort output.
    !  Based on: D. E. Knuth: The Art of Computer Programming
    !            vol. 3: Sorting and Searching, Addison Wesley, 1973
    !            ex. 5.2.12
    !-------------------------------------------------------------
        use types_names
        integer(ip), intent(in)    :: n
        real(rp),    intent(inout) :: x(*)
        integer(ip), intent(inout) :: i1(*)
        integer(ip), intent(inout) :: i2(*)
        integer(ip), intent(inout) :: iaux(0:*) 
        integer(ip) :: lswap, lp, k, isw1, isw2
        real(rp)    :: swap
    !-------------------------------------------------------------
        lp = iaux(0)
        k  = 1
        do 
            if ((lp == 0).or.(k>n)) exit
            do 
                if (lp >= k) exit
                lp = iaux(lp)
            end do
            swap     = x(lp)
            x(lp)    = x(k)
            x(k)     = swap
            isw1     = i1(lp)
            i1(lp)   = i1(k)
            i1(k)    = isw1
            isw2     = i2(lp)
            i2(lp)   = i2(k)
            i2(k)    = isw2
            lswap    = iaux(lp)
            iaux(lp) = iaux(k)
            iaux(k)  = lp
            lp = lswap 
            k  = k + 1
        enddo
        return
    end subroutine reorder_numeric_coo_from_link_list


    subroutine reorder_symbolic_coo_from_link_list(n,i1,i2,iaux)
    !-------------------------------------------------------------
    !  Reorder (an) input vector(s) based on a list sort output.
    !  Based on: D. E. Knuth: The Art of Computer Programming
    !            vol. 3: Sorting and Searching, Addison Wesley, 1973
    !            ex. 5.2.12
    !-------------------------------------------------------------
        use types_names
        integer(ip), intent(in)    :: n
        integer(ip), intent(inout) :: i1(*)
        integer(ip), intent(inout) :: i2(*)
        integer(ip), intent(inout) :: iaux(0:*) 
        integer(ip) :: lswap, lp, k, isw1, isw2
    !-------------------------------------------------------------
        lp = iaux(0)
        k  = 1
        do 
            if ((lp == 0).or.(k>n)) exit
            do 
                if (lp >= k) exit
                lp = iaux(lp)
            end do
            isw1     = i1(lp)
            i1(lp)   = i1(k)
            i1(k)    = isw1
            isw2     = i2(lp)
            i2(lp)   = i2(k)
            i2(k)    = isw2
            lswap    = iaux(lp)
            iaux(lp) = iaux(k)
            iaux(k)  = lp
            lp = lswap 
            k  = k + 1
        enddo
        return
    end subroutine reorder_symbolic_coo_from_link_list

end module sparse_matrix_utils_names
