program test_list_iterators

USE IR_precision, only: str
USE types_names
USE memor_names
USE list_types_names

implicit none

# include "debug.i90"

    type(list_t)             :: list
    type(list_iterator_t)    :: iterator
    integer(ip)              :: i, j
    integer(ip)              :: np = 5
    integer(ip)              :: nl = 15
    integer(ip), allocatable :: p(:)
    integer(ip), allocatable :: l(:)


    call meminit()


    call memalloc(np + 1, p, __FILE__, __LINE__)
    call memalloc(nl, l, __FILE__, __LINE__)

    p = (/1,3,3,5,10,16/)
    l = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)

    call list%create(np)
    call list%sum_to_pointer_index(index=1, value=2)
    call list%sum_to_pointer_index(index=3, value=2)
    call list%sum_to_pointer_index(index=4, value=5)
    call list%sum_to_pointer_index(index=5, value=6)
    call list%calculate_header()
    call list%allocate_list_from_pointer()
    iterator = list%create_iterator()
    i = 0

print*, '!----------------------------------------------------------------- '
print*, '!< FULL ITERATOR INCREASING LOOP - set_current'
print*, '!----------------------------------------------------------------- '
    do while(.not. iterator%is_upper_bound())
        i = i+1
        call iterator%set_current(l(i))
        assert(iterator%get_current() == l(i))
        call iterator%next()
    enddo

    call list%print(6)

print*, '!----------------------------------------------------------------- '
print*, '!< INDEX ITERATOR INCREASING LOOP - get_current'
print*, '!----------------------------------------------------------------- '

    j = 0
    do i=1, list%get_num_pointers()
        iterator = list%create_iterator(i)
        print*, 'Pointer: ', trim(str(no_sign=.true., n=i)), &
                ', components:', trim(str(no_sign=.true., n=iterator%get_size()))
        do while(.not. iterator%is_upper_bound())
            print*, '  current item: ', trim(str(no_sign=.true., n=iterator%get_current())), &
                    ', remaining items in iterator:', trim(str(no_sign=.true., n=iterator%get_distance_to_upper_bound()))
            j=j+1
            assert(iterator%get_current()==l(j))
            call iterator%next()
        enddo
    enddo

print*, '!----------------------------------------------------------------- '
print*, '!< INDEX ITERATOR DECREASING LOOP - get_current'
print*, '!----------------------------------------------------------------- '

    j = nl+1
    do i=list%get_num_pointers(),1,-1
        iterator = list%create_iterator(i)
        call iterator%end()
        print*, 'Pointer: ', trim(str(no_sign=.true., n=i)), &
                ', components:', trim(str(no_sign=.true., n=iterator%get_size()))
        do while(.not. iterator%is_lower_bound())
            print*, '  current item: ', trim(str(no_sign=.true., n=iterator%get_current())), &
                    ', remaining items in iterator:', trim(str(no_sign=.true., n=iterator%get_distance_to_lower_bound()))
            j=j-1
            assert(iterator%get_current()==l(j))
            call iterator%previous()
        enddo
    enddo

print*, '!----------------------------------------------------------------- '
print*, '!< RANGE ITERATOR INCREASING LOOP - get_current'
print*, '!----------------------------------------------------------------- '

    j = p(2)
    iterator = list%create_iterator(2,4)
    print*, 'Pointer from: ', trim(str(no_sign=.true., n=2)), ' to: ',trim(str(no_sign=.true., n=4)), &
            ', components:', trim(str(no_sign=.true., n=iterator%get_size()))
    do while(.not. iterator%is_upper_bound())
        print*, '  current item: ', trim(str(no_sign=.true., n=iterator%get_current())), &
                ', remaining items in iterator:', trim(str(no_sign=.true., n=iterator%get_distance_to_upper_bound()))
        assert(iterator%get_current() == l(j))
        j = j+1
        call iterator%next()
    enddo

print*, '!----------------------------------------------------------------- '
print*, '!< REACH SOME COMPONENTS FROM CURRENT POSITION'
print*, '!----------------------------------------------------------------- '
    iterator = list%create_iterator(4)
    print*, 'Pointer: ', trim(str(no_sign=.true., n=4)), &
            ', components:', trim(str(no_sign=.true., n=iterator%get_size()))
    do i=iterator%get_distance_to_upper_bound()-1, 1, -1
        print*, '  current item: ', trim(str(no_sign=.true., n=iterator%get_current())), &
                ', remaining items in iterator:', trim(str(no_sign=.true., n=iterator%get_distance_to_upper_bound())), &
                ', reached item:', iterator%get_from_current(i)
        assert(l(p(4)+i) == iterator%get_from_current(i))
    enddo

    call list%reinit(np)
    call list%sum_to_pointer_index(index=1, value=2)
    call list%sum_to_pointer_index(index=3, value=2)
    call list%sum_to_pointer_index(index=4, value=5)
    call list%sum_to_pointer_index(index=5, value=6)
    call list%calculate_header()
    call list%allocate_list_from_pointer()
    iterator = list%create_iterator()
    i = 0

print*, '!----------------------------------------------------------------- '
print*, '!< FULL ITERATOR INCREASING LOOP - set_current (after %reinit(n) )'
print*, '!----------------------------------------------------------------- '
    do while(.not. iterator%is_upper_bound())
        i = i+1
        call iterator%set_current(l(i))
        assert(iterator%get_current() == l(i))
        call iterator%next()
    enddo

    call list%print(6)

    call memfree(p, __FILE__, __LINE__)
    call memfree(l, __FILE__, __LINE__)
    call list%free()
    call memstatus()

end program test_list_iterators
