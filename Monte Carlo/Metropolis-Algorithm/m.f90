program try_loop_random
    implicit none

    call random_ge
end program try_loop_random

subroutine xc(x, y, z, r)
    real:: x,y,z,r

    z = x + y + r
    print *, "Nilai z = ", z

end 

subroutine random_ge
    integer::i
    real:: x, y, z, r
    x = 2
    y = 1


    do i = 1, 10
        call random_number(r)
        print *, "nilai random = ",r

        call xc(x, y, z, r)
    end do


end subroutine