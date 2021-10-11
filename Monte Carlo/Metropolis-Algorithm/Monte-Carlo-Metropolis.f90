program Montecarlo_Metropolis    
    implicit none
    integer::L_size = 8, J=1
    real, dimension(8,8)::n

    call initialize(L_size, n)
    call claster(n, L_size, J)    

end program Montecarlo_Metropolis

subroutine initialize(L, n)        !L = Lattice Size
    implicit none
    integer::L, row, col
    real, dimension(L,L)::lattice, n 
    real::r   

    lattice = 1
    do col=1, L    
        do row = 1, L
            call random_seed()
            call random_number(r)
            if(r.gt.0.5) then
                lattice(col,row) = -1.0
            end if
        end do
        
    end do

    n = lattice
    write(*, 1) lattice
    1       format(8f6.1)

end subroutine initialize

subroutine claster(gl, L, j)
    implicit none
    real, dimension(L,L)::gl
    integer::x,y,top,bottom,right,left, j,L
    real::Delta_E, E_0, E_1

    print *, '--------- From Cluster Subroutine------------'
    ! call rand_generator(L, x, y)

    do y=1, L
        print *, "Row ", y
        do x=1, L
            print *, "Column ", x
            if(x.eq.1) then
                left = L
            else 
                left = x-1
            end if
            if(x.eq.L) then
                right = 1
            else 
                right = x+1
            end if
            if(y.eq.1) then
                top = L
            else 
                top = y-1
            end if
            if(y.eq.L) then
                bottom = 1
            else 
                bottom = y+1
            end if
        
            E_0 = j*gl(x,y)*(gl(left, y) + gl(x, top) + gl(right, y) + gl(x, bottom))
            E_1 = j*(-gl(x,y))*(gl(left, y) + gl(x, top) + gl(right, y) + gl(x, bottom))
        
            Delta_E = E_1 - E_0
        
            if (Delta_E.le.0) then
                call flip_spin(gl, L, x, y)
            else
                call recheck_energi(Delta_E, gl, L, x, y)
            end if

            print *, " "
        end do
            print *, " "
            print *, " "
            print *, " "
    end do

    ! print *, 'x = ', x, 'y = ', y
    ! print *, 'Energi Before =', E_0
    ! print *, 'Energi After =', E_1
    ! print *, 'Energi Delta =', Delta_E
    ! print *, 'N1(Left) = ', gl(left,y), ' N2(Top) = ', gl(x, top)
    ! print *, 'N3(Right) = ', gl(right, y), ' N4(Bottom) = ', gl(x, bottom)
        
end subroutine

subroutine flip_spin(gl, L, x, y)
    real, dimension(L,L)::gl
    integer::x,y,L

    print *, "Flip Spin"
    
    gl(x,y) = - gl(x,y)
    write(*, 2) gl
    2       format(8f6.1)
end subroutine

subroutine recheck_energi(del_E, gl, L, x, y)
    real, dimension(L,L)::gl
    real::r, p
    integer::x,y,L,T = 5, k = 1

    print *, "Recheck Energi"

    call random_seed
    call random_number(r)

    p = exp(-del_E/(k*T))

    if(p.gt.r) then
        call flip_spin(gl, L, x, y)
    else
        print *, "Flip-Spin Canceled"
    end if    

end subroutine

subroutine rand_generator(L, x, y)
    implicit none
    integer:: L, x, y
    real rand0, rand1

    call random_number(rand0)
    call random_number(rand1)

    x = int(1 + rand0 + rand0*(L-2))
    y = int(1 + rand1 + rand1*(L-2))

end subroutine