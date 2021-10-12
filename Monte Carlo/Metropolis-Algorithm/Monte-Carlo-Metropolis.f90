program Montecarlo_Metropolis    
    implicit none
    integer,parameter::L_size = 8, J=1, mcs = 10000, k =1
    integer::iter
    real::f_temp, l_temp, energy_tot, Beta
    real,parameter::T = 5, minT = 0.5, delT = 0.1
    real, dimension(8,8)::n
    real, allocatable, dimension(:)::all_energy

    !Initialize lattice configuration
    call initialize(L_size, n)

    !Temperature Loop
    f_temp = T
    l_temp = minT
    do
        if(f_temp.gt.l_temp) then
            !Metropolis-Algorithm
            print *, "Temperature =", f_temp
            Beta = k*f_temp
            do iter = 1, mcs
                call energy_position(n, L_size, J, energy_tot, Beta)    
            end do
        else
            exit
        end if

        f_temp = f_temp - delT
        all_energy = [all_energy, energy_tot]
        print *, "All Energy = ", all_energy
    end do

end program Montecarlo_Metropolis

subroutine initialize(L, n)        !L = Lattice Size
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

subroutine energy_position(gl, L, j, energy_tot, Beta)
    real, dimension(L,L)::gl
    integer::x,y,top,bottom,right,left, j,L
    real::Delta_E, E_0, E_1, energy, energy_tot, Beta

    energy = 0
    do y=1, L
        do x=1, L
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
        
            ! E_0 = -j*gl(x,y)*(gl(left, y) + gl(x, top) + gl(right, y) + gl(x, bottom))
            ! E_1 = -j*(-gl(x,y))*(gl(left, y) + gl(x, top) + gl(right, y) + gl(x, bottom))
       
            Delta_E = -2*gl(x,y)
            call energy_total(Delta_E, energy)
        
            if (Delta_E.le.0) then
                call flip_spin(gl, L, x, y)
            else
                call recheck_energi(Delta_E, gl, L, x, y, Beta)
            end if
        end do
    end do
    
    energy_tot = energy/2.0

end subroutine

subroutine energy_total(dE, energy)
    real:: dE, energy
    energy = energy + dE
end subroutine

subroutine flip_spin(gl, L, x, y)
    real, dimension(L,L)::gl
    integer::x,y,L

    gl(x,y) = - gl(x,y)
    ! write(*, 2) gl
    ! 2       format(8f6.1)
end subroutine

subroutine recheck_energi(del_E, gl, L, x, y, Beta)
    real, dimension(L,L)::gl
    real::r, p, del_E, Beta  !Beta = k*T
    integer::x,y,L


    call random_seed
    call random_number(r)

    p = exp(-del_E/(Beta))

    if(p.gt.r) then
        call flip_spin(gl, L, x, y)
    end if    

end subroutine

subroutine rand_generator(L, x, y)
    integer:: L, x, y
    real rand0, rand1

    call random_number(rand0)
    call random_number(rand1)

    x = int(1 + rand0 + rand0*(L-2))
    y = int(1 + rand1 + rand1*(L-2))

end subroutine