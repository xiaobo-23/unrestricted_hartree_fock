!  02/22/2023
!  Unrestricted Hartree-Fock solutions to the Hubbard model at finite-temperature in a grand canonical ensemble.
program IHF
        
    ! use compute_energy_and_density  
    implicit none
    ! Lattice parameters including lattice sizes etc.
    integer, parameter :: precision=8
    integer, parameter :: Nx=16, Ny=8
    integer, parameter :: N = Nx * Ny
    integer, parameter ::  iran=20000
    integer, parameter :: number_of_trial_states=6
      
    ! integer :: ix, iy 
    integer :: xplus(N), xminus(N) 
    integer :: yplus(N), yminus(N)
    integer :: i, j, jprime
    
    ! Single-particle Hamiltonian, eigenvalues and eignvectors etc.
    real(kind=precision) :: Hup(N,N), evalup(N), evecup(N,N)
    real(kind=precision) :: Hdn(N,N), evaldn(N), evecdn(N,N) 
    real(kind=precision) :: Hupcopy(N,N), Hdncopy(N,N)
    real(kind=precision), dimension(N) :: nup
    real(kind=precision), dimension(N) :: ndn
    real(kind=precision) :: E, Free_Energy
    real(kind=precision) :: newnup(N), newndn(N)

    ! Electron hopping amplitudes 
    real(kind=precision), parameter :: tx=1.0d0
    real(kind=precision), parameter :: ty=1.0d0
    real(kind=precision), parameter :: tprime=-0.2d0 
    real(kind=precision), parameter :: scale=0.5d0
    
    ! Spin pinning fields which break SU(2) symmetry in the original Hamiltonian
    real(kind=precision), parameter :: pin_amplitude=0.2d0          
   
    ! Tune the chemical potential mu
    real(kind=precision) :: mu
    real(kind=precision) :: mu_min, mu_max
    real(kind=precision) :: mu_delta
    real(kind=precision) :: rho_relaxed
    real(kind=precision), parameter :: tolerance=1.0d-8

    real(kind=precision), parameter :: rho=0.875d0
    real(kind=precision), parameter :: beta=5
    real(precision), parameter :: U=2.8 
    
    integer, parameter :: Nit=3500
    integer :: it 

    ! Speed up convergence using annealing and simple mixing
    real(kind=precision), parameter :: annealing_amp=0.7d0    
    real(kind=precision), parameter :: relax = 0.7d0 

    integer :: INFO
    integer, parameter :: N6 = 6 * N
    real(kind=precision) :: work(N6)

    
    ! Store temporary physical observables in the optimization process
    integer :: index
    integer :: trial_index
    integer :: optimal_index
    real(kind=precision) :: totalF(number_of_trial_states)
    real(kind=precision) :: upElectron(N, number_of_trial_states) 
    real(kind=precision) :: dnElectron(N, number_of_trial_states)
    ! real(kind=precision) :: totalE(Nit, number_of_trial_states)

! **********************************************************************
! **********************************************************************

    open(unit=68, status='replace', file='Data/Benchmark/Pinning/L16W8/output_tprime-0p2_U2p8_beta5.dat')
    open(unit=69, status='replace', file='Data/Benchmark/Pinning/L16W8/spin_density_tprime-0p2_U2p8_beta5.dat')
    open(unit=70, status='replace', file='Data/Benchmark/Pinning/L16W8/free_energy_tprime-0p2_U2p8_beta5.dat')

    do trial_index = 1, number_of_trial_states
        !************************************************************************************************************
        !  Construct the non-interacting Hamiltonian
        !************************************************************************************************************
        call set_up_periodic_nearest_neighbor_list(Nx, Ny, N, xplus, xminus, yplus, yminus)
        if (pin_amplitude > 1E-8) then
            call set_up_open_and_periodic_hamiltonian(Nx, Ny, N, tx, ty, tprime, xplus, xminus, yplus, yminus, Hup, Hdn)
            print *, 'Using spin pinning fields on one edge of the cylinder', pin_amplitude
        else
            call set_up_periodic_hamiltonian(Nx, Ny, N, tx, ty, tprime, xplus, xminus, yplus, yminus, Hup, Hdn)
        end if
        ! print *, 'The non-interacting spin-up matrix'
        ! do i = 1, 2
        !     print *, Hup(i, :)
        ! enddo    

        !************************************************************************************************************
        !  Initialize the distribution of spin-up and spin-down electrons
        !************************************************************************************************************
        print *, 'Initial electron distribution'
        if (trial_index == number_of_trial_states) then
            call initialize_electron_distribution_AFM(Nx, N, rho, nup, ndn)
        else
            call initialize_electron_distribution_stripe(Nx, N, rho, nup, ndn, trial_index)
        end if
        
        !************************************************************************************************************
        !  Set up the interacting part of the Hamiltonian
        !************************************************************************************************************
        ! print *, 'Edge pinning fields'
        call set_up_interaction(Nx, Ny, N, U, pin_amplitude, nup, ndn, Hup, Hdn)
        ! print *, 'The spin-up matrix'
        ! do i = 1, 5
        !     print *, Hup(i, :)
        ! enddo

        ! Initialize the chemical potential
        mu = 0.0d0
        write (68, "()") 
        write (68, "()")
        write (68, "('# of trial state = ', i6)") trial_index
        write (68, "()")
        write (68, "()")

        write (70, "()") 
        write (70, "()")
        write (70, "('# of trial state = ', i6)") trial_index
        write (70, "()")
        write (70, "()")

        !************************************************************************************************************
        !  Solving the Hartree-Fock equations self-consistently
        !************************************************************************************************************
        do it = 1, Nit
            ! print "(' iteration = ', i6)", it
            INFO = 0

            ! Copy the original Hamiltonian matrix and diagonalize it
            call generate_matrix_copy(N, Hup, Hdn, Hupcopy, Hdncopy)
            call DSYEV('V', 'U', N, Hupcopy, N, evalup, WORK, 6*N, INFO)
            call DSYEV('V', 'U', N, Hdncopy, N, evaldn, WORK, 6*N, INFO)
            ! print *, 'The eigenvalues of spin-up matrix'
            ! print *, evalup

            ! print *, 'The eigenvalues of spin-down matrix'
            ! print *, evaldn
            
            ! Obtain the eigenvectors of the original Hamiltonian
            do i=1,N
                do j=1,N
                    evecup(i,j)=Hupcopy(i,j)
                    evecdn(i,j)=Hdncopy(i,j)
                enddo
            enddo
            
            call compute_energy(it, N, beta, mu, U, nup, ndn, evalup, evaldn, E)
            call compute_density(N, beta, mu, evalup, evaldn, evecup, evecdn, newnup, newndn)
            
            mu_min = -10.0d0
            mu_max =  10.0d0
            mu_delta = 0.0d0
            call chemical_potential_tuning(mu_min, mu_max, rho, tolerance, N, beta, mu, evalup, evaldn,& 
            evecup, evecdn, mu_delta)
            mu = mu + mu_delta
            call compute_density(N, beta, mu, evalup, evaldn, evecup, evecdn, newnup, newndn)

            ! Write temporary data e.g. energy per site & density histogram in the output file
            write (68, "('it, E/N = ', i6, f16.8)") it, E/dfloat(N)
            if (mod(50 * it, Nit) == 0) then
                write (68, "()")
                write (68, "(' old/new densities up and dn' )")
                write (68, "(' iteration =  ', i6)") it
                do i = 1, N
                    write(68, "(i6, 2f12.6, 2f12.6)") i, nup(i), newnup(i), ndn(i), newndn(i)
                enddo
            end if
            
            ! Compute the free energy per site
            call compute_free_energy(N, beta, mu, rho, U, evalup, evaldn, nup, ndn, Free_Energy)
            write(70, "('it, F/N=', i6, f16.8)") it, Free_Energy
            ! print '(i6, f12.6)', it, Free_Energy

            ! Use annealing to speed up the convergen in the first half of calculation & use simple mixing for every step except for the last step
            if (it /= Nit) then
                ! Add random noise on each site to speed up the convergence
                if ((mod(it, 50) == 0) .and. (it < 2000)) then
                    call annealing(N, rho, annealing_amp, newnup, newndn)
                endif

                ! Add simple mixing with a coefficient \alpha
                call simple_mixing(N, relax, U, rho_relaxed, Hup, Hdn, nup, ndn, newnup, newndn)
            end if 
        enddo
    enddo 

    do i = 1, N
        write(69, '(i6, 2f12.6)') i, nup(i), ndn(i)
    enddo

    close(68)
    close(69)
    close(70)

contains
    ! Split a real number into integer and floating parts
    subroutine split(input_number, amplifier, tmp_int, tmp_float)
        real(kind=precision), intent(in) :: input_number
        integer, intent(in) :: amplifier
        integer, intent(out) :: tmp_int
        integer, intent(out) :: tmp_float

        tmp_int = int(input_number + 1.0d-8)
        tmp_float = mod(nint(input_number * amplifier), amplifier)
    end subroutine split


    ! Set up the neighboring list
    subroutine set_up_periodic_nearest_neighbor_list(tmp_Nx, tmp_Ny, tmp_Nsites,&
        tmp_xplus, tmp_xminus, tmp_yplus, tmp_yminus)
            integer :: ix, iy
            integer, intent(in) :: tmp_Nx, tmp_Ny, tmp_Nsites
            integer, intent(out) :: tmp_xplus(tmp_Nsites), tmp_xminus(tmp_Nsites)
            integer, intent(out) :: tmp_yplus(tmp_Nsites), tmp_yminus(tmp_Nsites)
    
            do ix = 1, tmp_Nx 
                tmp_xplus(ix) = ix + 1
                tmp_xminus(ix) = ix - 1
            end do
            tmp_xplus(tmp_Nx) = 1 
            tmp_xminus(1) = tmp_Nx 
    
            do iy = 1, tmp_Ny
                tmp_yplus(iy) = iy + 1
                tmp_yminus(iy) = iy - 1
            end do 
            tmp_yplus(tmp_Ny) = 1 
            tmp_yminus(1) = tmp_Ny
        end subroutine set_up_periodic_nearest_neighbor_list


    ! Set up the non-interacting part of the Hamilltonian using periodic boundary conditions in both x and y directions
    subroutine set_up_periodic_hamiltonian(tmp_Nx, tmp_Ny, tmp_Nsites, tmp_tx, tmp_ty, tmp_tprime,&
    tmp_xplus, tmp_xminus, tmp_yplus, tmp_yminus, tmp_Hup, tmp_Hdn)
        integer, intent(in) :: tmp_Nx, tmp_Ny, tmp_Nsites
        real(kind=precision), intent(in) :: tmp_tx, tmp_ty, tmp_tprime
        real(kind=precision), intent(out) :: tmp_Hup(tmp_Nsites, tmp_Nsites), tmp_Hdn(tmp_Nsites, tmp_Nsites)

        integer :: i, j, ix, iy
        integer, intent(in) :: tmp_xplus(tmp_Nsites), tmp_xminus(tmp_Nsites)
        integer, intent(in) :: tmp_yplus(tmp_Nsites), tmp_yminus(tmp_Nsites)

        do i = 1, N
            do j = 1, N
            tmp_Hup(i, j) = 0.0d0
            tmp_Hdn(i, j) = 0.0d0
            end do
            ! tmp_Hup(i, i) = -tmp_mu
            ! tmp_Hdn(i, i) = -tmp_mu
        end do

        do ix = 1, tmp_Nx
            do iy = 1, tmp_Ny
                i =  ix + (iy - 1) * tmp_Nx
                j =  tmp_xplus(ix) + (iy - 1) * tmp_Nx
                tmp_Hup(i,j) = -tmp_tx
                tmp_Hdn(i,j) = -tmp_tx

                j = tmp_xminus(ix) + (iy - 1) * tmp_Nx
                tmp_Hup(i,j) = -tmp_tx
                tmp_Hdn(i,j) = -tmp_tx

                j = ix + (tmp_yplus(iy) - 1) * tmp_Nx
                tmp_Hup(i,j) = -tmp_ty
                tmp_Hdn(i,j) = -tmp_ty

                j = ix + (tmp_yminus(iy) - 1) * tmp_Nx
                tmp_Hup(i,j) = -tmp_ty
                tmp_Hdn(i,j) = -tmp_ty

                jprime = tmp_xplus(ix) + (tmp_yplus(iy) - 1) * tmp_Nx 
                tmp_Hup(i,jprime) = -tmp_tprime
                tmp_Hdn(i,jprime) = -tmp_tprime

                jprime = tmp_xminus(ix) + (tmp_yplus(iy) - 1) * tmp_Nx
                tmp_Hup(i,jprime) = -tmp_tprime
                tmp_Hdn(i,jprime) = -tmp_tprime

                jprime = tmp_xminus(ix) + (tmp_yminus(iy) - 1) * tmp_Nx
                tmp_Hup(i,jprime) = -tmp_tprime
                tmp_Hdn(i,jprime) = -tmp_tprime

                jprime = tmp_xplus(ix) + (tmp_yminus(iy) - 1) * tmp_Nx
                tmp_Hup(i,jprime) = -tmp_tprime
                tmp_Hdn(i,jprime) = -tmp_tprime
            end do
        end do
    end subroutine set_up_periodic_hamiltonian


    ! Set up the non-interacting part of the Hamilltonian using open boundary condition in the x direction and periodic boundary condition in y direction
    subroutine set_up_open_and_periodic_hamiltonian(tmp_Nx, tmp_Ny, tmp_Nsites, tmp_tx, tmp_ty, tmp_tprime,&
        tmp_xplus, tmp_xminus, tmp_yplus, tmp_yminus, tmp_Hup, tmp_Hdn)
            integer, intent(in) :: tmp_Nx, tmp_Ny, tmp_Nsites
            real(kind=precision), intent(in) :: tmp_tx, tmp_ty, tmp_tprime
            real(kind=precision), intent(out) :: tmp_Hup(tmp_Nsites, tmp_Nsites), tmp_Hdn(tmp_Nsites, tmp_Nsites)
    
            integer :: i, j, ix, iy
            integer, intent(in) :: tmp_xplus(tmp_Nsites), tmp_xminus(tmp_Nsites)
            integer, intent(in) :: tmp_yplus(tmp_Nsites), tmp_yminus(tmp_Nsites)
    
            do i = 1, N
                do j = 1, N
                tmp_Hup(i, j) = 0.0d0
                tmp_Hdn(i, j) = 0.0d0
                end do
            end do
    
            do ix = 1, tmp_Nx
                do iy = 1, tmp_Ny
                    if (ix == 1) then
                        i =  ix + (iy - 1) * tmp_Nx
                        j =  tmp_xplus(ix) + (iy - 1) * tmp_Nx
                        tmp_Hup(i,j) = -tmp_tx
                        tmp_Hdn(i,j) = -tmp_tx
        
                        j = ix + (tmp_yplus(iy) - 1) * tmp_Nx
                        tmp_Hup(i,j) = -tmp_ty
                        tmp_Hdn(i,j) = -tmp_ty
        
                        j = ix + (tmp_yminus(iy) - 1) * tmp_Nx
                        tmp_Hup(i,j) = -tmp_ty
                        tmp_Hdn(i,j) = -tmp_ty
        
                        jprime = tmp_xplus(ix) + (tmp_yplus(iy) - 1) * tmp_Nx 
                        tmp_Hup(i,jprime) = -tmp_tprime
                        tmp_Hdn(i,jprime) = -tmp_tprime
        
                        jprime = tmp_xplus(ix) + (tmp_yminus(iy) - 1) * tmp_Nx
                        tmp_Hup(i,jprime) = -tmp_tprime
                        tmp_Hdn(i,jprime) = -tmp_tprime
                    else if (ix == tmp_Nx) then
                        i =  ix + (iy - 1) * tmp_Nx
                        j = tmp_xminus(ix) + (iy - 1) * tmp_Nx
                        tmp_Hup(i,j) = -tmp_tx
                        tmp_Hdn(i,j) = -tmp_tx
        
                        j = ix + (tmp_yplus(iy) - 1) * tmp_Nx
                        tmp_Hup(i,j) = -tmp_ty
                        tmp_Hdn(i,j) = -tmp_ty
        
                        j = ix + (tmp_yminus(iy) - 1) * tmp_Nx
                        tmp_Hup(i,j) = -tmp_ty
                        tmp_Hdn(i,j) = -tmp_ty
        
                        jprime = tmp_xminus(ix) + (tmp_yplus(iy) - 1) * tmp_Nx
                        tmp_Hup(i,jprime) = -tmp_tprime
                        tmp_Hdn(i,jprime) = -tmp_tprime
        
                        jprime = tmp_xminus(ix) + (tmp_yminus(iy) - 1) * tmp_Nx
                        tmp_Hup(i,jprime) = -tmp_tprime
                        tmp_Hdn(i,jprime) = -tmp_tprime
                    else
                        i =  ix + (iy - 1) * tmp_Nx
                        j =  tmp_xplus(ix) + (iy - 1) * tmp_Nx
                        tmp_Hup(i,j) = -tmp_tx
                        tmp_Hdn(i,j) = -tmp_tx
        
                        j = tmp_xminus(ix) + (iy - 1) * tmp_Nx
                        tmp_Hup(i,j) = -tmp_tx
                        tmp_Hdn(i,j) = -tmp_tx
        
                        j = ix + (tmp_yplus(iy) - 1) * tmp_Nx
                        tmp_Hup(i,j) = -tmp_ty
                        tmp_Hdn(i,j) = -tmp_ty
        
                        j = ix + (tmp_yminus(iy) - 1) * tmp_Nx
                        tmp_Hup(i,j) = -tmp_ty
                        tmp_Hdn(i,j) = -tmp_ty
        
                        jprime = tmp_xplus(ix) + (tmp_yplus(iy) - 1) * tmp_Nx 
                        tmp_Hup(i,jprime) = -tmp_tprime
                        tmp_Hdn(i,jprime) = -tmp_tprime
        
                        jprime = tmp_xminus(ix) + (tmp_yplus(iy) - 1) * tmp_Nx
                        tmp_Hup(i,jprime) = -tmp_tprime
                        tmp_Hdn(i,jprime) = -tmp_tprime
        
                        jprime = tmp_xminus(ix) + (tmp_yminus(iy) - 1) * tmp_Nx
                        tmp_Hup(i,jprime) = -tmp_tprime
                        tmp_Hdn(i,jprime) = -tmp_tprime
        
                        jprime = tmp_xplus(ix) + (tmp_yminus(iy) - 1) * tmp_Nx
                        tmp_Hup(i,jprime) = -tmp_tprime
                        tmp_Hdn(i,jprime) = -tmp_tprime
                    end if
                end do
            end do
        end subroutine set_up_open_and_periodic_hamiltonian

    
    ! Initialize the distribution of spin-up and spin-down electron densities
    subroutine initialize_electron_distribution_AFM(tmp_Nx, tmp_Nsites, density, tmp_nup, tmp_ndn)
        integer :: index, i, j
        integer, intent(in) :: tmp_Nx, tmp_Nsites
        real(kind=precision), intent(in) :: density
        real(kind=precision), intent(out) :: tmp_nup(tmp_Nsites), tmp_ndn(tmp_Nsites)

        do index = 1, tmp_Nsites
            i = mod(index, tmp_Nx)
            if (i == 0) then
                i = tmp_Nx
            end if
            j = int((index - 1) / Nx + 1)

            if (mod(i + j, 2) == 0) then
                tmp_nup(index) = 1.0d0 * density
                tmp_ndn(index) = 0.0d0 * density
            else
                tmp_nup(index) = 0.0d0 * density
                tmp_ndn(index) = 1.0d0 * density
            end if 
            print '(i6, 2f12.6)', index, tmp_nup(index), tmp_ndn(index) 
        end do
    end subroutine initialize_electron_distribution_AFM


    
    
    ! Initialize the distribution of spin-up and spin-down electron densities using a stripe pattern
    subroutine initialize_electron_distribution_stripe(tmp_Nx, tmp_Nsites, density, tmp_nup, tmp_ndn, tmp_period)
        integer :: index, i, j, tmp_shift
        integer, intent(in) :: tmp_Nx, tmp_Nsites, tmp_period
        real(kind=precision), intent(in) :: density
        real(kind=precision), intent(out) :: tmp_nup(tmp_Nsites), tmp_ndn(tmp_Nsites)

        do index = 1, tmp_Nsites
            i = mod(index, tmp_Nx)
            if (i == 0) then
                i = tmp_Nx
            end if
            j = int((index - 1) / Nx + 1)
            tmp_shift = int((index - 1) / tmp_period)

            if (mod(i + j + tmp_shift, 2) == 0) then
                tmp_nup(index) = 1.0d0 * density
                tmp_ndn(index) = 0.0d0 * density
            else
                tmp_nup(index) = 0.0d0 * density
                tmp_ndn(index) = 1.0d0 * density
            end if 
            print '(i6, 3f12.6)', index, tmp_nup(index), tmp_ndn(index), 0.5 * (tmp_nup(index) - tmp_ndn(index))
        end do
    end subroutine initialize_electron_distribution_stripe


    ! Set up the interacting part of the Hamiltonian and pinning fields etc.
    subroutine set_up_interaction(tmp_Nx, tmp_Ny, tmp_Nsites, tmp_U, pinning_strength,&
    tmp_nup, tmp_ndn, tmp_Hup, tmp_Hdn)
        integer :: ind, tmp_ind
        integer, intent(in) :: tmp_Nx, tmp_Ny, tmp_Nsites
        real(kind=precision), intent(in) :: tmp_U, pinning_strength
        real(kind=precision), intent(in) :: tmp_nup(tmp_Nsites), tmp_ndn(tmp_Nsites)
        real(kind=precision), intent(out) :: tmp_Hup(tmp_Nsites, tmp_Nsites), tmp_Hdn(tmp_Nsites, tmp_Nsites)

        do ind = 1, tmp_Nsites
            tmp_Hup(ind, ind) = tmp_Hup(ind, ind) + tmp_U * tmp_ndn(ind)
            tmp_Hdn(ind, ind) = tmp_Hdn(ind, ind) + tmp_U * tmp_nup(ind)
            ! print '(i6, 2f12.6)', ind, tmp_Hup(ind, ind), tmp_Hdn(ind, ind)
        end do

        if (pinning_strength > 1e-8) then
            do ind = 1, tmp_Ny
                tmp_ind = (ind - 1) * tmp_Nx + 1
                tmp_Hup(tmp_ind, tmp_ind) = tmp_Hup(tmp_ind, tmp_ind) + (-1)**(ind) * 0.5d0 * pinning_strength
                tmp_Hdn(tmp_ind, tmp_ind) = tmp_Hdn(tmp_ind, tmp_ind) + (-1)**(ind + 1) * 0.5d0 * pinning_strength 
                ! print '(i6, f12.6)', ind, tmp_Hup(tmp_ind, tmp_ind) - tmp_Hdn(tmp_ind, tmp_ind)
                ! print '(i6, 2f12.6)', ind, tmp_Hup(tmp_ind, tmp_ind), tmp_Hdn(tmp_ind, tmp_ind)
            end do
        end if          
    end subroutine set_up_interaction

    
    ! Generate a copy of the matrix
    subroutine generate_matrix_copy(tmp_Nsites, input_Hup, input_Hdn, output_Hup, output_Hdn)
        integer, intent(in) :: tmp_Nsites
        real(kind=precision), intent(in) :: input_Hup(tmp_Nsites, tmp_Nsites), input_Hdn(tmp_Nsites, tmp_Nsites)
        real(kind=precision), intent(out) :: output_Hup(tmp_Nsites, tmp_Nsites), output_Hdn(tmp_Nsites, tmp_Nsites)

        integer :: i, j

        do i = 1, tmp_Nsites
            do j = 1, tmp_Nsites
                output_Hup(i, j) = input_Hup(i, j)
                output_Hdn(i, j) = input_Hdn(i, j)
                ! print '(2f12.6)', output_Hup(i,j) - input_Hup(i,j), output_Hdn(i, j) - output_Hdn(i, j)
            enddo
        enddo
    end subroutine generate_matrix_copy

    
    ! Compute the total energy using the Fermi-Dirac distribution
    subroutine compute_energy(tmp_iteration, tmp_Nsites, tmp_beta, tmp_mu, tmp_U,&
    tmp_nup, tmp_ndn, tmp_evalup, tmp_evaldn, tmp_energy)
        integer, intent(in) :: tmp_iteration, tmp_Nsites
        real(kind=precision), intent(in) :: tmp_beta, tmp_mu, tmp_U
        real(kind=precision), intent(in) :: tmp_nup(tmp_Nsites), tmp_ndn(tmp_Nsites)
        real(kind=precision), intent(in) :: tmp_evalup(tmp_Nsites), tmp_evaldn(tmp_Nsites)
        real(kind=precision), intent(out) :: tmp_energy
        
        integer :: ind
        real(kind=precision) :: tmp_fermiup, tmp_fermidn

        tmp_energy = 0.0d0
        do ind = 1, tmp_Nsites
            tmp_fermiup = 1.0d0 / (exp(tmp_beta * (tmp_evalup(ind) - tmp_mu)) + 1.0d0)
            tmp_fermidn = 1.0d0 / (exp(tmp_beta * (tmp_evaldn(ind) - tmp_mu)) + 1.0d0)
            tmp_energy = tmp_energy + tmp_fermiup * tmp_evalup(ind) + tmp_fermidn * tmp_evaldn(ind)
            tmp_energy = tmp_energy - tmp_U * tmp_nup(ind) * tmp_ndn(ind)
        enddo
        ! print "('iteration, E/N = ', i6, f16.8)", tmp_iteration, tmp_energy / dfloat(tmp_Nsites)
        ! write (68, "('it, E/N =   ', i8, f16.8)") it, E/dfloat(N)
        ! totalE(it,index) = E/dfloat(N)
    end subroutine compute_energy


    ! Compute the electron density based on the eigenvalues and eigenvectors
    subroutine compute_density(tmp_Nsites, tmp_beta, tmp_mu,& 
    tmp_evalup, tmp_evaldn, tmp_evecup, tmp_evecdn, tmp_newnup, tmp_newndn)
        integer, intent(in) :: tmp_Nsites
        real(kind=precision), intent(in) :: tmp_beta, tmp_mu
        real(kind=precision), intent(in) :: tmp_evalup(tmp_Nsites), tmp_evaldn(tmp_Nsites)
        real(kind=precision), intent(in) :: tmp_evecup(tmp_Nsites, tmp_Nsites), tmp_evecdn(tmp_Nsites, tmp_Nsites)
        real(kind=precision), intent(out) :: tmp_newnup(tmp_Nsites), tmp_newndn(tmp_Nsites)

        integer :: ind1, ind2
        real(kind=precision) :: tmp_fermiup, tmp_fermidn

        do ind1 = 1, tmp_Nsites
            tmp_newnup(ind1) = 0.0d0
            tmp_newndn(ind1) = 0.0d0

            do ind2 = 1, tmp_Nsites
                tmp_fermiup = 1.0d0 / (exp(tmp_beta * (tmp_evalup(ind2) - tmp_mu)) + 1.0d0)
                tmp_newnup(ind1) = tmp_newnup(ind1) + tmp_fermiup * tmp_evecup(ind1, ind2) * tmp_evecup(ind1, ind2)

                tmp_fermidn = 1.0d0 / (exp(tmp_beta * (tmp_evaldn(ind2) - tmp_mu)) + 1.0d0)
                tmp_newndn(ind1) = tmp_newndn(ind1) + tmp_fermidn * tmp_evecdn(ind1, ind2) * tmp_evecdn(ind1, ind2)
            enddo
        enddo
    end subroutine compute_density

    
    ! Compute the free energy per site
    subroutine compute_free_energy(tmp_Nsites, tmp_beta, tmp_mu, tmp_rho, tmp_U, &
        tmp_evalup, tmp_evaldn, tmp_nup, tmp_ndn, tmp_free_energy)
        integer, intent(in) :: tmp_Nsites
        real(kind=precision), intent(in) :: tmp_beta
        real(kind=precision), intent(in) :: tmp_mu
        real(kind=precision), intent(in) :: tmp_rho
        real(kind=precision), intent(in) :: tmp_U
        real(kind=precision), intent(in) :: tmp_evalup(tmp_Nsites), tmp_evaldn(tmp_Nsites)
        real(kind=precision), intent(in) :: tmp_nup(tmp_Nsites), tmp_ndn(tmp_Nsites)
        real(kind=precision), intent(out) :: tmp_free_energy

        tmp_free_energy = 0.0d0
        do index = 1, tmp_Nsites
                tmp_free_energy = tmp_free_energy + log(1.0d0 + exp(-tmp_beta * (tmp_evalup(index) - tmp_mu)))
                tmp_free_energy = tmp_free_energy + log(1.0d0 + exp(-tmp_beta * (tmp_evaldn(index) - tmp_mu)))
                tmp_free_energy = tmp_free_energy + tmp_beta * tmp_U * tmp_nup(index) * tmp_ndn(index)
        enddo
        tmp_free_energy = -1.0d0 / tmp_beta * tmp_free_energy
        tmp_free_energy = tmp_free_energy / dfloat(tmp_Nsites) + tmp_mu * tmp_rho
    end subroutine compute_free_energy

    
    ! Tune the chemical potential to reach the target density
    subroutine chemical_potential_tuning(lower_bound, upper_bound, target_density, tmp_tolerance, tmp_Nsites, tmp_beta, tmp_mu,& 
    tmp_evalup, tmp_evaldn, tmp_evecup, tmp_evecdn, tmp_delta_mu)
        integer, intent(in) :: tmp_Nsites
        integer :: ind1, ind2
        real(kind=precision), intent(in) :: tmp_tolerance, target_density
        real(kind=precision), intent(in) :: tmp_beta, tmp_mu
        real(kind=precision), intent(in) :: tmp_evalup(tmp_Nsites), tmp_evaldn(tmp_Nsites)
        real(kind=precision), intent(in) :: tmp_evecup(tmp_Nsites, tmp_Nsites), tmp_evecdn(tmp_Nsites, tmp_Nsites)
        real(kind=precision) :: lower_bound, upper_bound
        real(kind=precision), intent(out) :: tmp_delta_mu
        real(kind=precision) :: density_diff, tmp_density
        real(kind=precision) :: tmp_nup(tmp_Nsites), tmp_ndn(tmp_Nsites)
        real(kind=precision) :: tmp_fermiup, tmp_fermidn
        integer :: ind

        density_diff = 1.0d6
        do while (density_diff > tmp_tolerance)
            tmp_delta_mu = 0.5d0 * (lower_bound + upper_bound)
            tmp_density = 0.d0

            ! call compute_density(tmp_Nsites, tmp_beta, tmp_mu + tmp_delta_mu, tmp_evalup, tmp_evaldn, tmp_evecup, tmp_evecdn, tmp_nup, tmp_ndn)
            do ind1 = 1, tmp_Nsites
                tmp_nup(ind1) = 0.0d0
                tmp_ndn(ind1) = 0.0d0

                do ind2 = 1, tmp_Nsites
                    ! print '(i6)', ind2
                    tmp_fermiup = 1.0d0 / (exp(tmp_beta * (tmp_evalup(ind2) - (tmp_mu + tmp_delta_mu))) + 1.0d0)
                    tmp_nup(ind1) = tmp_nup(ind1) + tmp_fermiup * tmp_evecup(ind1, ind2) * tmp_evecup(ind1, ind2)

                    tmp_fermidn = 1.0d0 / (exp(tmp_beta * (tmp_evaldn(ind2) - (tmp_mu + tmp_delta_mu))) + 1.0d0)
                    tmp_ndn(ind1) = tmp_ndn(ind1) + tmp_fermidn * tmp_evecdn(ind1, ind2) * tmp_evecdn(ind1, ind2)
                enddo
            enddo

            do ind = 1, tmp_Nsites
                tmp_density = tmp_density + tmp_nup(ind) + tmp_ndn(ind)
            enddo

            tmp_density = tmp_density / dfloat(tmp_Nsites)
            density_diff = abs(tmp_density - target_density)
            ! print '(2f12.6)', tmp_delta_mu, tmp_density
            
            if (tmp_density > target_density) then
                upper_bound = tmp_delta_mu
            else if (tmp_density < target_density) then
                lower_bound = tmp_delta_mu
            endif
        enddo
    end subroutine chemical_potential_tuning

    
    ! Add random noise on each site to speed up the convergence
    subroutine annealing(tmp_Nsites, target_density, annealing_amplitude, tmp_newnup, tmp_newndn)
        integer, intent(in) :: tmp_Nsites
        real(kind=precision), intent(in) :: target_density
        real(kind=precision), intent(in) :: annealing_amplitude
        real(kind=precision), intent(out) :: tmp_newnup(tmp_Nsites), tmp_newndn(tmp_Nsites)
        real(kind=precision) :: tmp_total_electrons, tmp_shift
        real(kind=precision) :: tmp1, tmp2   ! Random numbers added to the densities on each site

        tmp_total_electrons = 0.0d0
        tmp_shift = 0.0d0

        do index = 1, tmp_Nsites
                call random_number(tmp1)
                call random_number(tmp2)
                tmp1 = 2 * tmp1 - 1.0d0
                tmp2 = 2 * tmp2 - 1.0d0
                tmp_newnup(index) = tmp_newnup(index) + annealing_amplitude * tmp1
                tmp_newndn(index) = tmp_newndn(index) + annealing_amplitude * tmp2

                ! Rescale spin-resolved electron densities on each site 
                if (tmp_newnup(index) > 1.0d0) then
                    tmp_newnup(index) = 1.0d0
                else if (tmp_newnup(index) < 0.0d0) then
                    tmp_newnup(index) = 0.0d0
                endif 

                if (tmp_newndn(index) > 1.0d0) then
                    tmp_newndn(index) = 1.0d0
                else if (tmp_newndn(index) < 0.0d0) then
                    tmp_newndn(index) = 0.0d0
                endif

                tmp_total_electrons = tmp_total_electrons + tmp_newnup(index) + tmp_newndn(index)
        enddo

        tmp_shift = (tmp_total_electrons - tmp_Nsites * target_density) / tmp_Nsites
        do index = 1, tmp_Nsites
                tmp_newnup(index) = tmp_newnup(index) - tmp_shift
                tmp_newndn(index) = tmp_newndn(index) - tmp_shift
        enddo
    end subroutine annealing

    
    ! Simple mixing of the new and old densities
    subroutine simple_mixing(tmp_Nsites, tmp_relax, tmp_U, tmp_rho, tmp_Hup, tmp_Hdn, &
        tmp_nup, tmp_ndn, tmp_newnup, tmp_newndn)
        integer, intent(in) :: tmp_Nsites
        real(kind=precision), intent(in) :: tmp_relax, tmp_U
        real(kind=precision), intent(out) :: tmp_rho
        real(kind=precision), intent(out) :: tmp_nup(tmp_Nsites), tmp_ndn(tmp_Nsites)
        real(kind=precision), intent(out) :: tmp_newnup(tmp_Nsites), tmp_newndn(tmp_Nsites)
        real(kind=precision), intent(out) :: tmp_Hup(tmp_Nsites, tmp_Nsites), tmp_Hdn(tmp_Nsites, tmp_Nsites)

        tmp_rho = 0.0d0
        do index = 1, tmp_Nsites
                tmp_newnup(index) = tmp_relax * tmp_newnup(index) + (1.0d0 - tmp_relax) * tmp_nup(index)
                tmp_newndn(index) = tmp_relax * tmp_newndn(index) + (1.0d0 - tmp_relax) * tmp_ndn(index)
                tmp_Hup(index, index) = tmp_Hup(index, index) + tmp_U * (tmp_newndn(index) - tmp_ndn(index))
                tmp_Hdn(index, index) = tmp_Hdn(index, index) + tmp_U * (tmp_newnup(index) - tmp_nup(index))
                tmp_nup(index) = tmp_newnup(index)
                tmp_ndn(index) = tmp_newndn(index)
                tmp_rho = tmp_rho + tmp_nup(index) + tmp_ndn(index)
        enddo
        tmp_rho = tmp_rho / dble(tmp_Nsites)
    end subroutine simple_mixing



end program IHF