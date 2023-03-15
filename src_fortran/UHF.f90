!  02/22/2023
!  Unrestricted Hartree-Fock solutions to the Hubbard model at finite-temperature in a grand canonical ensemble.
program IHF
        
      implicit none
      ! use iso_fortran_env

      integer, parameter :: precision=8
      integer, parameter :: Nx=4,  Ny=4
      integer, parameter :: N = Nx * Ny
      integer, parameter :: iran=20000

      
      integer :: ix, iy 
      integer :: xplus(N), xminus(N) 
      integer :: yplus(N), yminus(N)
      integer :: i, j, jprime
      integer :: tmp_i, tmp_j, site
      integer :: tmp
      ! integer :: tmp_sign, tmp_wave, tmp_loc
      ! integer :: tmp

      real(kind=precision) :: Hup(N,N), evalup(N), evecup(N,N)
      real(kind=precision) :: Hupcopy(N,N), Huptmp(N,N)
      real(kind=precision) :: Hdn(N,N), evaldn(N), evecdn(N,N) 
      real(kind=precision) :: Hdncopy(N,N), Hdntmp(N,N)
      real(kind=precision), dimension(N) :: nup
      real(kind=precision), dimension(N) :: ndn
      real(kind=precision) :: E, Free_Energy, fermiup, fermidn 
      real(kind=precision) :: tmp_energy
      real(kind=precision) :: newnup(N), newndn(N)
      real(kind=precision) :: tmpnup(N), tmpndn(N)

      ! real*8 tempup(N,N),checkidup(N,N)
      ! real*8 tempdn(N,N),checkiddn(N,N)

!******************************************************************
!     Add next-nearest neighbor hopping into the code
!******************************************************************
      real(kind=precision), parameter :: tx=1.0d0
      real(kind=precision), parameter :: ty=1.0d0
      real(kind=precision), parameter :: scale=0.5d0
      real(kind=precision) :: mu
      real(kind=precision) :: ran2

!*****************************************************************
!     Add pinning field into the code
!     Add electron density
!*****************************************************************
      real(kind=precision), parameter :: pin_amplitude=0.2d0
      real(kind=precision), parameter :: annealing_amp=0.7d0
      real(kind=precision), parameter :: rho=0.875d0
      real(kind=precision) :: tmp_total_electrons
      real(kind=precision) :: tmp_shift
      real(kind=precision) :: tmp1, tmp2
      logical if_pinning_field
!*****************************************************************
!*****************************************************************


!*****************************************************************
!    Add declarations for the chemical tuning part
!*****************************************************************
      real(kind=precision) :: mu_min, mu_max
      real(kind=precision) :: mu_delta
      real(kind=precision) :: rho_delta
      real(kind=precision) :: rho_tmp
      real(kind=precision) :: rho_relaxed
      real(kind=precision), parameter :: tolerance=1.0D-8
!*****************************************************************
!*****************************************************************


!***********************************************************************
!   Add temperature bounds and steps
!***********************************************************************
      real(precision), parameter :: beta_min=3.0d0
      real(precision), parameter :: beta_max=3.0d0
      real(precision), parameter :: beta_step=1.0d0
      integer :: beta_int, beta_float, beta_prev_int, beta_prev_float
      real(precision) :: beta, beta_prev
! ***********************************************************************
! ***********************************************************************


! ***********************************************************************
!      Add effective U bounds and steps
! ***********************************************************************
      real(precision), parameter :: Umin=3.6d0
      real(precision), parameter :: Umax=3.6d0
      real(precision), parameter :: Ustep=0.1d0
      integer :: Uint, Ufloat
      real(precision) :: U
! ***********************************************************************
! ***********************************************************************


!*****************************************************************
!     Add next nearest-neighbor hopping t' bounds and steps
!*****************************************************************
      real(precision), parameter :: tprime_min=-0.2d0
      real(precision), parameter :: tprime_max=-0.2d0
      real(precision), parameter :: tprime_step=0.05d0
      integer :: tprime_int, tprime_float
      real(precision) :: tprime
!*****************************************************************
!*****************************************************************

!***********************************************************************
!     Initialize different domain sizes as the starting point of 
!     Hartree-Fock solution
!***********************************************************************
      integer :: phase_shift
      integer :: domain_size
      integer :: index
      integer :: tmpInd
      integer :: optimal_index
      integer, parameter :: number_of_trial_states=1

!***********************************************************************
!***********************************************************************
      integer :: it 
      integer, parameter :: Nit = 500
      real(kind=precision), parameter :: relax = 0.7d0 
      integer :: INFO
      integer, parameter :: N6 = 6 * N
      real(kind=precision) :: work(N6)
      character(len = :), allocatable :: formatting_string
      character(len = 100) outname
      character(len = 100) input_name
      character(len = 100) spin
      character(len = 100) free_energy_string
      character(len = 100) grand_free_energy_string
      character(len = 100) e_hist
      character(len = 100) e_density
      character(len = :), allocatable :: tprime_string
      character(len = :), allocatable :: tmp_string

! **********************************************************************
! **********************************************************************
      real(kind=precision) :: totalE(Nit, number_of_trial_states)
      real(kind=precision) :: totalF(Nit, number_of_trial_states)
      real(kind=precision) :: totalGF(Nit, number_of_trial_states)
      real(kind=precision) :: upElectron(N, number_of_trial_states) 
      real(kind=precision) :: dnElectron(N, number_of_trial_states)
! **********************************************************************
! **********************************************************************

!! Generate output files
      do tprime = tprime_min, tprime_max, tprime_step
          call split(tprime, 100, tprime_int, tprime_float)
          
          do U = Umin, Umax, Ustep
              call split(U, 10, Uint, Ufloat)
              ! print *, Uint, Ufloat
          
              do beta = beta_min, beta_max, beta_step
                  call split(beta, 10, beta_int, beta_float)
                  ! print *, beta_int, beta_float
              

            !   if (beta_int>=10 .and. abs(tprime_float)>=10) then
            !       formatting_string = "(a, a, i1, a, i2, a, i1, a, i1, a, i2, a, i1, a)"
            !   else if (beta_int>=10 .and. abs(tprime_float)<10) then
            !     formatting_string = "(a, a, i1, a, i1, a, i1, a, i1, a, i2, a, i1, a)"
            !   else if (beta_int<10 .and. abs(tprime_float)>=10) then
            !     formatting_string = "(a, a, i1, a, i2, a, i1, a, i1, a, i1, a, i1, a)"
            !   else
            !     formatting_string = "(a, a, i1, a, i1, a, i1, a, i1, a, i1, a, i1, a)"
            !   endif

            !   if (tprime >= 0) then
            !     tprime_string = 'tprime'
            !   else
            !     tprime_string = 'tprime-'
            !   end if

            !   if (abs(tprime_float) >= 10) then
            !     tmp_string = 'p'
            !   else
            !     tmp_string = 'p0'
            !   end if

            ! write(outname,  formatting_string) 'Data/p02_test/tprime-0p2/rho0p875/OBC/L16W4/Results_', &
            !       tprime_string, tprime_int, tmp_string, abs(tprime_float), '_U', Uint, 'p', Ufloat, '_beta', beta_int, 'p', beta_float, '.dat'

            ! write(spin,  formatting_string) 'Data/p02_test/tprime-0p2/rho0p875/OBC/L16W4/Spin_', &
            !       tprime_string, tprime_int, tmp_string, abs(tprime_float), '_U', Uint, 'p', Ufloat, '_beta', beta_int, 'p', beta_float, '.dat' 

            ! write(e_hist,  formatting_string) 'Data/p02_test/tprime-0p2/rho0p875/OBC/L16W4/Energy_', &
            !       tprime_string, tprime_int, tmp_string, abs(tprime_float), '_U', Uint, 'p', Ufloat, '_beta', beta_int, 'p', beta_float, '.dat' 

            ! write(e_density, formatting_string) 'Data/p02_test/tprime-0p2/rho0p875/OBC/L16W4/Charge_', &
            !       tprime_string, tprime_int, tmp_string, abs(tprime_float), '_U', Uint, 'p', Ufloat, '_beta', beta_int, 'p', beta_float, '.dat'

            ! write(free_energy_string, formatting_string) 'Data/p02_test/tprime-0p2/rho0p875/OBC/L16W4/Free_Energy_', &
            !       tprime_string, tprime_int, tmp_string, abs(tprime_float), '_U', Uint, 'p', Ufloat, '_beta', beta_int, 'p', beta_float, '.dat'

            ! write(grand_free_energy_string, formatting_string) 'Data/p02_test/tprime-0p2/rho0p875/OBC/L16W4/Grand_Free_Energy_', &
            !       tprime_string, tprime_int, tmp_string, abs(tprime_float), '_U', Uint, 'p', Ufloat, '_beta', beta_int, 'p', beta_float, '.dat'

            ! open(unit=68, status='replace', file=outname)
            ! open(unit=69, status='replace', file=spin)
            ! open(unit=70, status='replace', file=e_hist)
            ! open(unit=71, status='replace', file=e_density)
            ! open(unit=80, status='replace', file=free_energy_string)
            ! open(unit=81, status='replace', file=grand_free_energy_string)

! !*******************************************************************
! !     No input file needed in this verion.
! !*******************************************************************
!       if (pin_amplitude > 1.0d-8) then
!         if_pinning_field = .true.
!       else
!         if_pinning_field = .false.
!       endif
!       write (*, *) if_pinning_field, pin_amplitude
!       mu=0.d0

!       write (68, *) ' ********************************************* '
!       write (68, *) ' Inhomogeneous Hartree-Fock'
!       write (68, *) ' Version 7 '
!       write (68, *) ' ********************************************* '
!       write (68, *) ' Lattice size: Nx,  Ny '
!       write (68, "(2i12)")  Nx, Ny
!       write (68, *) ' Nearest Hopping: tx, ty '
!       write (68, "(2f12.6)") tx, ty
!       write (68, *) ' Next-Nearest Hopping tprime'
!       write (68, "(f12.6)") tprime
!       write (68, *) ' Effective Coulumb Interaction '
!       write (68, "(f12.6)") U
!       write (68, *) ' Original Chemical Potential '
!       write (68, "(f12.6)") mu
!       write (68, *) ' Inverse Temperature'
!       write (68, "(f12.6)") beta
!       write (68, *) ' Pinning Field Strength '
!       write (68, "(f12.6)") pin_amplitude
!       write (68, *) ' Random Initialization Scale '
!       write (68, "(f12.6)") scale
!       write (68, *) ' Random Number Seed'
!       write (68, "(i12)") iran
!       write (68, *) ' Number of Iterations '
!       write (68, "(i12)") Nit
!       write (68, *) ' Relaxation Rate '
!       write (68, "(f12.6)") relax
!       write (68, *) ' If Pinning Is Turned On'
!       write (68, *)   if_pinning_field
!       write (68, *) ' Output File Name '
!       write (68, "(a60)") outname

!       do index = 1, number_of_trial_states
!           write(*, *) index



      !************************************************************************************************************
      !  Construct the non-interacting Hamiltonian
      !************************************************************************************************************
      call set_up_periodic_nearest_neighbor_list(Nx, Ny, N, xplus, xminus, yplus, yminus)
      call set_up_periodic_hamiltonian(Nx, Ny, N, tx, ty, tprime, xplus, xminus, yplus, yminus, Hup, Hdn)
      ! print *, Hup(:, 1)
      ! print *, ""
      ! print *, Hup(:, 2)
      ! print *, ""
      ! print *, Hdn(:, 3)

      !************************************************************************************************************
      !  Initialize the distribution of spin-up and spin-down electrons
      !************************************************************************************************************
      print *, 'Initial electron distribution'
      call initialize_electron_distribution_AFM(Nx, N, rho, nup, ndn)

      !************************************************************************************************************
      !  Set up the interacting part of the Hamiltonian
      !************************************************************************************************************
      print *, 'Edge pinning fields'
      call set_up_interaction(Nx, Ny, N, U, pin_amplitude, nup, ndn, Hup, Hdn)


      !************************************************************************************************************
      !  Perform the self-consistent iterations
      !************************************************************************************************************
      do it = 1, 100
          print "(' iteration = ', i6)", it
          ! write (6, "(' iteration= ', i6)") it

          INFO = 0

          call generate_matrix_copy(N, Hup, Hdn, Hupcopy, Hdncopy)
          call DSYEV('V', 'U', N, Hupcopy, N, evalup, WORK, 6*N, INFO)
          call DSYEV('V', 'U', N, Hdncopy, N, evaldn, WORK, 6*N, INFO)
          
          do i=1,N
              do j=1,N
                  evecup(i,j)=Hupcopy(i,j)
                  evecdn(i,j)=Hdncopy(i,j)
              enddo
          enddo

          call compute_energy(it, N, beta, mu, U, nup, ndn, evalup, evaldn, E)
          call compute_density(N, beta, mu, evalup, evaldn, evecup, evecdn, newnup, newndn)

          write (68, "('it, E/N = ', i6, f16.8)") it, E
          totalE(it, index) = E
          if (mod(20 * it, Nit) == 0) then
              write (68, "()")
              write (68, "(' old/new densities up and dn' )")
              write (68, "(' iteration =  ', i6)") it
              do i = 1, N
                  write(*, "(i6, 2f12.6, 2f12.6)") i, nup(i), newnup(i), ndn(i), newndn(i)
              enddo
          end if
 

          mu_min = -10d0
          mu_max = 10d0
          mu_delta = 0d0
          call chemical_potential_tuning(mu_min, mu_max, rho, tolerance, N, beta, mu, evalup, evaldn,& 
            evecup, evecdn, mu_delta)
          mu = mu + mu_delta
          call compute_density(N, beta, mu, evalup, evaldn, evecup, evecdn, newnup, newndn)
      enddo


!     SET UP HAMILTONIAN
!     SET UP ELECTRON HOPPING
!******************************************************************
!     Periodic boundary condition (PBC) in both x and y directions
!******************************************************************
      ! if (if_pinning_field == .false.) then
      !     do ix = 1, Nx
      !         do iy = 1, Ny
      !             i =  ix + (iy - 1) * Nx
      !             j =  xplus(ix) + (iy - 1) * Nx
      !             Hup(i,j) = -tx
      !             Hdn(i,j) = -tx

      !             j = xminus(ix) + (iy - 1) * Nx
      !             Hup(i,j) = -tx
      !             Hdn(i,j) = -tx

      !             j = ix + (yplus(iy) - 1) * Nx
      !             Hup(i,j) = -ty
      !             Hdn(i,j) = -ty

      !             j = ix + (yminus(iy) - 1) * Nx
      !             Hup(i,j) = -ty
      !             Hdn(i,j) = -ty

      !             jprime = xplus(ix) + (yplus(iy) - 1) * Nx 
      !             Hup(i,jprime) = -tprime
      !             Hdn(i,jprime) = -tprime

      !             jprime = xminus(ix) + (yplus(iy) - 1) * Nx
      !             Hup(i,jprime) = -tprime
      !             Hdn(i,jprime) = -tprime

      !             jprime = xminus(ix) + (yminus(iy) - 1) * Nx
      !             Hup(i,jprime) = -tprime
      !             Hdn(i,jprime) = -tprime

      !             jprime = xplus(ix) + (yminus(iy) - 1) * Nx
      !             Hup(i,jprime) = -tprime
      !             Hdn(i,jprime) = -tprime
      !         end do
      !     end do
! !******************************************************************
! !     Open boundary condition (OBC) in x direction
! !     Periodic boundary condition (PBC) in y direction
! !******************************************************************
!       else if (if_pinning_field == .true.) then
!           do ix = 1, Nx
!               do iy = 1, Ny
!                   if (ix == 1) then
!                       i = ix + (iy - 1) * Nx
!                       j = xplus(ix) + (iy - 1) * Nx
!                       Hup(i,j) = -tx
!                       Hdn(i,j) = -tx

!                       j = ix + (yplus(iy) - 1) * Nx
!                       Hup(i,j) = -ty
!                       Hdn(i,j) = -ty

!                       j = ix + (yminus(iy) - 1) * Nx
!                       Hup(i,j) = -ty
!                       Hdn(i,j) = -ty

!                       jprime = xplus(ix) + (yplus(iy) - 1) * Nx
!                       Hup(i,jprime) = -tprime
!                       Hdn(i,jprime) = -tprime

!                       jprime = xplus(ix) + (yminus(iy) - 1) * Nx
!                       Hup(i,jprime) = -tprime
!                       Hdn(i,jprime) = -tprime
!                   else if (ix == Nx) then
!                       i = ix + (iy - 1) * Nx
!                       j = xminus(ix) + (iy - 1) * Nx
!                       Hup(i,j) = -tx
!                       Hdn(i,j) = -tx

!                       j = ix + (yplus(iy) - 1) * Nx
!                       Hup(i,j) = -ty
!                       Hdn(i,j) = -ty

!                       j = ix + (yminus(iy) - 1) * Nx
!                       Hup(i,j) = -ty
!                       Hdn(i,j) = -ty

!                       jprime = xminus(ix) + (yplus(iy) - 1) * Nx
!                       Hup(i,jprime) = -tprime
!                       Hdn(i,jprime) = -tprime

!                       jprime = xminus(ix) + (yminus(iy) - 1) * Nx
!                       Hup(i,jprime) = -tprime
!                       Hdn(i,jprime) = -tprime
!                   else
!                       i = ix + (iy - 1) * Nx
!                       j = xplus(ix) + (iy - 1) * Nx
!                       Hup(i,j) = -tx
!                       Hdn(i,j) = -tx

!                       j = xminus(ix) + (iy - 1) * Nx
!                       Hup(i,j) = -tx
!                       Hdn(i,j) = -tx

!                       j = ix + (yplus(iy) - 1) * Nx
!                       Hup(i,j) = -ty
!                       Hdn(i,j) = -ty

!                       j = ix + (yminus(iy) - 1) * Nx
!                       Hup(i,j) = -ty
!                       Hdn(i,j) = -ty

!                       jprime = xplus(ix) + (yplus(iy) - 1) * Nx
!                       Hup(i,jprime) = -tprime
!                       Hdn(i,jprime) = -tprime

!                       jprime = xminus(ix) + (yplus(iy) - 1) * Nx
!                       Hup(i,jprime) = -tprime
!                       Hdn(i,jprime) = -tprime

!                       jprime = xminus(ix) + (yminus(iy) - 1) * Nx
!                       Hup(i,jprime) = -tprime
!                       Hdn(i,jprime) = -tprime

!                       jprime = xplus(ix) + (yminus(iy) - 1) * Nx
!                       Hup(i,jprime) = -tprime
!                       Hdn(i,jprime) = -tprime
!                   end if
!               end do
!           end do
!       end if
! !******************************************************************
! !******************************************************************


!******************************************************************
!     Check the Hamiltoian setup
!******************************************************************
!       write(*, *) ' Initial Hamiltonian for spin-up electrons '
!       do i = 1, N
!           write(*, "(16f8.3)") (Hup(i,j), j=1,N)
!       end do

!       write(*, *) ' Initial Hamiltonian for spin-down electrons '
!       do i = 1, N
!           write(*, "(16f8.3)") (Hdn(i,j), j=1,N)
!       end do
!******************************************************************
!******************************************************************
! !     INITIALIZE DENSITIES
!       write (68, *) ' '
!       write (68, *) ' Initial Electron Density '  

!*****************************************************************
!     Random initial distribution used for test
!*****************************************************************   
!        do i = 1, N
!            nup(i) = scale*ran2(iran)
!            ndn(i) = scale*ran2(iran)
!            write (68, "(i6, 2f12.6)") i,nup(i),ndn(i)
!        end do
!*****************************************************************
!*****************************************************************

!*****************************************************************
!     Antiferromagnetic order initialization
!*****************************************************************
!       do i = 1, N
!           tmp_i = mod(i, Nx)
!           tmp_j = int((i-1)/Nx+1)
! !           write(*, *) tmp_j
!            if (mod(tmp_i+tmp_j, 2)==0) then
!              nup(i) = 1.0d0*rho
!              ndn(i) = 0.0d0*rho
!            else
!              nup(i) = 0.0d0*rho
!              ndn(i) = 1.0d0*rho
!            end if
!           write (68, "(i6, 2f12.6)") i, nup(i), ndn(i)
!       end do
!*****************************************************************
!*****************************************************************


! !***********************************************************************
! !     Read initial electron distribution from an input file
! !***********************************************************************
!             beta_prev = beta - beta_step
!             beta_prev_int = int(beta_prev + 1.0d-8)
!             beta_prev_float = mod(nint(beta_prev * 10), 10)

!             if (beta_prev_int>=10 .and. abs(tprime_float)>=10) then
!               formatting_string = "(a, a, i1, a, i2, a, i1, a, i1, a, i2, a, i1, a)"
!             else if (beta_prev_int>=10 .and. abs(tprime_float)<10) then
!               formatting_string = "(a, a, i1, a, i1, a, i1, a, i1, a, i2, a, i1, a)"
!             else if (beta_prev_int<10 .and. abs(tprime_float)>=10) then
!               formatting_string = "(a, a, i1, a, i2, a, i1, a, i1, a, i1, a, i1, a)"
!             else
!               formatting_string = "(a, a, i1, a, i1, a, i1, a, i1, a, i1, a, i1, a)"
!             endif

!             write(input_name,  formatting_string) 'Data/p02_test/tprime-0p2/rho0p875/OBC/L16W4/Spin_', &
!                   tprime_string, tprime_int, tmp_string, abs(tprime_float), '_U', Uint, 'p', Ufloat, '_beta', beta_prev_int, 'p', beta_prev_float, '.dat'
!             write(*, *) input_name

!             open(unit=72,    status='old',    file=input_name)
!             do i = 1, N
!               read(72, *) tmp, nup(i), ndn(i)
!               write(*, *) tmp, nup(i), ndn(i)
!             end do
!             close(72)
! !***********************************************************************
! !***********************************************************************


! **********************************************************************
!     Initialize Stripe Order Phase with Different Domain Size
! **********************************************************************
!       if (index - number_of_trial_states < 1.0d-8) then
!           write(*, *) "************************************************"
!           write(*, *) "wavelength = ", index, 
!      1    "stripe order initialization"
!           write(*, *) if_pinning_field, pin_amplitude
!           write(*, *) "************************************************"
!             do i = 1, N
!                 tmp_i = mod(i, Nx)
!                 tmp_j = int((i - 1)/Nx) + 1
!                 phase_shift = int((i - 1)/index)
!                 if (mod(tmp_i+tmp_j+phase_shift, 2)==1) then
!                   nup(i) = 1.0d0*rho
!                   ndn(i) = 0.0d0
!                 else
!                   nup(i) = 0.0d0
!                   ndn(i) = 1.0d0*rho
!                 end if
!                 write(*, *) i, nup(i)-ndn(i)
!             end do
!       else
!             do i = 1, N
!                 tmp_i = mod(i, Nx)
!                 tmp_j = int((i-1)/Nx+1)
!                 if (mod(tmp_i+tmp_j, 2)==0) then
!                     nup(i) = 1.0d0*rho
!                     ndn(i) = 0.0d0*rho
!                 else
!                     nup(i) = 0.0d0*rho
!                     ndn(i) = 1.0d0*rho
!                 end if
!                 write (*, *) i, nup(i)-ndn(i)
!                 write (68, "(i6, 2f12.6)") i, nup(i), ndn(i)
!             end do
!       end if
! **********************************************************************
! **********************************************************************


!*******************************************************************************************
!     DBG:  CHECK SYMMETRIC

!      do 180 i=1,N
!      do 170 j=1,N
!          if (Hup(i,j).ne.Hup(j,i)) then
!              write (6,*) 'Hup symmetry violation at'
!              write (6,990) i,j,Hup(i,j),Hup(j,i)
!              write (6,*) 'stopping'
!              stop
!          endif
!          if (Hdn(i,j).ne.Hdn(j,i)) then
!              write (6,*) 'Hdn symmetry violation at'
!              write (6,990) i,j,Hdn(i,j),Hdn(j,i)
!              write (6,*) 'stopping'
!              stop
!          endif
!170   continue
!180   continue
!990   format(2i6,2f12.6)
!********************************************************************************************


!     SELF CONSISTENT LOOP
!       do it=1, Nit
!       write (6, "(' iteration= ', i6)") it

!     DBG:
!      write (68,*) ' '
!      write (68,*) ' eigenvalues of Hup,Hdn '
!      do 200 i=1,N
!          write (68,991) i,evalup(i),evaldn(i)
! 200   continue
!      write (68,*) ' '
!      write (68,*) ' eigenvectors of Hup,Hdn '
!      do 220 j=1,N
!      do 210 i=1,N
!          write (68,992) i,j,evecup(i,j),evecdn(i,j)
! 210   continue
! 220   continue
! 991   format(i6,2f16.8)
! 992   format(2i6,2f16.8)

!     DBG 2:  DIAGONALIZE

!      write (68,*) ' '
!      write (68,*) ' Check D = S H S  '
!      do 250 i=1,N
!      do 240 j=1,N
!           tempup(i,j)=0.d0
!           tempdn(i,j)=0.d0
!           do 230 k=1,N
!                tempup(i,j)=tempup(i,j)+Hup(i,k)*evecup(k,j)
!                tempdn(i,j)=tempdn(i,j)+Hdn(i,k)*evecdn(k,j)
! 230        continue
! 240   continue
! 250   continue
!      do 280 i=1,N
!      do 270 j=1,N
!           checkidup(i,j)=0.d0
!           checkiddn(i,j)=0.d0
!           do 260 k=1,N
!                checkidup(i,j)=checkidup(i,j)+evecup(k,i)*tempup(k,j)
!                checkiddn(i,j)=checkiddn(i,j)+evecdn(k,i)*tempdn(k,j)
! 260        continue
!      write (68,992) i,j,checkidup(i,j),checkiddn(i,j)
! 270   continue
! 280   continue


      ! write(*, *) mu
! ***********************************************************************
!     Add Perturbation to electron densities
! ***********************************************************************
      if ((mod(it, 50) .eq. 0) .and. (it .le. 300)) then
            tmp_total_electrons=0.0d0
            tmp_shift = 0.0d0
            tmp1 = 0.0d0
            tmp2 = 0.0d0
            do i = 1, N
                call random_number(tmp1)
                call random_number(tmp2)
                tmp1 = 2 * tmp1 - 1.0d0
                tmp2 = 2 * tmp2 - 1.0d0
                newnup(i) = newnup(i) + annealing_amp * tmp1
                newndn(i) = newndn(i) + annealing_amp * tmp2

                if (newnup(i) .gt. 1.0d0) then
                    newnup(i) = 1.0d0
                else if (newnup(i) .lt. 0.0d0) then
                    newnup(i) = 0.0d0
                end if

                if (newndn(i) .gt. 1.0d0) then
                    newndn(i) = 1.0d0
                else if (newndn(i) .lt. 0.0d0) then
                    newndn(i) = 0.0d0
                end if

                tmp_total_electrons = tmp_total_electrons + newnup(i)
                tmp_total_electrons = tmp_total_electrons + newndn(i)
            end do

            tmp_shift = (tmp_total_electrons - N * rho) / N

            do i = 1, N
                newnup(i) = newnup(i) - tmp_shift
                newndn(i) = newndn(i) - tmp_shift
            end do
      end if

      ! do i=1,N
      !     newnup(i) = relax*newnup(i) + (1.d0-relax)*nup(i)
      !     newndn(i) = relax*newndn(i) + (1.d0-relax)*ndn(i)
      !     Hup(i,i) = Hup(i,i) + U * ( newndn(i) - ndn(i) )
      !     Hdn(i,i) = Hdn(i,i) + U * ( newnup(i) - nup(i) )
      !     nup(i) = newnup(i)
      !     ndn(i) = newndn(i)
      ! end do
      
      ! rho_relaxed=0.d0
      ! do i = 1, N
      !   rho_relaxed = rho_relaxed + nup(i) + ndn(i)
      ! end do
      ! rho_relaxed = rho_relaxed / dble(N)
      ! write(71, "(i6, f12.6)") it, rho_relaxed

! **********************************************************************
!     Compute the grand free energy
! **********************************************************************
      Free_Energy = 0.0d0
      do i = 1, N
          Free_Energy=Free_Energy+log(1.0d0+exp(-beta*(evalup(i)-mu))) 
          Free_Energy=Free_Energy+log(1.0d0+exp(-beta*(evaldn(i)-mu)))
          Free_Energy=Free_Energy+beta*U*nup(i)*ndn(i) 
      end do 
      Free_Energy = -1.0d0 / beta * Free_Energy

      totalGF(it,index)=Free_Energy/dfloat(N)
      totalF(it,index)=Free_Energy/dfloat(N)+mu*rho
      write(*, *) E/dfloat(N), Free_Energy/dfloat(N)+mu*rho, mu
! **********************************************************************
! **********************************************************************
      end do
      
      do i = 1, N
        upElectron(i,index) = nup(i)
        dnElectron(i,index) = ndn(i)
      end do

!     Ending point for the different trial states
      end do

!***********************************************************************
!     Select and store the configuration that has the lowest
!     grand free energy.
!***********************************************************************
      if (number_of_trial_states-1<1E-8) then
        do i=1,N
          write(69, "(i6, 2f12.6)") i, upElectron(i,1), dnElectron(i,1)
        end do

        do i=1,Nit
          write(70, "(i6, f16.12)") i, totalE(i, 1)
          write(80, "(i6, f16.12)") i, totalF(i, 1)
          write(81, "(i6, f16.12)") i, totalGF(i, 1)
        end do
      else
        optimal_index = 1
        do index = 2,number_of_trial_states
          if (totalF(Nit, index) < totalF(Nit, optimal_index)) then
            optimal_index = index
          end if
        end do

        do i=1,N
          write(69, "(i6, 2f12.6)") i, upElectron(i, optimal_index), dnElectron(i, optimal_index)
        end do

        do i=1,Nit
          write(70, "(i6, f16.12)") i, totalE(i, optimal_index)
          write(80, "(i6, f16.12)") i, totalF(i, optimal_index)
          write(81, "(i6, f16.12)") i, totalGF(i, optimal_index)
        end do
      end if


      do index=1,number_of_trial_states
        do i=1,N
          write(69, "(i6, 2f12.6)") i, upElectron(i,index),  
     1    dnElectron(i,index)
        end do
      end do

      do tmpInd = 1, number_of_trial_states
        do i=1,Nit
          write(70, "(i6, f16.12)") i, totalE(i, tmpInd)
          write(80, "(i6, f16.12)") i, totalF(i, tmpInd)
          write(81, "(i6, f16.12)") i, totalGF(i, tmpInd)
        end do
      end do

!     Ending point for the effective inverse temperature loop
      end do 

!     Ending point for the effective Coulomb interaction loop      
      end do

!     Ending point for the tprime loop
      end do

contains
      subroutine split(input_number, amplifier, tmp_int, tmp_float)
            real(kind=precision), intent(in) :: input_number
            integer, intent(in) :: amplifier
            integer, intent(out) :: tmp_int
            integer, intent(out) :: tmp_float

            tmp_int = int(input_number + 1.0d-8)
            tmp_float = mod(nint(input_number * amplifier), amplifier)
      end subroutine split

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

      
      
      ! Initialize the spin-up and spin-down electron distributions 
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
              print '(i6, 2f12.6)', ind, tmp_Hup(ind, ind), tmp_Hdn(ind, ind)
          end do

          if (pinning_strength > 1e-8) then
              do ind = 1, tmp_Ny
                  tmp_ind = (ind - 1) * tmp_Nx + 1
                  tmp_Hup(tmp_ind, tmp_ind) = tmp_Hup(tmp_ind, tmp_ind) + (-1)**(ind) * 0.5d0 * pinning_strength
                  tmp_Hdn(tmp_ind, tmp_ind) = tmp_Hdn(tmp_ind, tmp_ind) + (-1)**(ind + 1) * 0.5d0 * pinning_strength 
                  ! print '(i6, f12.6)', ind, tmp_Hup(tmp_ind, tmp_ind) - tmp_Hdn(tmp_ind, tmp_ind)
                  print '(i6, 2f12.6)', ind, tmp_Hup(tmp_ind, tmp_ind), tmp_Hdn(tmp_ind, tmp_ind)
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
              tmp_energy = tmp_energy + tmp_fermiup + tmp_fermidn
              tmp_energy = tmp_energy - tmp_U * tmp_nup(ind) * tmp_ndn(ind)
          enddo
          print "('iteration, E/N = ', i6, f16.8)", tmp_iteration, tmp_energy / dfloat(tmp_Nsites)
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



      ! Tune the chemical potential 
      subroutine chemical_potential_tuning(lower_bound, upper_bound, target_density,& 
        tmp_tolerance, tmp_Nsites, tmp_beta, tmp_mu,& 
        tmp_evalup, tmp_evaldn, tmp_evecup, tmp_evecdn, tmp_delta_mu)
          integer, intent(in) :: tmp_Nsites
          real(kind=precision), intent(in) :: tmp_tolerance, target_density
          real(kind=precision), intent(in) :: tmp_beta, tmp_mu
          real(kind=precision), intent(in) :: tmp_evalup(tmp_Nsites), tmp_evaldn(tmp_Nsites)
          real(kind=precision), intent(in) :: tmp_evecup(tmp_Nsites, tmp_Nsites), tmp_evecdn(tmp_Nsites, tmp_Nsites)
          real(kind=precision), intent(out) :: lower_bound, upper_bound, tmp_delta_mu
          real(kind=precision) :: density_diff, tmp_density
          real(kind=precision) :: tmp_nup(tmp_Nsites), tmp_ndn(tmp_Nsites)
          integer :: ind

          density_diff = 1.0d6
          do while (density_diff > tmp_tolerance)
              tmp_delta_mu = 0.5d0 * (lower_bound + upper_bound)
              tmp_density = 0.d0

              call compute_density(tmp_Nsites, tmp_beta, tmp_mu + tmp_delta_mu,& 
                tmp_evalup, tmp_evaldn, tmp_evecup, tmp_evecdn, tmp_nup, tmp_ndn)

              do ind = 1, tmp_Nsites
                  tmp_density = tmp_density + tmp_nup(ind) + tmp_ndn(ind)
              enddo

              tmp_density = tmp_density / dfloat(tmp_Nsites)
              density_diff = abs(tmp_density - target_density)
              print '(f12.6)', tmp_density
              
              if (tmp_density > target_density) then
                  upper_bound = tmp_delta_mu
              else if (tmp_density < target_density) then
                  lower_bound = tmp_delta_mu
              end if
          enddo
      end subroutine chemical_potential_tuning

      subroutine density_mixture(tmp_Nsites, tmp_nup, tmp_ndn, tmp_newnup, tmp_newndn, tmp_U, tmp_Hup, tmp_Hdn, tmp_relax)
            integer, intent(in) :: tmp_Nsites
            real(kind=precision), intent(in) :: tmp_U, tmp_relax
            real(kind=precision), intent(in) :: tmp_nup(tmp_Nsites), tmp_ndn(tmp_Nsites)
            real(kind=precision), intent(in) :: tmp_newnup(tmp_Nsites), tmp_newndn(tmp_Nsites)
            real(kind=precision), intent(out) :: tmp_Hup(tmp_Nsites, tmp_Nsites), tmp_Udn(tmp_Nsites, tmp_Nsites)

            integer :: ind
            do i=1,N
          newnup(i) = relax*newnup(i) + (1.d0-relax)*nup(i)
          newndn(i) = relax*newndn(i) + (1.d0-relax)*ndn(i)
          Hup(i,i) = Hup(i,i) + U * ( newndn(i) - ndn(i) )
          Hdn(i,i) = Hdn(i,i) + U * ( newnup(i) - nup(i) )
          nup(i) = newnup(i)
          ndn(i) = newndn(i)
      end do
      
      rho_relaxed=0.d0
      do i = 1, N
        rho_relaxed = rho_relaxed + nup(i) + ndn(i)
      end do
      rho_relaxed = rho_relaxed / dble(N)
      write(71, "(i6, f12.6)") it, rho_relaxed
      end subroutine density_mixture
end program IHF