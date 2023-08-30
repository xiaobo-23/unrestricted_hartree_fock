!     01/25/2023
!     UHF FOR THE HUBBARD MODEL WITH NONZERO t'.

!     VERSION 4 USES LAPACK SO WE CAN CHECK LARGE LATTICES BETTER.
!     Nx AND Ny NO LONGER INPUT BUT SET IN PARAMETER STATEMENTS
!     [x]  MATCHES VERSION 3 *BEFORE* GOING TO LAPACK
!     [x]  COMPILES -Wall *BEFORE* GOING TO LAPACK
!     [x]  RUNS -fbounds-check *BEFORE* GOING TO LAPACK
!     [x]  MATCHES VERSION 3 *AFTER* GOING TO LAPACK
!                   8x8 lattice     16x16 lattice     32x32 lattice
!                    1000 its        200 its
!          ih3.f:    9.92 sec        129.1 sec
!          ih4.f:    7.85 sec         17.9 sec          1413 sec 
!     [ ]  AGREES WITH HIRSCH PHASE DIAGRAM.
!          https://journals.aps.org/prb/pdf/10.1103/PhysRevB.31.4403
!          WITH NNN HOPPING:
!          https://journals.aps.org/prb/pdf/10.1103/PhysRevB.35.3359

!     VERSION 3 WILL HAVE SELF CONSISTENCY.  
!     [x]  CONVERGES U NONZERO 
!     [x]  HALF-FILLING AT mu=0
!     [x]  AF AT mu=0
!         8x8 lattice U=4 mu=2  beta=1.37
!           1    0.569172    0.569177       0.430828    0.430823
!           2    0.430828    0.430823       0.569172    0.569177
!           3    0.569172    0.569177       0.430828    0.430823
!           4    0.430828    0.430823       0.569172    0.569177
!                      [...]

!     VERSION 2 CALCULATION OF ENERGY AND RECALCULATION OF DENSITY
!          NO LOOP
!     [x]  ENERGY CORRECT AT t=U=0  mu=-1.3 beta=0.7 Nx=8 Ny=4
!           E =       23.87838646
!             =   2 * 32 * 1.3 * [ 1 / (exp(1.3*0.7) + 1) ]
!     [x]  CONVERGES SINGLE ITERATION t=U=0
!             mu=-1.3 beta=0.7 Nx=8 Ny=4
!           old/new densities up and dn
!               1    0.096076    0.287000       0.010154    0.287000
!               2    0.017264    0.287000       0.009535    0.287000
!                      [...]
!              [ 1 / (exp(1.3*0.7) + 1) ] = 0.287000
!     [x]  CONVERGES SINGLE ITERATION t=1 U=0 TO HOMOGENEOUS DENSITY
!             tx=ty=1.2  mu=-1.3 beta=0.7 Nx=8 Ny=4
!                   old/new densities up and dn
!              1    0.096076    0.355289       0.010154    0.355289
!              2    0.017264    0.355289       0.009535    0.355289
!                      [...]
!             tx=ty=1.2  mu=-0.0 beta=0.7 Nx=8 Ny=4
!                   old/new densities up and dn
!              1    0.096076    0.500000       0.010154    0.500000
!              2    0.017264    0.500000       0.009535    0.500000
!                      [...]

!     VERSION 1 HAS NO SELF CONSISTENCY.  
!        JUST SETS UP MATRIX AND DIAGONALIZES ONCE
!     [x]  WRITTEN
!     [x]  COMPILES -Wall
!     [x] RUNS -fbounds-check 
!         CHECKED AT NP=256  Nx=Ny=16 
!     GIVES SENSIBLE RESULTS
!     [x]  EIGENVALUES CORRECT AT U=0
!           6   6   1.0  0.1   0.0  10.0   1.2      0.6
!           Nx  Ny   tx   ty    U     mu   beta    scale
!             eigenvalues of Hup,Hdn
!                1    -12.20000000    -12.20000000
!                2    -12.10000000    -12.10000000
!                3    -12.10000000    -12.10000000
!                4    -11.90000000    -11.90000000
!                5    -11.90000000    -11.90000000
!                6    -11.80000000    -11.80000000
!                7    -11.20000000    -11.20000000
!                      [...]
!               35     -7.90000000     -7.90000000
!               36     -7.80000000     -7.80000000
!     [x]  EIGENVALUES CORRECT AT U=0
!           8   4   1.0  0.1   0.0  10.0   1.2      0.6
!           Nx  Ny   tx   ty    U     mu   beta    scale
!             eigenvalues of Hup,Hdn
!                1    -12.20000000    -12.20000000
!                2    -12.00000000    -12.00000000
!                3    -12.00000000    -12.00000000
!                4    -11.80000000    -11.80000000
!                5    -11.61421356    -11.61421356
!     [x]  CAN DIAGONALIZE WITH EIGENVECTOR MATRIX
!           6   6   1.0  0.1   0.0  10.0   1.2      0.6
!           Nx  Ny   tx   ty    U     mu   beta    scale
!              Check D = S H S
!                  1     1    -12.20000000    -12.20000000
!                  1     2      0.00000000      0.00000000
!                  1     3      0.00000000      0.00000000
!                  1     4      0.00000000      0.00000000
!                      [...]
!                  2     1      0.00000000      0.00000000
!                  2     2    -12.14142136    -12.14142136
!                  2     3      0.00000000      0.00000000
!                      [...]
!     [x]  CAN DIAGONALIZE WITH EIGENVECTOR MATRIX
!           8   4   1.0  0.1   0.0  10.0   1.2      0.6
!           Nx  Ny   tx   ty    U     mu   beta    scale
!              Check D = S H S
!                      [...]
!                  5     4     -0.00000000     -0.00000000
!                  5     5    -11.61421356    -11.61421356
!                  5     6     -0.00000000     -0.00000000
!                      [...]
!     [x]  EIGENVALUES CORRECT AT t=0  U=3  mu=-10
!          LARGEST  ndn = 0.194618
!          LARGEST  Eup = -9.416147 = -10+3(0.194618)
!          LARGEST  nup = 0.197427
!          LARGEST  Edn = -9.407720 = -10+3(0.197427)
      
      program IHF

      use iso_fortran_env
      implicit none

      integer, parameter :: precision=8
      integer, parameter :: Nx=16, Ny=4
      integer, parameter :: N = Nx * Ny
      integer, parameter ::  iran=20000
!     integer nrot,NP
      
      integer :: ix, iy 
      integer :: xplus(N), xminus(N) 
      integer :: yplus(N), yminus(N)
      integer :: i, j, jprime
      integer :: tmp_i, tmp_j, site
      integer :: tmp_sign, tmp_wave, tmp_loc
      integer :: tmp

!     integer k
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

!     real*8 tempup(N,N),checkidup(N,N)
!     real*8 tempdn(N,N),checkiddn(N,N)

!******************************************************************
!     Add next-nearest neighbor hopping into the code
!******************************************************************
      real(kind=precision), parameter :: tx=1.0d0
      real(kind=precision), parameter :: ty=1.0d0
      real(kind=precision), parameter :: scale=0.5d0
      real(kind=precision) :: mu
      real(kind=precision) :: ran2
!*****************************************************************
!*****************************************************************


!*****************************************************************
!     Add pinning field into the code
!     Add electron density
!*****************************************************************
      real(kind=precision), parameter :: pin_amplitude=0.0d0
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
      real(kind=precision) :: rho1, rho2
      real(kind=precision), parameter :: tolerance=1.0d-8
!*****************************************************************
!*****************************************************************


C***********************************************************************
C     Add temperature bounds and steps
C***********************************************************************
      real(precision), parameter :: beta_min=50.0d0
      real(precision), parameter :: beta_max=50.0d0
      real(precision), parameter :: beta_step=0.2d0
      integer :: beta_int, beta_float, beta_prev_int, beta_prev_float
      real(precision) :: beta, beta_prev
C***********************************************************************
C***********************************************************************


C***********************************************************************
C     Add effective U bounds and steps
C***********************************************************************
      real(precision), parameter :: Umin=2.8d0
      real(precision), parameter :: Umax=2.8d0
      real(precision), parameter :: Ustep=0.1d0
      integer :: Uint, Ufloat
      real(precision) :: U
C***********************************************************************
C***********************************************************************


!*****************************************************************
!     Add next nearest-neighbor hopping t' bounds and steps
!*****************************************************************
      real(precision), parameter :: tprime_min=0d0
      real(precision), parameter :: tprime_max=0d0
      real(precision), parameter :: tprime_step=0.05d0
      integer :: tprime_int, tprime_float
      real(precision) :: tprime
!*****************************************************************
!*****************************************************************

C***********************************************************************
C     Initialize different domain sizes as the starting point of 
C     Hartree-Fock solution
C***********************************************************************
      integer :: phase_shift
      integer :: domain_size
      integer :: index
      integer :: tmpInd
      integer :: optimal_index
      integer, parameter :: number_of_trial_states=1

!***********************************************************************
!***********************************************************************
      integer :: it 
      integer, parameter :: Nit = 3500
      real(kind=precision), parameter :: relax = 0.7d0 
      integer :: INFO
      integer, parameter :: N6 = 6 * N
      real(kind=precision) :: work(N6)
      character(len = :), allocatable :: formatting_string
      character(len = 100) outname
      character(len = 100) input_name
      character(len = 100) spin
      character(len = 100) free_energy_str
      character(len = 100) grand_str
      character(len = 100) e_hist
      character(len = 100) e_density
      character(len = 100) eigen
      character(len = 100) path
      character(len = 100) parameter_string
      character(len = :), allocatable :: tprime_string
      character(len = :), allocatable :: tmp_string

C **********************************************************************
C **********************************************************************
      real(kind=precision) :: totalE(Nit, number_of_trial_states)
      real(kind=precision) :: totalF(Nit, number_of_trial_states)
      real(kind=precision) :: totalGF(Nit, number_of_trial_states)
      real(kind=precision) :: upElectron(N, number_of_trial_states) 
      real(kind=precision) :: dnElectron(N, number_of_trial_states)
C **********************************************************************
C **********************************************************************

!********************************************************************
!     Output files.
!     Spin density file can be used as CP-QMC input file directly.
!********************************************************************
      path='Data/PBC/tprime0/rho0p875/p0/L16W4/'

      do tprime = tprime_min, tprime_max, tprime_step
        tprime_int = int(tprime + 1.0d-7)
        tprime_float = mod(nint(tprime * 100), 100)
          
          do U = Umin, Umax, Ustep
            Uint = int(U + 1.0d-7)
            Ufloat = mod(nint(U * 10), 10)
            do beta = beta_min, beta_max, beta_step
C             do beta = beta_max, beta_min, beta_step
              write(*, *) U, Ustep, beta, beta_step
              beta_int = int(beta + 1.0d-7)
              beta_float = mod(nint(beta * 10), 10)

              if (beta_int>=10 .and. abs(tprime_float)>=10) then
                formatting_string = 
     1          "(a, i1, a, i2, a, i1, a, i1, a, i2, a, i1, a)"
              else if (beta_int>=10 .and. abs(tprime_float)<10) then
                formatting_string = 
     1          "(a, i1, a, i1, a, i1, a, i1, a, i2, a, i1, a)"
              else if (beta_int<10 .and. abs(tprime_float)>=10) then
                formatting_string = 
     1         "(a, i1, a, i2, a, i1, a, i1, a, i1, a, i1, a)"
              else
                formatting_string = 
     1          "(a, i1, a, i1, a, i1, a, i1, a, i1, a, i1, a)"
              endif

              if (tprime >= 0) then
                tprime_string = 'tprime'
              else
                tprime_string = 'tprime-'
              end if

              if (abs(tprime_float) >= 10) then
                tmp_string = 'p'
              else
                tmp_string = 'p0'
              end if

        write(parameter_string, formatting_string) 
     1  tprime_string, tprime_int, tmp_string, abs(tprime_float),
     1  '_U', Uint, 'p', Ufloat, '_beta', beta_int, 'p', 
     1  beta_float, '.dat' 


        ! Optimization log
        outname = trim(path)//trim('Results_')//trim(parameter_string) 
        open(unit=68, status='replace', file=outname)

        ! Converged spin density
        spin=trim(path)//trim("Spin_")//trim(parameter_string)
        open(unit=69, status='replace', file=spin)
            
        ! Energy histogram
        e_hist=trim(path)//trim("Energy_")//trim(parameter_string)
        open(unit=70, status='replace', file=e_hist)

        ! Converged charge density
        e_density=trim(path)//trim("Charge_")//trim(parameter_string)
        open(unit=71, status="replace", file=e_density)

        ! Free energy histogram
        free_energy_str=trim(path)//trim("FE_")//trim(parameter_string)
        open(unit=80, status='replace', file=free_energy_str)

        ! Gran potential histogram
        grand_str=trim(path)//trim("GE_")//trim(parameter_string)
        open(unit=81, status='replace', file=grand_str)

        ! Converged eigenvalues
        eigen=trim(path)//trim("Eigen_")//trim(parameter_string)
        open(unit=82, status='replace', file=eigen)

!*******************************************************************
!     No input file needed in this verion.
!*******************************************************************
      if (pin_amplitude > 1.0d-8) then
        if_pinning_field = .true.
      else
        if_pinning_field = .false.
      endif
      mu=0.d0

      write (68, *) ' ********************************************* '
      write (68, *) ' Inhomogeneous Hartree-Fock'
      write (68, *) ' Version 7 '
      write (68, *) ' ********************************************* '
      write (68, *) ' Lattice size: Nx,  Ny '
      write (68, "(2i12)")  Nx, Ny
      write (68, *) ' Nearest Hopping: tx, ty '
      write (68, "(2f12.6)") tx, ty
      write (68, *) ' Next-Nearest Hopping tprime'
      write (68, "(f12.6)") tprime
      write (68, *) ' Effective Coulumb Interaction '
      write (68, "(f12.6)") U
      write (68, *) ' Original Chemical Potential '
      write (68, "(f12.6)") mu
      write (68, *) ' Inverse Temperature'
      write (68, "(f12.6)") beta
      write (68, *) ' Pinning Field Strength '
      write (68, "(f12.6)") pin_amplitude
      write (68, *) ' Random Initialization Scale '
      write (68, "(f12.6)") scale
      write (68, *) ' Random Number Seed'
      write (68, "(i12)") iran
      write (68, *) ' Number of Iterations '
      write (68, "(i12)") Nit
      write (68, *) ' Relaxation Rate '
      write (68, "(f12.6)") relax
      write (68, *) ' If Pinning Is Turned On'
      write (68, *)   if_pinning_field
      write (68, *) ' Output File Name '
      write (68, "(a60)") outname

      do index = 1, number_of_trial_states
          write(*, *) 'Number of trial states is:', index

!     ZERO HAMILTONIAN.  ASSIGN -mu TO DIAGONAL.
!****************************************************************
!     Construct the Hamiltonian; extrapolate -mu from the matrix
!     Add pinning field
!     Fix boundary conditions based on the pinning field 
!     Tune chemical potential during each iteration
!****************************************************************
      do i=1,N
          do j=1,N
              Hup(i,j) = 0.d0
              Hdn(i,j) = 0.d0
          end do
C            Hup(i,i) = -mu
C            Hdn(i,i) = -mu
C            Hup(i,i) = -0.5*U
C            Hdn(i,i) = -0.5*U
      end do


!     IMPLEMENT PERIODIC BOUNDARY CONDITION
!     INDIRECT ADDRESSES
      do ix = 1, Nx
          xplus(ix) = ix + 1
          xminus(ix) = ix - 1
      end do
      xplus(Nx) = 1
      xminus(1)  = Nx

      do iy = 1, Ny
          yplus(iy) = iy + 1
          yminus(iy) = iy - 1
      end do
      yplus(Ny) = 1
      yminus(1)  = Ny


!     SET UP HAMILTONIAN
!     SET UP ELECTRON HOPPING
!******************************************************************
!     Periodic boundary condition (PBC) in both x and y directions
!******************************************************************
      if (if_pinning_field == .false.) then
          do ix = 1, Nx
              do iy = 1, Ny
                  i =  ix + (iy - 1) * Nx
                  j =  xplus(ix) + (iy - 1) * Nx
                  Hup(i,j) = -tx
                  Hdn(i,j) = -tx

                  j = xminus(ix) + (iy - 1) * Nx
                  Hup(i,j) = -tx
                  Hdn(i,j) = -tx

                  j = ix + (yplus(iy) - 1) * Nx
                  Hup(i,j) = -ty
                  Hdn(i,j) = -ty

                  j = ix + (yminus(iy) - 1) * Nx
                  Hup(i,j) = -ty
                  Hdn(i,j) = -ty

                  jprime = xplus(ix) + (yplus(iy) - 1) * Nx 
                  Hup(i,jprime) = -tprime
                  Hdn(i,jprime) = -tprime

                  jprime = xminus(ix) + (yplus(iy) - 1) * Nx
                  Hup(i,jprime) = -tprime
                  Hdn(i,jprime) = -tprime

                  jprime = xminus(ix) + (yminus(iy) - 1) * Nx
                  Hup(i,jprime) = -tprime
                  Hdn(i,jprime) = -tprime

                  jprime = xplus(ix) + (yminus(iy) - 1) * Nx
                  Hup(i,jprime) = -tprime
                  Hdn(i,jprime) = -tprime
              end do
          end do
!******************************************************************
!     Open boundary condition (OBC) in x direction
!     Anti-Periodic boundary condition (APBC) in y direction
!******************************************************************
      else if (if_pinning_field == .true.) then
          do ix = 1, Nx
              do iy = 1, Ny
                  if (ix == 1) then
                      i = ix + (iy - 1) * Nx
                      j = xplus(ix) + (iy - 1) * Nx
                      Hup(i,j) = -tx
                      Hdn(i,j) = -tx

                      if (iy == 1) then
                        j=ix+(yplus(iy)-1)*Nx
                        Hup(i,j)=-ty
                        Hdn(i,j)=-ty

                        j=ix+(yminus(iy)-1)*Nx
                        Hup(i,j)=ty
                        Hdn(i,j)=ty

                        jprime=xplus(ix)+(yplus(iy)-1)*Nx
                        Hup(i,jprime)=-tprime
                        Hdn(i,jprime)=-tprime

                        jprime=xplus(ix)+(yminus(iy)-1)*Nx
                        Hup(i,jprime)=tprime
                        Hdn(i,jprime)=tprime
                      else if (iy == Ny) then
                        j=ix+(yplus(iy)-1)*Nx
                        Hup(i,j)=ty
                        Hdn(i,j)=ty

                        j=ix+(yminus(iy)-1)*Nx
                        Hup(i,j)=-ty
                        Hdn(i,j)=-ty

                        jprime=xplus(ix)+(yplus(iy)-1)*Nx
                        Hup(i,jprime)=tprime
                        Hdn(i,jprime)=tprime

                        jprime=xplus(ix)+(yminus(iy)-1)*Nx
                        Hup(i,jprime)=-tprime
                        Hdn(i,jprime)=-tprime
                      else
                        j=ix+(yplus(iy)-1)*Nx
                        Hup(i,j)=-ty
                        Hdn(i,j)=-ty

                        j=ix+(yminus(iy)-1)*Nx
                        Hup(i,j)=-ty
                        Hdn(i,j)=-ty

                        jprime=xplus(ix)+(yplus(iy)-1)*Nx
                        Hup(i,jprime)=-tprime
                        Hdn(i,jprime)=-tprime

                        jprime=xplus(ix)+(yminus(iy)-1)*Nx
                        Hup(i,jprime)=-tprime
                        Hdn(i,jprime)=-tprime
                      endif
                  else if (ix == Nx) then
                      i = ix + (iy - 1) * Nx
                      j = xminus(ix) + (iy - 1) * Nx
                      Hup(i,j) = -tx
                      Hdn(i,j) = -tx

                      if (iy==1) then
                        j = ix + (yplus(iy) - 1) * Nx
                        Hup(i,j) = -ty
                        Hdn(i,j) = -ty

                        j = ix + (yminus(iy) - 1) * Nx
                        Hup(i,j) = ty
                        Hdn(i,j) = ty

                        jprime = xminus(ix) + (yplus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime

                        jprime = xminus(ix) + (yminus(iy) - 1) * Nx
                        Hup(i,jprime) = tprime
                        Hdn(i,jprime) = tprime
                      else if (iy==Ny) then
                        j = ix + (yplus(iy) - 1) * Nx
                        Hup(i,j)=ty
                        Hdn(i,j)=ty

                        j = ix + (yminus(iy) - 1) * Nx
                        Hup(i,j) = -ty
                        Hdn(i,j) = -ty

                        jprime = xminus(ix) + (yplus(iy) - 1) * Nx
                        Hup(i,jprime)=tprime
                        Hdn(i,jprime)=tprime

                        jprime = xminus(ix) + (yminus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime
                      else
                        j = ix + (yplus(iy) - 1) * Nx
                        Hup(i,j) = -ty
                        Hdn(i,j) = -ty

                        j = ix + (yminus(iy) - 1) * Nx
                        Hup(i,j) = -ty
                        Hdn(i,j) = -ty

                        jprime = xminus(ix) + (yplus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime

                        jprime = xminus(ix) + (yminus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime
                      endif
                  else
                      i = ix + (iy - 1) * Nx
                      j = xplus(ix) + (iy - 1) * Nx
                      Hup(i,j) = -tx
                      Hdn(i,j) = -tx

                      j = xminus(ix) + (iy - 1) * Nx
                      Hup(i,j) = -tx
                      Hdn(i,j) = -tx

                      if (iy==1) then
                        j = ix + (yplus(iy) - 1) * Nx
                        Hup(i,j) = -ty
                        Hdn(i,j) = -ty

                        j = ix + (yminus(iy) - 1) * Nx
                        Hup(i,j)=ty
                        Hdn(i,j)=ty

                        jprime = xplus(ix) + (yplus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime

                        jprime = xminus(ix) + (yplus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime

                        jprime = xminus(ix)+(yminus(iy)-1)*Nx
                        Hup(i,jprime)=tprime
                        Hdn(i,jprime)=tprime

                        jprime = xplus(ix) + (yminus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime
                      else if (iy==Ny) then
                        j = ix + (yplus(iy) - 1) * Nx
                        Hup(i,j)=ty
                        Hdn(i,j)=ty

                        j = ix + (yminus(iy) - 1) * Nx
                        Hup(i,j) = -ty
                        Hdn(i,j) = -ty

                        jprime = xplus(ix) + (yplus(iy) - 1) * Nx
                        Hup(i,jprime)=tprime
                        Hdn(i,jprime)=tprime

                        jprime = xminus(ix) + (yplus(iy) - 1) * Nx
                        Hup(i,jprime)=tprime
                        Hdn(i,jprime)=tprime

                        jprime = xminus(ix) + (yminus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime

                        jprime = xplus(ix) + (yminus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime
                      else
                        j = ix + (yplus(iy) - 1) * Nx
                        Hup(i,j) = -ty
                        Hdn(i,j) = -ty

                        j = ix + (yminus(iy) - 1) * Nx
                        Hup(i,j) = -ty
                        Hdn(i,j) = -ty

                        jprime = xplus(ix) + (yplus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime

                        jprime = xminus(ix) + (yplus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime

                        jprime = xminus(ix) + (yminus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime

                        jprime = xplus(ix) + (yminus(iy) - 1) * Nx
                        Hup(i,jprime) = -tprime
                        Hdn(i,jprime) = -tprime
                      endif
                  end if
              end do
          end do
      end if
!******************************************************************
!******************************************************************


!******************************************************************
!     Check the Hamiltoian setup
!******************************************************************
C       write(*, *) ' Initial Hamiltonian for spin-up electrons '
C       do i = 1, N
C           write(*, "(16f8.3)") (Hup(i,j), j=1,N)
C       end do

C       write(*, *) ' Initial Hamiltonian for spin-down electrons '
C       do i = 1, N
C           write(*, "(16f8.3)") (Hdn(i,j), j=1,N)
C       end do
!******************************************************************
!******************************************************************
      print *, 'The non-interacting spin-up matrix'
      do i = 1, 5
        print *, Hup(i, :)
      enddo   


!     INITIALIZE DENSITIES
      write (68, *) ' '
      write (68, *) ' Initial Electron Density '  

!*****************************************************************
!     Random initial distribution used for test
!*****************************************************************   
C        do i = 1, N
C            nup(i) = scale*ran2(iran)
C            ndn(i) = scale*ran2(iran)
C            write (68, "(i6, 2f12.6)") i,nup(i),ndn(i)
C        end do
!*****************************************************************
!*****************************************************************

!*****************************************************************
!     Antiferromagnetic order initialization
!*****************************************************************
      do i = 1, N
        tmp_i = mod(i, Nx)
        tmp_j = int((i-1)/Nx+1)
        if (mod(tmp_i+tmp_j, 2)==0) then
            nup(i) = 1.0d0*rho
            ndn(i) = 0.0d0*rho
        else
            nup(i) = 0.0d0*rho
            ndn(i) = 1.0d0*rho
        end if
        write (68, "(i6, 2f12.6)") i, nup(i), ndn(i)
        print '(i6, 2f12.6)', i, nup(i), ndn(i)
      enddo
!*****************************************************************
!*****************************************************************


C***********************************************************************
C     Read initial electron distribution from an input file
C***********************************************************************
C         open(unit=72,    status='old',    file="input.dat")
C         do i = 1, N
C           read(72, *) tmp, nup(i), ndn(i)
C C           write(*, *) tmp, nup(i), ndn(i)
C         end do
C         close(72)
C***********************************************************************
C***********************************************************************


C***********************************************************************
C     Read initial electron distribution from an input file
C***********************************************************************
C             beta_prev = beta - beta_step
C             beta_prev_int = int(beta_prev + 1.0d-8)
C             beta_prev_float = mod(nint(beta_prev * 10), 10)

C             if (beta_prev_int>=10 .and. abs(tprime_float)>=10) then
C               formatting_string = 
C      1        "(a, a, i1, a, i2, a, i1, a, i1, a, i2, a, i1, a)"
C             else if (beta_prev_int>=10 .and. abs(tprime_float)<10) then
C               formatting_string = 
C      1        "(a, a, i1, a, i1, a, i1, a, i1, a, i2, a, i1, a)"
C             else if (beta_prev_int<10 .and. abs(tprime_float)>=10) then
C               formatting_string = 
C      1        "(a, a, i1, a, i2, a, i1, a, i1, a, i1, a, i1, a)"
C             else
C               formatting_string = 
C      1        "(a, a, i1, a, i1, a, i1, a, i1, a, i1, a, i1, a)"
C             endif

C             write(input_name,   formatting_string) 
C      1      'Data/MF/PBC/tprime-0p2/rho1/p0/L16W4/Spin_', 
C      1      tprime_string, tprime_int, tmp_string, abs(tprime_float), 
C      1      '_U', Uint, 'p', Ufloat, '_beta', beta_prev_int, 'p', 
C      1      beta_prev_float, '.dat'
C             write(*, *) input_name

C             open(unit=72,    status='old',    file=input_name)
C             do i = 1, N
C               read(72, *) tmp, nup(i), ndn(i)
C               write(*, *) tmp, nup(i), ndn(i)
C             end do
C             close(72)
C***********************************************************************
C***********************************************************************


C C **********************************************************************
C C     Initialize Stripe Order Phase with Different Domain Size
C C **********************************************************************
C       if (index - number_of_trial_states < 1.0d-8) then
C           write(*, *) "************************************************"
C           write(*, *) "wavelength = ", index, 
C      1    "stripe order initialization"
C           write(*, *) if_pinning_field, pin_amplitude
C           write(*, *) "************************************************"
C             do i = 1, N
C                 tmp_i = mod(i, Nx)
C                 tmp_j = int((i - 1)/Nx) + 1
C                 phase_shift = int((i - 1)/index)
C                 if (mod(tmp_i+tmp_j+phase_shift, 2)==1) then
C                   nup(i) = 1.0d0*rho
C                   ndn(i) = 0.0d0
C                 else
C                   nup(i) = 0.0d0
C                   ndn(i) = 1.0d0*rho
C                 end if
C                 write(*, *) i, nup(i)-ndn(i)
C             end do
C       else
C             do i = 1, N
C                 tmp_i = mod(i, Nx)
C                 tmp_j = int((i-1)/Nx+1)
C                 if (mod(tmp_i+tmp_j, 2)==0) then
C                     nup(i) = 1.0d0*rho
C                     ndn(i) = 0.0d0*rho
C                 else
C                     nup(i) = 0.0d0*rho
C                     ndn(i) = 1.0d0*rho
C                 end if
C                 write (*, *) i, nup(i)-ndn(i)
C                 write (68, "(i6, 2f12.6)") i, nup(i), ndn(i)
C             end do
C       end if
C C **********************************************************************
C C **********************************************************************

C **********************************************************************
C **********************************************************************    
!     INCLUDE IN IHF TERMS IN H 
      do i = 1, N                      
           Hup(i,i) = Hup(i,i) + U * ndn(i)
           Hdn(i,i) = Hdn(i,i) + U * nup(i)
      end do
!***********************************************************************
!     Add pinning fields into the lattice on the left edge
!***********************************************************************
      do j = 1, Ny
          site = (j-1)*Nx + 1
          Hup(site,site) = Hup(site,site) 
     1     + (-1)**(j) * 0.5d0 * pin_amplitude
          Hdn(site,site) = Hdn(site,site) 
     1     + (-1)**(j+1) * 0.5d0 * pin_amplitude
      end do

      print *, 'The spin-up matrix'
      do i = 1, 5
        print *, Hup(i, :)
      enddo   


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

!     SELF CONSISTENT LOOP
      do it=1, Nit
C       write (6, "(' iteration= ', i6)") it

!     JACOBI/DSYEV DESTROY H, SO WORK WITH COPIES
      do i=1,N
          do j=1,N
              Hupcopy(i,j) = Hup(i,j)
              Hdncopy(i,j) = Hdn(i,j)
          end do
      end do

!      call jacobi(Hupcopy,N,N,evalup,evecup,nrot)
!      CALL SORT4(N,N,evalup,evecup)
!      call jacobi(Hdncopy,N,N,evaldn,evecdn,nrot)
!      CALL SORT4(N,N,evaldn,evecdn)

!        ARGUMENT 1:   N:  eigenvalues only,  V:eigenvectors
!        ARGUMENT 2:   U:  stored as upper tri,  L: stored as lower tri
!        INFO=0 SORTS EIGENVALUES.  LATER INFO TELLS OF PROBLEMS.
!        INPUT MATRIX TO BE DIAGONALIZED H IS OVERWRITTEN TO
!          CONTAIN EIGENVECTORS.

      INFO=0
      call DSYEV('V','U',N,Hupcopy,N,evalup,WORK,6*N,INFO)
      call DSYEV('V','U',N,Hdncopy,N,evaldn,WORK,6*N,INFO)
      do i=1,N
        do j=1,N
             evecup(i,j)=Hupcopy(i,j)
             evecdn(i,j)=Hdncopy(i,j)
        end do
      end do


!     COMPUTE ENERGY
      E = 0.0d0
      do i=1,N
          fermiup = 1.0d0 / (exp(beta*(evalup(i)-mu)) + 1.0d0)
          fermidn = 1.0d0 / (exp(beta*(evaldn(i)-mu)) + 1.0d0)
          E = E + fermiup*evalup(i) + fermidn*evaldn(i)
          E = E - U * nup(i) * ndn(i) 
      end do
      write (68, "('it, E/N =   ', i8, f16.8)") it, E/dfloat(N)
      totalE(it,index) = E/dfloat(N)

C **********************************************************************
C     Compute the grand free energy
C **********************************************************************
      Free_Energy = 0.0d0
      do i = 1, N
          Free_Energy=Free_Energy+log(1.0d0+exp(-beta*(evalup(i)-mu))) 
          Free_Energy=Free_Energy+log(1.0d0+exp(-beta*(evaldn(i)-mu)))
          Free_Energy=Free_Energy+beta*U*nup(i)*ndn(i) 
      end do 
      Free_Energy = -1.0d0 / beta * Free_Energy

      totalGF(it,index)=Free_Energy/dfloat(N)
      totalF(it,index)=Free_Energy/dfloat(N)+mu*rho
C       write(*, *) E/dfloat(N), Free_Energy/dfloat(N)+mu*rho, mu
C **********************************************************************
C **********************************************************************

!     RECOMPUTE DENSITIES FOR SELF-CONSISTENCY
      if (mod(20*it,Nit).eq.0) then
         write (68, "()")
         write (68, "(' old/new densities up and dn' )")
         write (68, "(' iteration =  ', i6)") it
      end if


      do i=1,N
          newnup(i) = 0.0d0
          newndn(i) = 0.0d0
          do j=1,N 
              fermiup = 1.0d0 / (exp(beta*(evalup(j)-mu)) + 1.0d0)
              newnup(i) = newnup(i) + fermiup*evecup(i,j)*evecup(i,j)
              fermidn = 1.0d0 / (exp(beta*(evaldn(j)-mu)) + 1.0d0)
              newndn(i) = newndn(i) + fermidn*evecdn(i,j)*evecdn(i,j)
          end do
          
          if (mod(50*it, Nit).eq.0) then
              write(68, "(i6, 2f12.6, 2f12.6)")  i, nup(i), newnup(i), 
     1         ndn(i), newndn(i)
          end if
      end do


      rho1=0d0
      do i = 1, N
        rho1 = rho1 + newnup(i) + newndn(i)
      end do
      rho1 = rho1 / dble(N)

!     RESET DIAGONALS OF H
!******************************************************************
!     Add the chemical potential tuning part
!******************************************************************
      mu_min = -10.0d0
      mu_max =  10.0d0
      rho_delta = 1.0d6 

      do while (rho_delta .gt. tolerance)
        mu_delta=(mu_min+mu_max)/2.0d0
        rho_tmp=0.0d0
        do i = 1, N
          tmpnup(i)=0.0d0
          tmpndn(i)=0.0d0
          do j = 1, N
            fermiup=1.0d0/(exp(beta*(evalup(j)-(mu+mu_delta)))+ 1.0d0)
            tmpnup(i) = tmpnup(i) + fermiup*evecup(i,j)*evecup(i,j)
            fermidn=1.0d0/(exp(beta*(evaldn(j)-(mu+mu_delta)))+ 1.0d0)
            tmpndn(i) = tmpndn(i) + fermidn*evecdn(i,j)*evecdn(i,j)
          enddo
        enddo

        do i = 1, N
          rho_tmp = rho_tmp + tmpnup(i) + tmpndn(i)
        end do
        rho_tmp = rho_tmp / dble(N)
        rho_delta = abs(rho_tmp - rho)

        if (rho_tmp .gt. rho) then 
          mu_max = mu_delta
        else if (rho_tmp < rho) then 
          mu_min = mu_delta
        end if
C         print *, rho_delta, tolerance
      end do

!***********************************************************************
!     Tune the electron density based on the new chemical potential
!***********************************************************************     
      mu=mu+mu_delta
      do i = 1, N
        newnup(i) = 0.d0
        newndn(i) = 0.d0
        do j = 1, N
          fermiup = 1.d0 / ( exp(beta*(evalup(j)-mu)) + 1.d0 )
          newnup(i) = newnup(i) + fermiup*evecup(i,j)*evecup(i,j)
          fermidn = 1.d0 / ( exp(beta*(evaldn(j)-mu)) + 1.d0 )
          newndn(i) = newndn(i) + fermidn*evecdn(i,j)*evecdn(i,j)
        end do
      end do

      rho2=0d0
      do i = 1, N
        rho2 = rho2 + newnup(i) + newndn(i)
      end do
      rho2 = rho2 / dble(N)
      write(71, "(i6, 3f12.6)") it, mu, rho1, rho2

C       write(*, *) mu
!***********************************************************************
!     Add Perturbation to electron densities
!***********************************************************************
      if ((mod(it, 50) .eq. 0) .and. (it .le. 2000)) then
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

      do i=1,N
           newnup(i) = relax*newnup(i) + (1.d0-relax)*nup(i)
           newndn(i) = relax*newndn(i) + (1.d0-relax)*ndn(i)
           Hup(i,i) = Hup(i,i) + U * ( newndn(i) - ndn(i) )
           Hdn(i,i) = Hdn(i,i) + U * ( newnup(i) - nup(i) )
           nup(i) = newnup(i)
           ndn(i) = newndn(i)
      end do

      end do
      
      do i = 1, N
        upElectron(i,index) = nup(i)
        dnElectron(i,index) = ndn(i)
      end do

C     Ending point for the different trial states
      end do

C***********************************************************************
C     Select and store the configuration that has the lowest
C     grand free energy.
C***********************************************************************
      if (number_of_trial_states-1<1E-8) then
        do i=1,N
          write(69, "(i6, 2f12.6)") i, upElectron(i,1), 
     1    dnElectron(i,1)
          write(82, "(i6, 2f12.6)") i, evalup(i), evaldn(i)
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
          write(69, "(i6, 2f12.6)") i, upElectron(i, optimal_index), 
     1    dnElectron(i, optimal_index)
        end do

        do i=1,Nit
          write(70, "(i6, f16.12)") i, totalE(i, optimal_index)
          write(80, "(i6, f16.12)") i, totalF(i, optimal_index)
          write(81, "(i6, f16.12)") i, totalGF(i, optimal_index)
        end do
      end if


C       do index=1,number_of_trial_states
C         do i=1,N
C           write(69, "(i6, 2f12.6)") i, upElectron(i,index),  
C      1    dnElectron(i,index)
C         end do
C       end do

C       do tmpInd = 1, number_of_trial_states
C         do i=1,Nit
C           write(70, "(i6, f16.12)") i, totalE(i, tmpInd)
C           write(80, "(i6, f16.12)") i, totalF(i, tmpInd)
C           write(81, "(i6, f16.12)") i, totalGF(i, tmpInd)
C         end do
C       end do

C     Ending point for the effective inverse temperature loop
      end do 

C     Ending point for the effective Coulomb interaction loop      
      end do

C     Ending point for the tprime loop
      end do

      end program IHF

C       subroutine reset_chemical_potential(tmp_min, tmp_max, value)
C         implicit none
C         real(kind=precision), intent(in)  :: value
C         real(kind=precision), intent(out) :: tmp_min, tmp_max

C         tmp_min = -value
C         tmp_max = value
C       end subroutine
