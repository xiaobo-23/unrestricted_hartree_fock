    ! 06/29/2023
    ! Simplify the Hartree-Fock code to benchmark agasint Stewart's code
    module compute_energy_and_density
        implicit none
        contains

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
           print "('iteration, E/N = ', i6, f16.8)", tmp_iteration, tmp_energy / dfloat(tmp_Nsites)
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

    end module compute_energy_and_density