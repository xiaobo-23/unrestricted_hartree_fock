# !/usr/bin/env python
# 
# written by Bo Xiao <bxiao@flatironinsitute.org>
# 01/30/2023
#
#

'''
    Unrestricted Hartree-Fock method for 2D Hubbard Model
'''




import glob
import sys
import os
import numpy as np; # np.set_printoptions(threshold=sys.maxsize)
import math
import cmath
# import matplotlib.pyplot as plt
# import datetime; import random


#
# Set up parameters and hamiltonain
#


t=1.0                           # Nearest-neighbor hopping
tprime=-0.2                      # Next-nearest neighbor hopping
pinning=0.2                     # Edge pinning fields
Ueff=2.8                        # Effective Coulomb interaction
mu=0                            # Chemical potential to control particle number in a grand canonical ensemble
rho=0.875                       # Electron density
alpha=0.7                       # Factor to mix electron density from steps (n-2) and (n-1) at step n
beta=20                         # Inverse temperature beta=1/T
beta_list=[2,10]
Lx=4                            # Length of a cylinder 
Ly=4                            # Width of a cylinder
Niteration=500                  # Default number of iterations for the self-consistent convergence


# # NN_hopping=1.0; NNN_hopping=0.0
# Lx=16; Ly=4; N_sites=Lx*Ly; Nup=int(N_sites/2); Ndn=int(N_sites/2); N_electron=N_sites
# pinning=0.2; effective_potential=2.55; chemical_potential=0.
# rho=0.875; # chemical_upper_limit=2.; chemical_lower_limit=0.; global rho_diff
# delta=0; alpha=0.7; beta_list=[6] 
# # convergence_residual=1e5
# # criteria=10e-8 
# Niteration=500


# Set up the Hamiltonian for 
def generate_hamiltonian(t1, t2, lx, ly, pinning_field, effective_u, environment, spin_index, mu):
    num_sites=lx*ly
    matrix=np.zeros((num_sites, num_sites), dtype=np.float64)

    '''
        Need to be generalized. 
        lx >= 2 and ly >= 2 in this version.
    '''

    for i in range(lx):
        for j in range(ly):
            index=i+j*lx
            i_r=np.mod(i+1, lx)
            i_l=np.mod(i-1+lx, lx)
            j_u=np.mod(j-1+ly, ly)
            j_d=np.mod(j+1, ly)

            ## Get index for nearest neighbors
            index1=i_r+j*lx
            index2=i+j_u*lx
            index3=i_l+j*lx
            index4=i+j_d*lx

            ## Get index for next-nearest neighbors
            index5=i_r+j_u*lx
            index6=i_l+j_u*lx
            index7=i_l+j_d*lx
            index8=i_r+j_d*lx

            #
            # Periodic boundary 
            #
            matrix[index, index2]+=-t1
            matrix[index, index4]+=-t1

            ## Open boundary condition in x direction
            if i_r==i+1:
                matrix[index, index1]+=-t1
                matrix[index, index5]+=-t2
                matrix[index, index8]+=-t2

            if i_l==i-1:
                matrix[index, index3]+=-t1
                matrix[index, index6]+=-t2
                matrix[index, index7]+=-t2

            if i==0:
                matrix[index, index]+=0.5*pinning_field*pow(-1, i+j+spin_index)

            # matrix[index, index]+=effective_u*environment[index]-.5*effective_u
            matrix[index, index]+=effective_u*environment[index]
            Hermitian=matrix.conj().T

    try:
        np.max(np.abs(matrix-Hermitian)) <1e-8
    except ValueError:
        print("Oops!  The input tight-binding Hamiltonian is not Hermitian.  Please check carefully ...")
    # print("matrix difference", np.max(np.abs(matrix-Hermitian)))
    return matrix


# ## Initialize the density distributions for spin-up and spin-down electrons
# def init_electron_density(width, height, doping, density_up, density_dn):
#     # wavelength=int(2*int(1./doping))
#     wavelength=16
#     num_waves=int(width/wavelength)
#     # print(num_waves)
#     # print(wavelength)
#     nodal_points=np.zeros([2*num_waves])
#     amplitude=0.2

#     spin_up=np.zeros([width*height])
#     spin_dn=np.zeros([width*height])

#     for index in range(num_waves):
#         nodal_points[2*index]=int(index*wavelength+wavelength/4)
#         nodal_points[2*index+1]=int(index*wavelength+int(0.75*wavelength))
#     # print(nodal_points)

#     for tmp2 in range(height):
#         for tmp1 in range(width):
#             site_index=tmp1+tmp2*width
#             switch=False
#             for wave_index in range(num_waves):
#                 if(tmp1 >= nodal_points[2*wave_index] and tmp1 < nodal_points[2*wave_index+1]):
#                     switch=True
#             # print(switch)
#             if (switch==True):
#                 spin_up[site_index]=density_up*1./(width*height)-amplitude*pow(-1, tmp1+tmp2)
#                 spin_dn[site_index]=density_dn*1./(width*height)+amplitude*pow(-1, tmp1+tmp2)
#             else:
#                 spin_up[site_index]=density_up*1./(width*height)+amplitude*pow(-1, tmp1+tmp2)
#                 spin_dn[site_index]=density_dn*1./(width*height)-amplitude*pow(-1, tmp1+tmp2)
#             # print((-1)**(tmp1+tmp2)*(spin_up[site_index]-spin_dn[site_index]))
#     return spin_up, spin_dn

# ## Initialize the density distribution of spin-up and spin-down 
# ## electrons using a random distribution
# def init_electron_density_random(width, height, scale, doping):
#     spin_up=np.zeros([width*height])
#     spin_dn=np.zeros([width*height])

#     for i in range(width*height):
#         spin_up[i]=scale*random.uniform(0, 1)
#         spin_dn[i]=scale*random.uniform(0, 1)
#     return spin_up, spin_dn

#
# Initialize the spin-up and spin-down electrons using antiferromagnetic (AFM) order
#

def init_electron_density_antiferromagnetic(length, width, electron_density):
    spin_up=np.zeros([length * width])
    spin_dn=np.zeros([length * width])

    for j in range(width):
        for i in range(length):
            index=i+j*length
            if (i+j)%2==0:
                spin_up[index]=electron_density
                spin_dn[index]=0
            else:
                spin_up[index]=0
                spin_dn[index]=electron_density
    return spin_up, spin_dn



# ## Initialize the density distribution of spin-up and spin-down electrons
# ## by reading the input file
# def init_electron_input(filename, width, height):
#     '''Initialize the electron distribution from the input file'''
#     spin_up=np.zeros([width*height])
#     spin_dn=np.zeros([width*height])

#     input_matrix=np.loadtxt(filename)
#     if(input_matrix.shape[0]==width*height):
#         spin_up[:]=input_matrix[:, 1]
#         spin_dn[:]=input_matrix[:, 2]
#         return spin_up, spin_dn
#     else:
#         print("*********************************************************************************")
#         print("Error Message: ")
#         print("The dimensionality of the input spin densities doesn't match the parameters!!")
#         print("*********************************************************************************")
#         sys.exit()

#     # for index in range(input_matrix.shape[0]):
#     #     spin_up[index]=input_matrix[index, 1]
#     #     spin_dn[index]=input_matrix[index, 2]
#     # return spin_up, spin_dn



# Calculate the Fermi-Dirac function
def Fermi_Dirac(inverse_temp, energy, chem_potential):
    return 1./(1.+np.exp(inverse_temp*(energy-chem_potential)))


## Define the Fermi_Dirac function for the grand canonical ensemble
# def Fermi_Dirac(energy, inverse_temp):
#     # energy=np.array(energy, dtype=np.float128)
#     return 1./(1.+np.exp(inverse_temp*energy))


## Compute the electron density matrix
def solve_Hamiltonian(tmp_matrix, inverseT, chem):
    value_tmp, vector_tmp=np.linalg.eigh(tmp_matrix, UPLO='L')
    density_tmp=np.zeros([vector_tmp.shape[0], vector_tmp.shape[1]])
    
    for i in range(vector_tmp.shape[0]):
        for  j in range(vector_tmp.shape[1]):
            for eigen_index, eigen in enumerate(value_tmp):
                density_tmp[i,j]+=vector_tmp[i, eigen_index]*vector_tmp[j, eigen_index]\
                    *Fermi_Dirac(inverseT, eigen, chem)
    return value_tmp, density_tmp


## Set up and solve the Hamiltonian
def set_up_and_solve_Hamiltonian(t1, t2, length, width, pinning_field, Hubbard_U, \
    input_density, spin_index, inverseT, chem):
    tmp_matrix=generate_hamiltonian(t1, t2, length, width, pinning_field, \
        Hubbard_U, input_density, spin_index, chem)
    print(tmp_matrix)
    value_tmp, vector_tmp=np.linalg.eigh(tmp_matrix, UPLO='L')
    density_tmp=np.zeros([vector_tmp.shape[0], vector_tmp.shape[1]])

    for i in range(vector_tmp.shape[0]):
        for  j in range(vector_tmp.shape[1]):
            for eigen_index, eigen in enumerate(value_tmp):
                density_tmp[i,j]+=vector_tmp[i, eigen_index]*vector_tmp[j, eigen_index]\
                    *Fermi_Dirac(inverseT, eigen, chem)
    return value_tmp, vector_tmp, density_tmp


## Compute the L2 or L1 norm based on the choice
def compute_norm(array1, array2, length):
    if array1.shape[0]==length and array2.shape[0]==length:
        # print(array1.shape)
        # print("Yeah!")
        result=0
        for tmp_index in range(length):
            result+=math.pow(array1[tmp_index]-array2[tmp_index], 2)
        result=np.sqrt(result)
        result/=(1.*length)
    return result


## Compute the energy based on eigenvalues and eigenvectors
def compute_energy(eigen_up, eigen_dn, inv_temp, chemical):
    total_energy=0.
    for index in range(eigen_up.shape[0]):
        fermiup=1./(1.+np.exp(inv_temp*(eigen_up[index]-chemical)))
        fermidn=1./(1.+np.exp(inv_temp*(eigen_dn[index]-chemical)))
        total_energy+=fermiup*eigen_up[index]+fermidn*eigen_dn[index]
    return total_energy


## Tune chemical self-consistently to find the corresponding density
# def tune_chemical(t1, t2, width, height, pinning_field, Hubbard_U, up, dn, inv_temp, filling):
#     filling_diff=1e6
#     upper_limit=2.
#     lower_limit=-2.
#     chem_tmp=0

#     while(abs(filling_diff)>1e-8):
#         chem_tmp=(upper_limit+lower_limit)/2.
#         matrix_up=hamiltonian(t1, t2, width, height, \
#             pinning_field, Hubbard_U, dn, 1, chem_tmp)
#         matrix_dn=hamiltonian(t1, t2, width, height, \
#             pinning_field, Hubbard_U, up, 0, chem_tmp)
#         evalup, up_density=solve_Hamiltonian(matrix_up, inv_temp, chemical_potential)
#         evaldn, dn_density=solve_Hamiltonian(matrix_dn, inv_temp, chemical_potential)

#         rho_total=(up_density.diagonal().sum()+dn_density.diagonal().sum())/N_sites
#         filling_diff=rho_total-filling

#         if rho_total>filling:
#             lower_limit=chem_tmp
#         elif rho_total<filling:
#             upper_limit=chem_tmp
#         # print(rho_total)
#         # print(chem_tmp)
#         # print("\n")

#     return chem_tmp


def tune_chemical(eval_up, evector_up, eval_dn, evector_dn, inverseT, filling):
    filling_diff=float('inf')
    upper_limit=2.
    lower_limit=-5.
    chem_tmp=0

    while(abs(filling_diff)>1e-8):
        chem_tmp=(upper_limit+lower_limit)/2.
        up_tmp=np.zeros([evector_up.shape[0], evector_up.shape[1]])
        dn_tmp=np.zeros([evector_dn.shape[0], evector_dn.shape[1]])

        for i in range(evector_up.shape[0]):
            for  j in range(evector_up.shape[1]):
                for eigen_index_up, eigen_up in enumerate(eval_up):
                    up_tmp[i, j]+=evector_up[i, eigen_index_up]*evector_up[j, eigen_index_up]\
                    *Fermi_Dirac(inverseT, eigen_up, chem_tmp)

        for k in range(evector_dn.shape[0]):
            for  l in range(evector_dn.shape[1]):
                for eigen_index_dn, eigen_dn in enumerate(eval_dn):
                    dn_tmp[k, l]+=evector_dn[k, eigen_index_dn]*evector_dn[l, eigen_index_dn]\
                    *Fermi_Dirac(inverseT, eigen_dn, chem_tmp)

        rho_total=(up_tmp.diagonal().sum()+dn_tmp.diagonal().sum())/N_sites
        # print(rho_total)
        filling_diff=rho_total-filling

        if rho_total>filling:
            upper_limit=chem_tmp
        elif rho_total<filling:
            lower_limit=chem_tmp
        # print(chem_tmp)
        # print(rho_total)
        # print(chem_tmp)
        # print("\n")

    return chem_tmp, up_tmp.diagonal(), dn_tmp.diagonal()


## Update electron density whenever there is a change
# def update_electron_density(t1, t2, width, height, pinning_field, Hubbard_U, \
#     spin_tmp, spin_index, inv_temp, chem):
#     matrix_tmp=hamiltonian(t1, t2, width, height, \
#             pinning_field, Hubbard_U, spin_tmp, spin_index, chem)
#     eval_tmp, up_tmp=solve_Hamiltonian(matrix_tmp, inv_temp)
#     return eval_tmp, up_tmp



## The main function to perform the calculation
def main():
    residual_np=np.zeros([len(beta_list), Niteration])
    E_list=np.zeros([len(beta_list), Niteration])
    spin_np=np.zeros([len(beta_list), Lx])

    for beta in beta_list:
        # up_electron, dn_electron=init_electron_density(Lx, Ly, delta, Nup, Ndn)
        # up_electron, dn_electron=init_electron_density_random(Lx, Ly, 1., delta)
        # print(up_electron)
        # print(dn_electron)

        #####################################################################
        ## Perform 1st round self-consistent procedure
        #####################################################################
        #up_electron, dn_electron=init_electron_density_random(Lx, Ly, 1., delta)
        up_electron, dn_electron=init_electron_density_antiferromagnetic(Lx, Ly, rho)
        # input_file="./Initial_Electron_Density.txt"
        # up_electron, dn_electron=init_electron_input(input_file, Lx, Ly)
        print("*********************************************************************")
        print("Initialize the electron density distribution")
        print("*********************************************************************")
        print(up_electron)
        print(dn_electron)
        print("\n")
        evalup, evectup, density_up=set_up_and_solve_Hamiltonian(t, tprime, Lx, Ly, \
            pinning, Ueff, dn_electron, 1, beta, mu)
        evaldn, evectdn, density_dn=set_up_and_solve_Hamiltonian(t, tprime, Lx, Ly, \
            pinning, Ueff, up_electron, 0, beta, mu)

        # print(density_up.diagonal().sum()/N_sites)
        # print(density_dn.diagonal().sum()/N_sites)
        # print("\n")
        # print("\n")

        # chem1, spin_up_density, spin_dn_density=tune_chemical(evalup, evectup, evaldn, evectdn, beta, rho)
        # print(spin_up_density.sum()/N_sites)
        # print(spin_dn_density.sum()/N_sites)
        # print("\n")
        # print("\n")


        # up_electron=np.vstack((up_electron, spin_up_density))
        # dn_electron=np.vstack((dn_electron, spin_dn_density))
        # print(up_electron[up_electron.shape[0]-1, :])
        # print(dn_electron[dn_electron.shape[0]-1, :])
        

        # ####################################################################
        # ## 2nd round self-consistent procedure
        # ####################################################################
        # evalup, evectup, density_up=set_up_and_solve_Hamiltonian(NN_hopping, NNN_hopping, Lx, Ly, \
        #     pinning, effective_potential, spin_dn_density, 1, beta, chem1)
        # evaldn, evectdn, density_dn=set_up_and_solve_Hamiltonian(NN_hopping, NNN_hopping, Lx, Ly, \
        #     pinning, effective_potential, spin_up_density, 0, beta, chem1)

        # print(density_up.diagonal().sum()/N_sites)
        # print(density_dn.diagonal().sum()/N_sites)
        # print("\n")
        # print("\n")

        # chem2, spin_up_density, spin_dn_density=tune_chemical(evalup, evectup, evaldn, evectdn, beta, rho)

        # print(spin_up_density.sum()/N_sites)
        # print(spin_dn_density.sum()/N_sites)
        # print("\n")
        # print("\n")


        # up_electron=np.vstack((up_electron, spin_up_density))
        # dn_electron=np.vstack((dn_electron, spin_dn_density))
        # print(up_electron[up_electron.shape[0]-1, :])
        # print(dn_electron[dn_electron.shape[0]-1, :])


        # ###############################################################################
        # ## Perform iterations self-consistently
        # ###############################################################################
        # residual_list=[]
        # energy_list=[]
        # mu=chem2

        # for iteration_index in range(int(Niteration)):
        #     residual=0
        #     E=0
        #     EP=0

        #     ## Update spin-up electron density distribution due to changes in spin-down mean-field
        #     dn_count=dn_electron.shape[0];  print(dn_count)
        #     mix_dn=alpha*dn_electron[dn_count-1, :]+(1-alpha)*dn_electron[dn_count-2, :]
        #     tmp_evalup, tmp_evectup, tmp_density_up=set_up_and_solve_Hamiltonian(NN_hopping, NNN_hopping, \
        #         Lx, Ly, pinning, effective_potential, mix_dn, 1, beta, mu)


        #     ## Update spin-down electron distribution due to changes in spin-up electron mean-field
        #     up_count=up_electron.shape[0];  print(up_count)
        #     mix_up=alpha*up_electron[dn_count-1, :]+(1-alpha)*up_electron[dn_count-2, :]
        #     tmp_evaldn, tmp_evectdn, tmp_density_dn=set_up_and_solve_Hamiltonian(NN_hopping, NNN_hopping, \
        #         Lx, Ly, pinning, effective_potential, mix_up, 0, beta, mu)

        #     E=compute_energy(tmp_evalup, tmp_evaldn, beta, mu)
        #     for index_tmp in range(tmp_density_up.shape[0]):
        #         # print(tmp_density_up.shape)
        #         EP+=effective_potential*tmp_density_up.diagonal()[index_tmp]\
        #         *tmp_density_dn.diagonal()[index_tmp]
        #     print(EP)
        #     print("\n")
        #     print("\n")
        #     E-=EP
        #     energy_list.append(E)
        #     # print(energy_list)
        #     # print("\n")
        #     # print("\n")
            
        #     chem3, spin_up_density, spin_dn_density=tune_chemical(tmp_evalup, tmp_evectup, \
        #         tmp_evaldn, tmp_evectdn, beta, rho)

        #     print(spin_up_density.sum()/N_sites)
        #     print(spin_dn_density.sum()/N_sites)
        #     print(chem3)
        #     print("\n")
        #     print("\n")
            
        #     # potential=0
        #     # for index in range(spin_up_density.shape[0]):
        #     #     potential+=effective_potential*spin_up_density[index]*spin_dn_density[index]
        #     # print(potential)
        #     # print("\n")
        #     # print("\n")

        #     ## Reset chemical potential and save the convergence history
        #     mu=chem3

        #     up_electron=np.vstack((up_electron, spin_up_density))
        #     dn_electron=np.vstack((dn_electron, spin_dn_density))

        #     residual+=compute_norm(up_electron[up_electron.shape[0]-1, :], \
        #         up_electron[up_electron.shape[0]-2, :], N_sites)
        #     residual+=compute_norm(dn_electron[dn_electron.shape[0]-1, :], \
        #         dn_electron[dn_electron.shape[0]-2, :], N_sites)
        #     residual_list.append(residual)

        #     # E=compute_energy(tmp_evalup, tmp_evaldn, beta, mu)
        #     # energy_list.append(E)
            
        #     # print(residual_list)
        
        # # Store output results at different inverse temperature  
        # spin_plot=np.zeros([Lx])
        # for index1 in range(Lx):
        #     # spin_plot[index1]=(-1)**(index1)*(up_electron[up_electron.shape[0]-1, index1]\
        #     #     -dn_electron[dn_electron.shape[0]-1, index1])/2.
        #     spin_plot[index1]=(up_electron[up_electron.shape[0]-1, index1]\
        #         -dn_electron[dn_electron.shape[0]-1, index1])/2.
        
        # residual_np[beta_list.index(beta), :]=residual_list
        # E_list[beta_list.index(beta), :]=energy_list
        # spin_np[beta_list.index(beta), :]=spin_plot

        # output_file=open("./spin_L16_beta6.txt", "a")
        # for index1 in range(Ly):
        #     for index2 in range(Lx):
        #         index=index2+index1*Lx
        #         output_file.write("{:>6d}".format(index)+"{:3s}".format(" ")\
        #         +"{:>20.16e}".format(up_electron[up_electron.shape[0]-1, index])+"{:3s}".format(" ")\
        #         +"{:>20.16e}".format(dn_electron[dn_electron.shape[0]-1, index])+"\n")
        # output_file.close()

        # output_file=open("./energy_L16_beta6.txt", "a")
        # for index in range(E_list.shape[1]):
        #     output_file.write("{:>6d}".format(index)+"{:3s}".format(" ")\
        #         +"{:>20.16e}".format(E_list[beta_list.index(beta), index])+"\n")
        # output_file.close()

        # output_file=open("./residual_L16_beta6.txt", "a")
        # for index in range(residual_np.shape[1]):
        #     output_file.write("{:>6d}".format(index)+"{:3s}".format(" ")\
        #         +"{:>20.16e}".format(residual_np[beta_list.index(beta), index])+"\n")
        # output_file.close()

main()