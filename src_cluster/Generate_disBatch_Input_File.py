##################################################################
## Generate input file for disBatch
## Run holoQUADS circuits in parallel
##################################################################

import numpy as np
import math

def generate_input_file(input_index, task_file):
    '''Generate corresponding folders and input files based on chemical potential'''
    
    fractional_part, _ = math.modf(input_index)
    
    if fractional_part > 1e-8:
        folder_name = "beta" + "{}".format(input_index) + "/"
    else:
        folder_name = "beta" + "{}".format(int(input_index)) + "/"
    task_file.write("cd " + folder_name \
        + " &&  UHF.o" + " &> UHF_" \
        + "{}".format(input_index) + ".log" + "\n")
    
def main():
    lower_bound=2.0
    upper_bound=50.1
    sample_list = np.round(np.arange(lower_bound, upper_bound, 0.1), 1)

    submit_file = open("Submission_U2.1", "a")
    for tmp in sample_list:
        print(tmp)
        generate_input_file(tmp, submit_file)
    submit_file.close()    

main()