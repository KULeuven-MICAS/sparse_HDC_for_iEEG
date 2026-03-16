import scipy.io as scp
from pathlib import Path
import numpy as np
import math
import random
import time
import HDC_functions as HDC


for iteration in [5]:
    first = 1
    generate_new_HVs_ = 0 #CHANGE!!!!

    shift_bind = 0
    if shift_bind:
        bind_str = 'shift binding'
    else:
        bind_str = 'segm shift binding'

    thin_first_bund = 1
    if thin_first_bund: #SWAPPED IN DATA
        thin_str = 'both thinned'
    else:
        thin_str = 'only temp thinned'

    for patient in [13,2,4,5,6,8,11,16]: #[13,2,4,5,6,8,11,16]
        D = 1024 #HV length
        LBP_length = 6 #number of bits
        p = 0.0078125 #density
        nb_segments = int(round(D*p))
        length_segment = int(round(1/p))

        class_bundling_p = 0.50
        if generate_new_HVs_ and first:
            generate_new_HVs = 1
            first = 0
        else: 
            generate_new_HVs = 0


        dict_pat_to_nb_seizures = {2:4, 4:4, 5:6, 6:2, 8:3, 11:2, 13:2, 16:2}
        dict_pat_to_channels = {2:64, 4:42, 5:59, 6:36, 8:61, 11:59, 13:98, 16:64}
        channels = dict_pat_to_channels[patient]
        nb_LBP = dict_pat_to_nb_seizures[patient]

        #generate_new_HVs = 0
        #thinning = 0 #0 means majority thinning, 1 is laiho thinning, 2 is random thinning, 3 is CDT thinning

        #mode = 2 #0 is dense, 1 is sparse with segmented shifting as binding, 2 is sparse with segmented shifting with thinning on both bundlings 
        #3 is sparse with binding by shifting whole EM vector by LBP value, 
        #4 is sparse with binding by shifting each segment by the LBP value,
        #5 is sparse with CDT-binding
        #Only mode 1 is supported with thinning 0!!!


        working_dir = Path(__file__).parent / 'no_backup' / ('Pat'+str(patient))


        # K = 32
        # random_perm_list = HDC.random_permutation_list(K)


        #IM and EM initialization
        if (generate_new_HVs):
            saving_dir = Path(__file__).parent / 'no_backup' / str(iteration)
            saving_dir.mkdir(exist_ok=True, parents=True)
            IM = HDC.compressed_IM_sparse(2**(LBP_length), nb_segments, length_segment)
            EM = HDC.EM_sparse(98, nb_segments, length_segment)
            #write to files:
            with open(saving_dir / 'IM.txt','w') as file: 
                for index in IM:
                    HV = IM[index]
                    file.writelines('[')
                    for hv_index in range(len(HV)):
                        file.writelines(str(HV[hv_index]) + ' ')
                    file.writelines(']'+'\n')

            with open(saving_dir / 'EM.txt','w') as file:
                for index in EM:
                    HV = EM[index]
                    for i in range(len(HV)):
                        file.writelines('[')
                        segment = HV[i]
                        for j in range(len(segment)):
                            file.writelines(str(segment[j]) + ' ')
                        file.writelines(']'+'\n')


        else: #read from files
            IM = {}
            vector = []
            file = open(Path(__file__).parent / 'no_backup' / str(iteration) / 'IM.txt','r')
            content = file.read()
            i = 1
            IM_counter = 0
            content = content.split("\n")
            content = content[0:-1]
            for element in content:
                temp = element[1:len(element)-2]
                temp_list = np.array(temp.split(" "))
                IM[IM_counter] = temp_list.astype(int)
                IM_counter += 1
            file.close()

            EM = {}
            vector = []
            file = open(Path(__file__).parent / 'no_backup' / str(iteration) / 'EM.txt','r')
            content = file.read()
            i = 1
            EM_counter = 0
            content = content.split("\n")
            for n in range(98):
                part_content = content[n*nb_segments:(n+1)*nb_segments]
                HV = []
                for element in part_content:
                    temp = element[1:len(element)-2]
                    temp_list = np.array(temp.split(" "))
                    # print(temp_list)
                    HV.append(temp_list.astype(int))
                EM[n] = np.array(HV)
                EM_counter += 1
            file.close()
            
        use_p_max_or_threshold = 0 #0 is p_max, 1 is threshold
        # p_max_or_threshold_list = [0.02,0.04,0.06,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.38,0.40,0.42,0.44,0.46,0.48,0.50,0.52,0.54,0.56,0.58,0.60,0.62,0.64,0.66,0.68,0.70,0.72,0.74,0.76,0.78,0.80,0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98]
        # p_max_or_threshold_list = [*range(256)]
        # p_max_or_threshold_list = [120,130,180]

        # p_max_or_threshold_list = [10,50,85]
        p_max_or_threshold_list = [*range(1,100)]
        for i in range(len(p_max_or_threshold_list)):
            p_max_or_threshold_list[i] = p_max_or_threshold_list[i]/100
        


        #Max number of halfseconds of a seizure for this patient:
        len_list = []
        for i in range(1, nb_LBP+1):
            LBP_address = working_dir / ('LBP_' + str(i) + '.mat')
            f1 = scp.loadmat(LBP_address)
            LBP_train = np.array(f1['LBP'])
            len_list.append(len(LBP_train))


        #First train for each LBP
        ictal_HV_list_list = np.zeros((len(p_max_or_threshold_list),nb_LBP,D), dtype=int)
        interictal_HV_list_list = np.zeros((len(p_max_or_threshold_list),nb_LBP,D), dtype=int)
        block_outputs_per_halfsecond_list_list = np.zeros((len(p_max_or_threshold_list),nb_LBP,math.ceil((max(len_list)/256)-360),D), dtype=int) #not doing last 180 seconds as is postictal state
        nb_of_halfseconds_list = []
        print(f'Patient {patient}:',flush=True)
        for i in range(1, nb_LBP+1):
            print('Starting LBP' + str(i) + ':',flush=True)
            LBP_address = working_dir / ('LBP_' + str(i) + '.mat')
            f1 = scp.loadmat(LBP_address)
            LBP_train = np.array(f1['LBP'])

            nb_of_halfseconds = math.ceil(len(LBP_train)*1.0/256-360)
            nb_of_halfseconds_list.append(nb_of_halfseconds)
            for k in range(nb_of_halfseconds):
                if (k % 5 == 0):
                    print('Progress: '+str(round(k/nb_of_halfseconds*100,2))+'%',flush=True)
                LBP_part = LBP_train[k*256:((k+1)*256),:]
                block_outer_ = [ [] for _ in range(len(p_max_or_threshold_list)) ]
                block_outer = []
                for k2 in range(len(LBP_part)):
                    LBP_row = LBP_part[k2]
                    block_inner = []
                    for electrode_nb in range(len(LBP_row)):
                        if shift_bind:
                            block_inner.append((HDC.perm((EM[electrode_nb]).flatten(), -int(LBP_row[electrode_nb])))) #shfit bind
                        else:
                            block_inner.append(HDC.binding_sparse_segm_shift_fast(-IM[LBP_row[electrode_nb]], EM[electrode_nb], D)) #segm shift bind
                    if thin_first_bund:
                        if (use_p_max_or_threshold == 1):
                            ValueError("Should use p_max!")
                        for j in range(len(p_max_or_threshold_list)):
                            block_outer_[j].append(HDC.bundle_sparse_time_ideal(block_inner, p_max_or_threshold_list[j], D))
                            # print('Density_inner:')
                            # print(p_max_or_threshold_list[j])
                            # print(HDC.similarity_sparse_fast(block_outer_[j][-1],[1]*D,D)/D)
                            # time.sleep(1)
                    else:
                        block_outer.append(HDC.bundle_sparse_space(block_inner))
                if (use_p_max_or_threshold == 0):
                    for j in range(len(p_max_or_threshold_list)):
                        if (thin_first_bund):
                            block_outputs_per_halfsecond_list_list[j][i-1][k] = HDC.bundle_sparse_time_ideal(block_outer_[j], p_max_or_threshold_list[j], D)
                            # print('Density_outer:')
                            # print(p_max_or_threshold_list[j])
                            # print(HDC.similarity_sparse_fast(block_outputs_per_halfsecond_list_list[j][i-1][k],[1]*D,D)/D)
                            # time.sleep(10)
                        else:
                            block_outputs_per_halfsecond_list_list[j][i-1][k] = HDC.bundle_sparse_time_ideal(block_outer, p_max_or_threshold_list[j], D)

                elif (use_p_max_or_threshold == 1):
                    for j in range(len(p_max_or_threshold_list)):
                        block_outputs_per_halfsecond_list_list[j][i-1][k] = HDC.bundle_sparse_time(block_outer, p_max_or_threshold_list[j], D)
                
            for j in range(len(p_max_or_threshold_list)):
                interictal_HV_list_list[j][i-1] = HDC.bundle_sparse_time_ideal(block_outputs_per_halfsecond_list_list[j][i-1][:360], class_bundling_p, D)
                ictal_HV_list_list[j][i-1] = HDC.bundle_sparse_time_ideal(block_outputs_per_halfsecond_list_list[j][i-1][360:], class_bundling_p, D)


            #write ictal and interictal HV
            for j in range(len(p_max_or_threshold_list)):
                if (use_p_max_or_threshold):
                    writing_dir = Path(__file__).parent / 'no_backup' / str(iteration) / ('Pat'+str(patient)) / (bind_str +' '+thin_str) / ('threshold '+str(p_max_or_threshold_list[j]))
                else:
                    writing_dir = Path(__file__).parent / 'no_backup' / str(iteration) / ('Pat'+str(patient)) / (bind_str +' '+thin_str) / ('p_max '+str(p_max_or_threshold_list[j]))

                writing_dir.mkdir(exist_ok=True, parents=True)
                folder_to_write = writing_dir / ('LBP_' + str(i))
                folder_to_write.mkdir(exist_ok=True)
                file = open(folder_to_write / 'ictal_HV.txt', 'w')
                file.writelines('[')
                for hv_index in range(D):
                    file.writelines(str(ictal_HV_list_list[j][i-1][hv_index]) + ' ')
                file.writelines(']'+'\n')
                file.close()
                file = open(folder_to_write / 'interictal_HV.txt', 'w')
                file.writelines('[')
                for hv_index in range(D):
                    file.writelines(str(interictal_HV_list_list[j][i-1][hv_index]) + ' ')
                file.writelines(']'+'\n')
                file.close()



        #Now test for every combination of LBPs
        print('Similarity Search in progress:')
        for i in range(1, nb_LBP+1):
            print('Starting LBP' + str(i) + ':')
            for j in range(len(p_max_or_threshold_list)):
                if (use_p_max_or_threshold):
                    writing_dir = Path(__file__).parent / 'no_backup' / str(iteration) / ('Pat'+str(patient)) / (bind_str +' '+thin_str) / ('threshold '+str(p_max_or_threshold_list[j]))
                else:
                    writing_dir = Path(__file__).parent / 'no_backup' / str(iteration) / ('Pat'+str(patient)) / (bind_str +' '+thin_str) / ('p_max '+str(p_max_or_threshold_list[j]))

                ictal_HV = ictal_HV_list_list[j][i-1]
                interictal_HV = interictal_HV_list_list[j][i-1]
                for i2 in range(1, nb_LBP+1):
                    sim_scores_interictal = []
                    sim_scores_ictal = []
                    classifications = [] #per halfsecond
                    for k in range(nb_of_halfseconds_list[i2-1]):
                        hv = block_outputs_per_halfsecond_list_list[j][i2-1][k]
                        sim_ictal = HDC.similarity_sparse_fast(hv,ictal_HV,D)
                        sim_scores_ictal.append(sim_ictal)
                        sim_interictal = HDC.similarity_sparse_fast(hv,interictal_HV,D)
                        sim_scores_interictal.append(sim_interictal)
                        if (sim_ictal > sim_interictal):
                            classifications.append(1)
                        else:
                            classifications.append(0)

                    #write classifications and sim_scores to file
                    folder_to_write_pre = writing_dir / ('LBP_' + str(i))
                    folder_to_write_pre.mkdir(exist_ok=True)
                    folder_to_write = folder_to_write_pre / ('LBP_' + str(i2))
                    folder_to_write.mkdir(exist_ok=True)
                    file = open(folder_to_write / 'classifications.txt', 'w')
                    for index in range(len(classifications)):
                        file.writelines(str(classifications[index]) + ', ')
                    file.close()
                    file = open(folder_to_write / 'sim_scores.txt', 'w')
                    for index in range(len(sim_scores_ictal)):
                        file.writelines('('+str(sim_scores_ictal[index])+','+str(sim_scores_interictal[index])+') ')
                    file.close()
