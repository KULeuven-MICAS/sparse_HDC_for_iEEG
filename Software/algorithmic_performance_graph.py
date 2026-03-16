import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

def str_to_list(string):
    string = string[0:len(string)-2]
    string = string.replace(",", "")
    lst = np.array(string.split(" ")).astype(int)
    return lst


#Do postprocessing on data as normal, get detection delay for each combination.
#Then instead of not reporting if not all detected, 
#now keep normalized detection rate and average detection delay of the ones that did detect.
def Postprocessing_better(patients, p_max_or_thr_list, p_max_or_thr_q, group, subpath):
    output_data_patient = []
    for patient in patients:
        dir = Path(__file__).parent.parent / 'no_backup' / subpath / f'Pat{patient}' / group
        output_data_patient.append({})
        for p_max_or_thr in p_max_or_thr_list:
            #find base directory
            if p_max_or_thr_q:
                #p_max
                directory = dir / ('p_max '+str(p_max_or_thr))
            else:
                #thr
                directory = dir / ('threshold '+str(p_max_or_thr))

            nb_seizures_dict = {2:4, 4:4, 5:6, 6:2, 8:3, 11:2, 13:2, 16:2} #patient linked to number of seizures
            nb_seizures = nb_seizures_dict[patient]
            width_window = 10 #10 halfseconds or 5 second window

            thresholds = [-9.5/10,-9.5/10,-7.5/10,-9.5/10,-9.5/10,-8.5/10,-9.5/10,-7.5/10,-8.5/10,-8.5/10,-9.5/10,-6.5/10,-8.5/10,-5.5/10,-9.5/10,-8.5/10] #same as original dense HDC code
            #for each patient
            threshold = thresholds[patient-1]

            detection_delays = []
            mispredicts = []
            for i in range(1, nb_seizures+1):
                for j in range(1, nb_seizures+1):
                    #read classification data
                    address = directory / ('LBP_' + str(i)) / ('LBP_' + str(j))
                    file = open(address / 'classifications.txt', 'r')
                    content = file.read()
                    file.close()

                    #process halfseconds results in real classification by postprocessing with sliding window
                    classification_post_proc_data = []
                    data = str_to_list(content)
                    first = 1
                    mispredict = 0
                    for k in np.arange(width_window, len(data)):
                        classification = ((np.mean(data[k-width_window:k])+threshold) >= 0)
                        classification_post_proc_data.append(classification)
                        if (classification) and (k < 360):
                            mispredict = 1
                            detection_delay = 0
                        elif (classification) and (first):
                            first = 0
                            detection_delay = (k-360)/2 #in seconds
                    if (first):
                        mispredict = 1
                        detection_delay = 0

                    #keep result for each test:
                    with open(address / 'class_post_proc.txt', 'w') as f:
                        f.write(str(classification_post_proc_data))
                        f.write('\n'+'Detection delay: '+str(detection_delay)+'s')
                        f.write('\n'+'Mispredict? '+str(mispredict))
                    if (mispredict==0):
                        detection_delays.append(detection_delay)
                    mispredicts.append(mispredict)
            
            #Average data over seizures of patient:
            with open(directory / 'result.txt', 'w') as f:
                prediction_rate = 1-np.mean(mispredicts)
                f.write(f'Prediction rate: {prediction_rate}\n')
                if (len(detection_delays)==0):
                    print(f'Warning: no right predictions for patient {patient}, p_max/thr {p_max_or_thr}!')
                else:
                    mean_delay = np.mean(detection_delays)
                    f.write(f'Mean detection delay: {mean_delay}')
                    #output data:
                    output_data_patient[-1][p_max_or_thr] = [prediction_rate, mean_delay]
    return output_data_patient


#Make graph that averages over the patients and displays average detection rate and delay over patients:
patients = [2,4,5,6,8,11,13,16]
p_max_or_thr_list_ = [*range(1,75)]
p_max_or_thr_list = [*range(1,75)]
for i in range(len(p_max_or_thr_list_)):
    p_max_or_thr_list[i] = p_max_or_thr_list_[i]/100

which = ['segm shift binding both thinned', 'segm shift binding only temp thinned', 'shift binding only temp thinned']

# 1. Dynamically find which subpaths actually exist to avoid hard crashes
base_dir = Path(__file__).parent.parent / 'no_backup'
possible_runs = ['1', '2', '3', '4']
subpaths_dict = {}

for group in which:
    valid_subpaths = []
    for subpath in possible_runs:
        # Check if the directory for the first patient exists to confirm the run is available
        check_dir = base_dir / subpath / f'Pat{patients[0]}' / group
        if check_dir.exists():
            valid_subpaths.append(subpath)
    subpaths_dict[group] = valid_subpaths

# Get data for every subpath individually
output_data_all = [] 
for group in which:
    group_subpath_data = []
    for subpath in subpaths_dict[group]:
        output_data = Postprocessing_better(patients, p_max_or_thr_list, 1, group=group, subpath=subpath)
        group_subpath_data.append(output_data)
    output_data_all.append(group_subpath_data)

# 2. Average the overall graph data across subpaths
output_data_list = [] 
for i_group in range(len(which)):
    group_averaged_patients = []
    for i_patient in range(len(patients)):
        patient_dict = {}
        for p_max in p_max_or_thr_list:
            sum_pred = 0
            sum_delay = 0
            count = 0
            for i_subpath in range(len(subpaths_dict[which[i_group]])):
                sub_data = output_data_all[i_group][i_subpath][i_patient]
                if p_max in sub_data:
                    sum_pred += sub_data[p_max][0]
                    sum_delay += sub_data[p_max][1]
                    count += 1
            if count > 0:
                patient_dict[p_max] = [sum_pred/count, sum_delay/count]
        group_averaged_patients.append(patient_dict)
    output_data_list.append(group_averaged_patients)

# 3. Select best data point per subpath, then average those BEST points
best_outputs = np.zeros((len(which), len(patients), 3))
for i_group in range(len(which)):
    for i_patient in range(len(patients)):
        sum_best_p_max = 0
        sum_best_pred = 0
        sum_best_delay = 0
        num_subpaths = len(subpaths_dict[which[i_group]])
        
        # Skip processing if no valid subpaths were found for this group to prevent division by zero
        if num_subpaths == 0:
            continue
            
        for i_subpath in range(num_subpaths):
            best_pred_rate = 0
            best_delay = 0
            best_p_max = 0
            sub_data = output_data_all[i_group][i_subpath][i_patient]
            
            for p_max in p_max_or_thr_list:
                if p_max in sub_data:
                    pred = sub_data[p_max][0]
                    delay = sub_data[p_max][1]
                    # Betterness criterion
                    if (pred > best_pred_rate) or ((pred == best_pred_rate) and (delay < best_delay)): 
                        best_p_max = p_max
                        best_pred_rate = pred
                        best_delay = delay
                        
            sum_best_p_max += best_p_max
            sum_best_pred += best_pred_rate
            sum_best_delay += best_delay
            
        # Average the bests across subpaths for this patient
        best_outputs[i_group][i_patient][0] = sum_best_p_max / num_subpaths
        best_outputs[i_group][i_patient][1] = sum_best_pred / num_subpaths
        best_outputs[i_group][i_patient][2] = sum_best_delay / num_subpaths

# 4. Average over best points for each group (over patients)
av_best_outputs = np.zeros((len(which),3))
for i_group in range(len(which)):
    sum_p_max = 0
    sum_pred_rate = 0
    sum_delay = 0
    for i_patient in range(len(patients)):
        sum_p_max += best_outputs[i_group][i_patient][0]
        sum_pred_rate += best_outputs[i_group][i_patient][1]
        sum_delay += best_outputs[i_group][i_patient][2]
    av_best_outputs[i_group][0] = sum_p_max/len(patients)
    av_best_outputs[i_group][1] = sum_pred_rate/len(patients)
    av_best_outputs[i_group][2] = sum_delay/len(patients)



#Average over patients:
av_outputs = np.zeros((len(which), len(p_max_or_thr_list), 2))
for i_group in range(len(which)):
    nb_patients_per_pmax = [0]*len(p_max_or_thr_list)
    for i_p_max in range(len(p_max_or_thr_list)):
        p_max = p_max_or_thr_list[i_p_max]
        sum_pred_rate = 0
        sum_mean_delay = 0
        for i_patient in range(len(patients)):
            if (p_max in output_data_list[i_group][i_patient].keys()):
                nb_patients_per_pmax[i_p_max] += 1
                sum_pred_rate += output_data_list[i_group][i_patient][p_max][0]
                sum_mean_delay += output_data_list[i_group][i_patient][p_max][1]
        if (nb_patients_per_pmax[i_p_max] == 0):
            av_outputs[i_group][i_p_max][0] = 0
            av_outputs[i_group][i_p_max][1] = 0
        else:
            av_outputs[i_group][i_p_max][0] = sum_pred_rate/nb_patients_per_pmax[i_p_max]
            av_outputs[i_group][i_p_max][1] = sum_mean_delay/nb_patients_per_pmax[i_p_max]


#Dense average from https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8584751
detec_delay_dense = {2: 15.1, 4: 34.5, 5: 20.9, 6: 6.3, 8: 13.2, 11: 7.0, 13: 10.0, 16: 32.3}
av_dense_detc_delay = 0
for patient in patients:
    av_dense_detc_delay += detec_delay_dense[patient]
av_dense_detc_delay /= len(patients)

#Make graphs:
save_dir = Path(__file__).parent.parent / 'Figures7'
save_dir.mkdir(exist_ok=True, parents=True)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

lwidth = 3

ax2.plot(p_max_or_thr_list_, av_outputs[0,:,0]*100, label='sparse baseline', color='tab:blue', linewidth=lwidth, zorder=1)
ax2.plot(p_max_or_thr_list_, av_outputs[0,:,0]*100, ':', label='sparse + CompIM', color='tab:orange', linewidth=lwidth, zorder=1)
ax2.plot(p_max_or_thr_list_, av_outputs[1,:,0]*100, label='sparse + CompIM + thin. opt.', color='tab:purple', linewidth=lwidth, zorder=1)
ax2.plot(p_max_or_thr_list_, av_outputs[2,:,0]*100, label='sparse + IM-free + thin. opt.', color='tab:green', linewidth=lwidth, zorder=1)

#Dense
ax2.plot(p_max_or_thr_list_, [100]*len(p_max_or_thr_list), 'r', label='dense baseline [2]', linewidth=lwidth, zorder=1)


#Average best
ax2.scatter(av_best_outputs[0][0]*100, av_best_outputs[0][1]*100, marker='s', c='tab:blue', edgecolors='gold', linewidth=1.5, s=75, label='personalized threshold: sparse baseline', zorder=3)
ax2.scatter(av_best_outputs[1][0]*100, av_best_outputs[1][1]*100, marker='*', c='tab:purple', edgecolors='gold', linewidth=1.5, s=200, label='personalized threshold: sparse + CompIM + thin. opt.', zorder=3)
ax2.scatter(av_best_outputs[2][0]*100, av_best_outputs[2][1]*100, marker='v', c='tab:green', edgecolors='gold', linewidth=1.5, s=90, label='personalized threshold: sparse + IM-free + thin. opt.', zorder=3)
print((av_best_outputs[0][0]*100, av_best_outputs[0][1]))
print((av_best_outputs[1][0]*100, av_best_outputs[1][1]))
print((av_best_outputs[2][0]*100, av_best_outputs[2][1])) 

ax2.set_xlabel('Density after bundling')
ax2.set_ylabel('%',labelpad=-5)
ax2.set_title('(b) Average detection accuracy')
ax2.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)

#Axis limits
ax2.set_xlim(0,60)
ax2.set_ylim(80,100.2)

# Detection Delay

# prediction rate == 0, then dont print any detection delay
p_max_or_thr_list_group = []
detection_delays_group = []
for i_group in range(len(which)):
    p_max_or_thr_list_group.append([])
    detection_delays_group.append([])
    for i_p_max in range(len(p_max_or_thr_list_)):
        if (av_outputs[i_group,i_p_max,0] != 0):
            p_max_or_thr_list_group[-1].append(p_max_or_thr_list_[i_p_max])
            detection_delays_group[-1].append(av_outputs[i_group,i_p_max,1])


ax1.plot(p_max_or_thr_list_group[0], detection_delays_group[0], color='tab:blue', linewidth=lwidth, zorder=1)
ax1.plot(p_max_or_thr_list_group[0], detection_delays_group[0], ':', color='tab:orange', linewidth=lwidth, zorder=1)
ax1.plot(p_max_or_thr_list_group[1], detection_delays_group[1], color='tab:purple', linewidth=lwidth, zorder=1)
ax1.plot(p_max_or_thr_list_group[2], detection_delays_group[2], color='tab:green', linewidth=lwidth, zorder=1)
ax1.plot(p_max_or_thr_list_, [av_dense_detc_delay]*len(p_max_or_thr_list_), 'r', linewidth=lwidth, zorder=1)

#Average best
ax1.scatter(av_best_outputs[0][0]*100, av_best_outputs[0][2], marker='s', c='tab:blue', edgecolors='gold', linewidth=1.5, s=75, zorder=3)
ax1.scatter(av_best_outputs[1][0]*100, av_best_outputs[1][2], marker='*', c='tab:purple', edgecolors='gold', linewidth=1.5, s=200, zorder=3)
ax1.scatter(av_best_outputs[2][0]*100, av_best_outputs[2][2], marker='v', c='tab:green', edgecolors='gold', linewidth=1.5, s=90, zorder=3)
print((av_best_outputs[0][0]*100, av_best_outputs[0][2]))
print((av_best_outputs[1][0]*100, av_best_outputs[1][2]))
print((av_best_outputs[2][0]*100, av_best_outputs[2][2])) 
print([av_dense_detc_delay]*len(p_max_or_thr_list_))


ax1.set_title('(a) Average detection delay')
ax1.set_xlabel('Density after bundling')
ax1.set_ylabel('Seconds', labelpad=0)

# 1. Increased the y-coordinate from 1.01 to 1.08 to push the legend up
fig.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.50, 1.02), frameon=False, fontsize=11, columnspacing=0.5)

ax1.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)

# Axis limits
ax1.set_xlim(0,60)
ax1.set_ylim(15,30)

# 2. Lowered the top boundary of the subplots from 0.85 to 0.78 to give the legend more room
plt.tight_layout(rect=[0, 0, 1, 0.78])
plt.subplots_adjust(wspace=0.14)

# Dynamically name the output based on what runs were found
if len(which) > 0 and len(subpaths_dict[which[0]]) > 0:
    run_str = "".join(subpaths_dict[which[0]])
else:
    run_str = "empty"

if (len(patients) == 1):
    fig.savefig(save_dir / f'detec_delay_{patients[0]}.pdf', bbox_inches='tight')
else:
    fig.savefig(save_dir / f'summary_runs_{run_str}.pdf', bbox_inches='tight')
    fig.savefig(save_dir / f'summary_runs_{run_str}.png', dpi=1000, bbox_inches='tight')