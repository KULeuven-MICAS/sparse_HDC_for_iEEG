import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import numpy as np
import math

# ==========================================
# --- FONT SIZE & STYLE CONFIGURATION ---
# ==========================================
plt.rc('font', family='DejaVu Math TeX Gyre')

plt.rcParams.update({
    'axes.titlesize': 26,     
    'axes.labelsize': 22,     
    'xtick.labelsize': 20,    
    'ytick.labelsize': 20,    
    'legend.fontsize': 22,    
    'axes.titlepad': -55      
})
# ==========================================
# ==========================================

def format_func(value, tick_number=None):
    if value == 0:
        return "0"
    num_thousands = 0 if abs(value) < 1000 else math.floor(math.log10(abs(value))/3)
    value = round(value / 1000**num_thousands, 2)
    return f'{value:g}'+' KMGTPEZY'[num_thousands]

# --- Helper Function for Improvement Arrows ---
def add_improvement_arrows(ax, totals):
    base_val = totals[0]
    # Draw a reference dashed line from the baseline bar across the plot
    ax.plot([0, len(totals)-1], [base_val, base_val], color='black', linestyle='--', alpha=0.4, zorder=1)
    
    for i in range(1, len(totals)):
        curr_val = totals[i]
        factor = base_val / curr_val
        
        ax.annotate(
            '',
            xy=(i, curr_val), 
            xytext=(i, base_val),
            arrowprops=dict(arrowstyle="->,head_length=1.2,head_width=0.6", 
                            color='#222222', lw=4.0, shrinkA=0, shrinkB=0),
            zorder=4
        )
        
        mid_y = (base_val + curr_val) / 2
        ax.text(i - 0.04, mid_y, f'{factor:.1f}×', 
                va='center', ha='right', fontsize=24, fontweight='bold',
                color='#222222',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=2),
                zorder=5)

# Define the data
categories = ['IM', 'One-hot decode', 'CompIM', 'Binding', 'Bund. (spatial)', 'Bund. (temporal)', 'Sim. search', 'Gen. overhead']
v1_total_area = np.array([189960.43, 128697.21, 96833.66, 44574, 23171])
detailed_categories = ['Dense\nbaseline', 'Sparse\nbaseline', 'Sparse\n+CompIM', 'Sparse\n+CompIM\n+thin opt.', 'Sparse\nIM-free\n+thin opt.']
v1_breakdown_area = np.array([
    [44.38,  1.86,  0,     0,     0],
    [0,      15.90, 0,     0,     0],
    [0,      0,     6.16,  13.38, 0],
    [1.79,   22.2,  14.05, 30.54, 5.42],
    [43.65,  44.9,  59.67, 12.62, 11.79],
    [7.76,   12.2,  16.21, 33.45, 64.352],
    [1.40,   1.15,  1.53,  5.25,  10.1],
    [1.02,   1.79,  2.37,  4.33,  8.33]
])

v1_total_energy = np.array([93.7, 21.6, 13.554, 10.0, 7.01])
v1_breakdown_energy = np.array([
    [16.5,  1.42,  0,     0,     0],
    [0,     23.7,  0,     0,     0],
    [0,     0,     8.06,  10.92, 0],
    [3.01,  27.6,  12.76, 17.29, 2.04],
    [66.1,  14,    22.31, 5.78,  4.78],
    [14.2,  32.3,  51.47, 64.19, 90.88],
    [0.108, 0.203, 0.32,  0.83,   1.17],
    [0.082, 0.777, 1.246, 0.985,  1.13]
])

area_colors = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3', '#937860', '#DA8BC3', '#8C8C8C']

fig = plt.figure(figsize=(24, 6))

gs1 = gridspec.GridSpec(1, 2)
gs1.update(wspace=0.10, hspace=0.1)  
ax1 = fig.add_subplot(gs1[0, 0])
ax2 = fig.add_subplot(gs1[0, 1])

width = 0.65 
x = np.arange(len(detailed_categories))

# --- Left Plot: Energy Breakdown (ax1) ---
bottom_energy = np.zeros(len(detailed_categories))
bar_handles = [] # <-- FIX: Create a list to store the pure bar handles

for i in range(len(v1_breakdown_energy)):
    bar = ax1.bar(x, v1_breakdown_energy[i] * v1_total_energy / 100, width, bottom=bottom_energy, 
                  color=area_colors[i], edgecolor='black', alpha=0.85, zorder=3)
    bar_handles.append(bar) # <-- FIX: Save the bar handle
    bottom_energy += v1_breakdown_energy[i] * v1_total_energy / 100

ax1.set_xticks(x)
ax1.set_xticklabels(detailed_categories, rotation=0)
ax1.tick_params(axis='y')
ax1.set_title('(a) Energy Breakdown', y=1.13)
ax1.set_ylabel('Energy/pred [nJ]')

add_improvement_arrows(ax1, v1_total_energy)
ax1.set_ylim(0, v1_total_energy[0] * 1.15) 
ax1.yaxis.set_major_locator(MultipleLocator(10))

# --- Right Plot: Area Breakdown (ax2) ---
bottom_area = np.zeros(len(detailed_categories))
for i in range(len(v1_breakdown_area)):
    ax2.bar(x, v1_breakdown_area[i] * v1_total_area / 100, width, bottom=bottom_area, 
            color=area_colors[i], edgecolor='black', alpha=0.85, zorder=3)
    bottom_area += v1_breakdown_area[i] * v1_total_area / 100

ax2.set_xticks(x)
ax2.set_xticklabels(detailed_categories, rotation=0)
ax2.tick_params(axis='y')
ax2.set_title('(b) Area Breakdown', y=1.13)
ax2.set_ylabel('Area [μm²]')
ax2.yaxis.set_major_formatter(plt.FuncFormatter(format_func))

add_improvement_arrows(ax2, v1_total_area)
ax2.set_ylim(0, v1_total_area[0] * 1.15)

# --- Clean up axes ---
for ax in [ax1, ax2]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=0.5, zorder=0)

# --- THE FIX: Pass bar_handles explicitly to the legend ---
fig.legend(bar_handles, categories, ncol=4, loc='upper center', bbox_to_anchor=(0.5, 1.08), 
           frameon=False, columnspacing=1.0)

for label in ax1.get_xticklabels() + ax1.get_yticklabels() + ax2.get_xticklabels() + ax2.get_yticklabels():
    label.set_alpha(0.85)

# Save the picture
plt.savefig("v2_area_energy_breakdown_J.png", bbox_inches='tight', dpi=300)
plt.savefig("v2_area_energy_breakdown_J.svg", bbox_inches='tight')
plt.savefig("v2_area_energy_breakdown_J.pdf", bbox_inches='tight', format='pdf')

plt.show()