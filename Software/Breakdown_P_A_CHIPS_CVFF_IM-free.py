import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import matplotlib.patches as mpatches
import numpy as np
import math

# ==========================================
# --- FONT SIZE & STYLE CONFIGURATION ---
# ==========================================
plt.rc('font', family='DejaVu Math TeX Gyre')

plt.rcParams.update({
    'axes.titlesize': 26,     # Font size for subplot titles (a) and (b)
    'axes.labelsize': 22,     # Font size for X/Y axis labels
    'xtick.labelsize': 20,    # Font size for X-axis tick labels (categories)
    'ytick.labelsize': 20,    # Font size for Y-axis tick labels (numbers)
    'legend.fontsize': 22,    # Font size for the shared legend
    'axes.titlepad': -45      # <-- Negative value pushes the title down
})
# ==========================================

def format_func(value, tick_number=None):
    if value == 0:
        return "0"
    num_thousands = 0 if abs(value) < 1000 else math.floor(math.log10(abs(value))/3)
    value = round(value / 1000**num_thousands, 2)
    return f'{value:g}'+' KMGTPEZY'[num_thousands]

# Define the area breakdown data for V1
categories = ['IM', 'One-hot decode', 'CompIM', 'Binding', 'Bund. (spatial)', 'Bund. (temporal)', 'Sim. search', 'Gen. overhead']

v1_total_area_CFF = np.array([23171, 26151, 26194, 25529, 24339, 22976, 22078])
detailed_categories = ['No FF', 'CFF|VFF\n= 2', 'CFF|VFF\n= 4', 'CFF|VFF\n= 8', 'CFF|VFF\n= 16', 'CFF|VFF\n= 32', 'CFF|VFF\n= 64']
v1_total_area_VFF = np.array([23171, 26285, 26099, 24824, 0, 0, 0])


v1_breakdown_area_CFF = np.array([
    [0,0,      0,     0,     0,     0,    0],
    [0,0,      0,     0,     0,     0,    0],
    [0,0,      0,     0,     0,     0,    0],
    [5.42,12,   15.5,  14.8,   13.3,  10,   6.73],
    [11.79,12.9,  9.89,  7.91, 6.57,  6.03, 6.04],
    [64.35,57.7,  57.5,  59,   61.9,  65.1, 67.8],
    [10.1,8.95,   8.93,  9.17,  9.62,  10.2, 10.6],
    [8.33,8.47,   8.19,  9.12,  8.59,  8.64, 8.86]
])

v1_breakdown_area_VFF = np.array([
    [0,0,      0,     0,     0,     0,    0],
    [0,0,      0,     0,     0,     0,    0],
    [0,0,      0,     0,     0,     0,    0],
    [5.42,14.1,  19.28,  18.03,  0,  0, 0],
    [11.79,7.83, 5.02,  2.83,  0,  0, 0],
    [64.35,59.66, 57.15, 59.64, 0,  0, 0],
    [10.1,8.91,  8.97,  9.43,  0,  0, 0],
    [8.33,9.50,  9.58,  10.1,  0,  0, 0]
])

# Define the energy consumption data for V1
v1_total_energy_CFF = np.array([7.01, 6.55, 10.1, 16.1, 27.2, 52.8, 99.3])
v1_total_energy_VFF = np.array([7.01, 7.83, 7.95, 11.4, 0, 0, 0])

v1_breakdown_energy_CFF = np.array([
    [0,0,     0,     0,     0,     0,     0],
    [0,0,     0,     0,     0,     0,     0],
    [0,0,     0,     0,     0,     0,     0],
    [2.04,19,    24.2,  23.4,  17.17,   18.3,  15.50],
    [4.78,42.34,  42.2,  43.46,  43.07,   44,     46],
    [90.9,25.7,  21.8,  20.1,  27.7,  26.4, 27.1],
    [1.17,2.27,   2.9,  3.68,  4.37,   4.5,  4.79],
    [11.3,10.5,   8.80,  9.40,  7.68,   6.82,  6.62]
])

v1_breakdown_energy_VFF = np.array([
    [0,0,     0,     0,     0,     0,     0],
    [0,0,     0,     0,     0,     0,     0],
    [0,0,     0,     0,     0,     0,     0],
    [2.04,14.77, 30.03, 6.72,  0,   0,  0],
    [4.78,8.63,  13.4,  5.77,  0,   0,     0],
    [90.9,63.01, 37.50, 57.7,  0,  0, 0],
    [1.17,1.93,  3.78,  8.15,  0,   0,  0],
    [11.3,11.7,  15.3,  21.7,  0,   0,  0]
])

# A modern, muted, colorblind-friendly palette
area_colors = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3', '#937860', '#DA8BC3', '#8C8C8C']

# Create the figure and subplots
fig = plt.figure(figsize=(24, 6))

gs1 = gridspec.GridSpec(1, 2)
gs1.update(wspace=0.10, hspace=0.1)  
ax1 = fig.add_subplot(gs1[0, 0])
ax2 = fig.add_subplot(gs1[0, 1])

# Narrower width to allow two bars per category
# Narrower width to allow two bars per category
width = 0.35 

# Multiply by a factor > 1 to increase the spacing between pairs
# Try 1.5 for moderate spacing, or 2.0 for wide spacing
x = np.arange(len(detailed_categories)) 

# --- Right Plot: Area Breakdown (ax2) ---
bottom_area_cff = np.zeros(len(detailed_categories))
bottom_area_vff = np.zeros(len(detailed_categories))

for i in range(len(categories)):
    # Plot CFF bars (Shifted left)
    ax2.bar(x - width/2, v1_breakdown_area_CFF[i] * v1_total_area_CFF / 100, width, 
            bottom=bottom_area_cff, color=area_colors[i], edgecolor='black', alpha=0.85, zorder=3)
    bottom_area_cff += v1_breakdown_area_CFF[i] * v1_total_area_CFF / 100

    # Plot VFF bars (Shifted right, adding hatch to distinguish)
    ax2.bar(x + width/2, v1_breakdown_area_VFF[i] * v1_total_area_VFF / 100, width, 
            bottom=bottom_area_vff, color=area_colors[i], edgecolor='black', alpha=0.85, zorder=3)
    bottom_area_vff += v1_breakdown_area_VFF[i] * v1_total_area_VFF / 100

ax2.set_xticks(x)
ax2.set_xticklabels(detailed_categories, rotation=0)
ax2.tick_params(axis='y')
ax2.set_title('(b) Area Breakdown', y=1.13)
ax2.set_ylabel('Area [μm²]')
ax2.yaxis.set_major_formatter(plt.FuncFormatter(format_func))

# --- Left Plot: Energy Breakdown (ax1) ---
bottom_energy_cff = np.zeros(len(detailed_categories))
bottom_energy_vff = np.zeros(len(detailed_categories))

for i in range(len(categories)):
    # Plot CFF bars (Shifted left)
    ax1.bar(x - width/2, v1_breakdown_energy_CFF[i] * v1_total_energy_CFF / 100, width, 
            bottom=bottom_energy_cff, color=area_colors[i], edgecolor='black', alpha=0.85, zorder=3)
    bottom_energy_cff += v1_breakdown_energy_CFF[i] * v1_total_energy_CFF / 100

    # Plot VFF bars (Shifted right, adding hatch to distinguish)
    ax1.bar(x + width/2, v1_breakdown_energy_VFF[i] * v1_total_energy_VFF / 100, width, 
            bottom=bottom_energy_vff, color=area_colors[i], edgecolor='black', alpha=0.85, zorder=3)
    bottom_energy_vff += v1_breakdown_energy_VFF[i] * v1_total_energy_VFF / 100

ax1.set_xticks(x)
ax1.set_xticklabels(detailed_categories, rotation=0)
ax1.tick_params(axis='y')
ax1.set_title('(a) Energy Breakdown', y=1.13)
ax1.set_ylabel('Energy/pred [nJ]')

# --- Force y-axis to tick every 10 units for the energy plot ---
log = False
if log:
    ax1.set_yscale('log')
else:
    ax1.yaxis.set_major_locator(MultipleLocator(10))

# --- Clean up axes (Despining & Grids) ---
for ax in [ax1, ax2]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=0.5, zorder=0)

# ==========================================
# --- LEGEND FIX ---
# ==========================================
# Main legend for the breakdown categories (explicitly mapped to colors, now perfectly centered)
category_patches = [mpatches.Patch(facecolor=area_colors[i], edgecolor='black', alpha=0.85, label=categories[i]) for i in range(len(categories))]
fig.legend(handles=category_patches, ncol=4, loc='lower center', bbox_to_anchor=(0.5, 0.95), 
           frameon=False, columnspacing=1.0)
# ==========================================

for label in ax1.get_xticklabels() + ax1.get_yticklabels() + ax2.get_xticklabels() + ax2.get_yticklabels():
    label.set_alpha(0.85)

# Save the picture
plt.savefig("v2_area_energy_breakdown_J_CVFF_sb.png", bbox_inches='tight', dpi=300)
plt.savefig("v2_area_energy_breakdown_J_CVFF_sb.svg", bbox_inches='tight')
plt.savefig("v2_area_energy_breakdown_J_CVFF_sb.pdf", bbox_inches='tight', format='pdf')

plt.show()