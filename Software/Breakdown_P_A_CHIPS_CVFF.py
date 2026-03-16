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
# ==========================================

def format_func(value, tick_number=None):
    if value == 0:
        return "0"
    num_thousands = 0 if abs(value) < 1000 else math.floor(math.log10(abs(value))/3)
    value = round(value / 1000**num_thousands, 2)
    return f'{value:g}'+' KMGTPEZY'[num_thousands]

# Define the area breakdown data for V1
categories = ['IM', 'One-hot decode', 'CompIM', 'Binding', 'Bund. (spatial)', 'Bund. (temporal)', 'Sim. search', 'Gen. overhead']

v1_total_area_CFF = np.array([44574, 35594, 31819, 28550, 25357, 23169, 21662])
v1_total_area_VFF = np.array([44574, 39468, 36786, 34877, 0, 0, 0])
detailed_categories = ['No FF','CFF|VFF\n= 2', 'CFF|VFF\n= 4', 'CFF|VFF\n= 8', 'CFF|VFF\n= 16', 'CFF|VFF\n= 32', 'CFF|VFF\n= 64']

v1_breakdown_area_CFF = np.array([
    [0,0,      0,     0,     0,     0,    0],
    [0,0,      0,     0,     0,     0,    0],
    [13.8,9,     6.43,  3.08,  0,     0,  0],
    [30.5,27.0,  26.8,  24.8,  19.0,  13.0, 8.51],
    [12.6,8.05,  4.17,  2.42,  0,     0,  0],
    [33.5,42.4,  47.4,  52.8,  59.4,  64.6, 69.0],
    [5.3,6.57,   7.36,  8.20,  9.23,  10.1, 10.8],
    [4.3,7.02,   7.85,  8.75,  12.4,  12.3, 10.5]
])

v1_breakdown_area_VFF = np.array([
    [0,0,     0,     0,     0, 0, 0],
    [0,0,     0,     0,     0, 0, 0],
    [13.8,16.2,  18.7,  21.0,  0, 0, 0],
    [30.5,24.6,  23.8,  20.6,  0, 0, 0],
    [12.6,7.13,  3.82,  2.02,  0, 0, 0],
    [33.5,39.7,  40.6,  42.5,  0, 0, 0],
    [5.3,5.93,   6.36,  6.71,  0, 0, 0],
    [4.3,6.42,   6.79,  7.15,  0, 0, 0]
])

# Define the energy consumption data for V1
v1_total_energy_CFF = np.array([10.5, 16.8, 21.7, 26.2, 37.0, 61.4, 82.1])
v1_total_energy_VFF = np.array([10.5, 19.5, 25.7, 33.0, 0, 0, 0])

v1_breakdown_energy_CFF = np.array([
    [0,0,     0,     0,     0,     0,     0],
    [0,0,     0,     0,     0,     0,     0],
    [10.7,30.4,  27.3,  24.69,  23,     19.87,   23.95],
    [16.7,42.86, 48.6,  47.34,  40.6,  34.27, 26.75],
    [5.52,8.78,  6,     5.83,   5.4,     5.2,     5.1],
    [65.5,11.8,  11.2,  13.36,  21.6,  25.47, 33.33],
    [0.8,0.88,  1.37,  2.26,   3.2,   3.87,  5.79],
    [0.77,5.26,  5.52,  6.52,  6.12,  11.0,  4.88]
])

v1_breakdown_energy_VFF = np.array([
    [0,0,     0,    0,     0, 0, 0],
    [0,0,     0,    0,     0, 0, 0],
    [10.7,24.7,  22.9, 20,  0, 0, 0],
    [16.7,36.6,  42.2, 37.3,  0, 0, 0],
    [5.52,6.37,  8.45, 8.88,  0, 0, 0],
    [65.5,26.3,  20.06, 26.9,  0, 0, 0],
    [0.8,0.774, 1.16, 1.9,  0, 0, 0],
    [0.77,5.28,  4.65, 5.06,  0, 0, 0]
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
width = 0.35 
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
plt.savefig("v2_area_energy_breakdown_J_CVFF.png", bbox_inches='tight', dpi=300)
plt.savefig("v2_area_energy_breakdown_J_CVFF.svg", bbox_inches='tight')
plt.savefig("v2_area_energy_breakdown_J_CVFF.pdf", bbox_inches='tight', format='pdf')

plt.show()