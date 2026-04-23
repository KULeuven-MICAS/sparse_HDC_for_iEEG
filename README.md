# iEEG Seizure Detection with a Sparse Hyperdimensional Computing Accelerator

Official repository for S. Cuyckens, R. Antonio, C. Fang, M. Verhelst, "iEEG Seizure Detection with a Sparse Hyperdimensional Computing Accelerator," *20th International Conference on PhD Research in Microelectronics and Electronics (PRIME)*, 2025.

<p align="center">
  <img src="https://img.shields.io/badge/Language-SystemVerilog%20%7C%20Python%20%7C%20MATLAB-blue" alt="Languages"/>
  <img src="https://img.shields.io/badge/Venue-PRIME%202025-green" alt="PRIME 2025"/>
  <img src="https://img.shields.io/badge/Topic-Sparse%20HDC%20%7C%20iEEG-orange" alt="Sparse HDC"/>
</p>

This repository contains the hardware accelerator designs (SystemVerilog) and software implementations (Python, MATLAB) for sparse HDC-based seizure detection from intracranial EEG.

## Repository Structure

```
Hardware/
├── DenseHDC.sv                    # Dense HDC baseline
├── SparseHDC.sv                   # Sparse HDC with CompIM and simplified spatial bundling [1]
├── Shift_bind_optimized.sv        # IM-free architecture with shift binding [2]
└── Segm_shift_bind_optimized.sv   # Segmented shift binding with channel folding [2]

Software/
├── HDC_functions.py               # Core HDC operations (binding, bundling, similarity)
├── even_more_efficient_eeg.py     # Training and inference pipeline
├── LBP_from_EEG.m                # LBP feature extraction from raw iEEG (MATLAB)
├── algorithmic_performance_graph.py     # Algorithmic performance figures
├── performance_graph_better_figures.py  # Publication-quality figures
├── Breakdown_P_A_CHIPS.py              # Power/area breakdown figures
├── Breakdown_P_A_CHIPS_CVFF.py         # Breakdown with channel/vector folding
├── Breakdown_P_A_CHIPS_CVFF_IM-free.py # Breakdown for IM-free variants
└── main.py
```

## Hardware

All designs are parameterized with D=1024 (hypervector dimension), 64 EEG channels, and 6-bit LBP codes. The four designs correspond to the dense baseline and the successive optimizations introduced in [1] and [2].

## Software

### Prerequisites

- **Python 3** with NumPy, SciPy, Matplotlib
- **MATLAB** for LBP feature extraction

### Data

The software expects LBP-extracted seizure data organized as:
```
no_backup/
├── Pat2/, Pat4/, Pat5/, Pat6/, Pat8/, Pat11/, Pat13/, Pat16/
│   ├── LBP_1.mat
│   ├── LBP_2.mat
│   └── ...
```

Use `LBP_from_EEG.m` to extract LBP features from raw iEEG `.mat` files. The main training and inference pipeline is in `even_more_efficient_eeg.py`.

## Citation

If you use this code in your work, please cite:

<details>
<summary>BibTeX</summary>
<p>

```bibtex
@INPROCEEDINGS{11203735,
  author={Cuyckens, Stef and Antonio, Ryan and Fang, Chao and Verhelst, Marian},
  booktitle={2025 20th International Conference on PhD Research in Microelectronics and Electronics (PRIME)}, 
  title={iEEG Seizure Detection with a Sparse Hyperdimensional Computing Accelerator}, 
  year={2025},
  pages={1-4},
  doi={10.1109/PRIME66228.2025.11203735}
}

@misc{Cuyckens_Antonio_Fang_Verhelst_2026, title={Hardware design optimization of a sparse hyperdimensional computing accelerator for IEEG seizure detection}, url={https://www.mdpi.com/2674-0729/5/2/10}, journal={MDPI}, publisher={Multidisciplinary Digital Publishing Institute}, author={Cuyckens, Stef and Antonio, Ryan and Fang, Chao and Verhelst, Marian}, year={2026}, month={Apr}} 
```

</p>
</details>
