# 🧬 Analysis Colorectal Cancer In Situ Sequencing Dataset with InSituPy

> ⚠️ **Note:** This repository is currently under active development. Some features or notebooks may not be fully functional yet. Contributions and feedback are welcome!

This repository contains Jupyter notebooks for the analysis of a colorectal cancer Xenium In Situ dataset. The primary goal is to demonstrate the capabilities and functionalities of [**InSituPy**](https://github.com/SpatialPathology/InSituPy).

## 📁 Repository Structure

```
├── notebooks/
│   ├── 01_data_loading.ipynb
│   ├── 02_preprocessing.ipynb
│   ├── 03_spatial_analysis.ipynb
│   └── 04_visualization.ipynb
├── data/ (not included)
├── README.md
└── sample_configuration.csv
```

## 🚀 Features

- Loading and exploring Xenium In Situ spatial transcriptomics data
- Preprocessing and quality control
- Spatial clustering and cell-type annotation
- Visualization of spatial gene expression patterns
- Demonstration of InSituPy's core functionalities

## 🛠️ Requirements

To run the demo notebooks, InSituPy has to be installed. For information on how to install it, see [here](https://github.com/SpatialPathology/InSituPy?tab=readme-ov-file#installation).

The data can be downloaded from [GEO]() and needs to be placed in the `data/` directory.

> ⚠️ **Note:** The Xenium output files on GEO are provided as `.tar.gz` archives. These must be extracted before use.

## 📊 Getting Started

1. Clone the repository:
   ```bash
   git clone https://github.com/jwrth/colorectalcancer-insitupy-demo.git
   cd colorectalcancer-insitupy-demo
   ```

2. Launch Jupyter Lab or Notebook:
   ```bash
   jupyter lab
   ```

3. Open the notebooks in the `notebooks/` folder and follow along.

## 📬 Contact

For questions or collaborations, feel free to open an issue or reach out via email.

---

This README was generated with the assistance of **Microsoft Copilot**.
