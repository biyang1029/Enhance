
# Enhance: Enhanced Heat Exchange Simulation Framework

## 🧩 Overview
**Enhance** is a C++-based numerical simulation framework designed for studying enhanced heat transfer performance in coaxial or tube-type heat exchangers.  
The code adopts a modular finite volume method (FVM) structure, enabling efficient and flexible simulation of heat transfer processes in geothermal or solar-assisted systems.

This repository is part of a postdoctoral research project focusing on **the optimization of enhanced heat exchange structures for medium–deep geothermal seasonal storage systems**.

---

## ⚙️ Features
- 🧠 **Finite Volume Solver** for transient 1D/2D/3D heat conduction and convection.
- 🔄 **Modular Structure** including:
  - Fluid flow and heat exchange module (inner and outer pipes)
  - Soil/rock thermal conduction module
  - Solar/heat pump coupling interface
  - Non-uniform grid generator and implicit solver
- 📊 **Parametric Simulation** of pipe geometry, mass flow rate, temperature, and thermal conductivity.
- ⚡ **High Performance**: optimized for large-scale computation with adjustable grid density.
- 🧾 **Configurable Inputs**: all physical and operational parameters can be modified through input files.

---

## 🧱 Code Structure
```
Enhance/
│
├── src/                          # Core C++ source files
│   ├── main.cpp                  # Main program entry (simulation control)
│   ├── dataconfig.cpp/.h         # Input reading, data initialization
│   ├── groundmodule.cpp/.h       # 3D soil heat transfer model (FVM solver)
│   ├── coaxialwell.cpp/.h        # Coaxial pipe heat transfer model
│   ├── heatpumpmodule.cpp/.h     # Heat pump performance and COP calculation
│   ├── solarmodule.cpp/.h        # Solar collector and auxiliary heating model
│   ├── simulationcontroller.cpp/.h # Overall system control and coupling logic
│   ├── progressviz.cpp/.h        # Console progress and result visualization
│   ├── logger.cpp/.h             # Logging and output management
│   └── utils.cpp/.h              # Math utilities, interpolation, matrix operations
│
├── include/                      # Header files (if separated from src/)
│
├── input/                        # User input and configuration files
│   ├── config.txt                # System parameter definitions
│   ├── weather_data.txt          # Meteorological data (optional)
│   └── soil_properties.txt       # Soil layer and thermal properties
│
├── output/                       # Simulation outputs
│   ├── temperature_field/        # 3D soil temperature snapshots
│   ├── performance_log.txt       # System performance summary
│   └── progress_log.txt          # Iteration and convergence information
│
├── data/                         # Experimental or validation datasets (optional)
│
├── scripts/                      # Post-processing or plotting scripts
│
├── README.md                     # Project documentation
└── .gitignore                    # Ignore rules for build cache, VS, binaries
```
---

## 🧮 Dependencies
- **C++17** or later  
- **CMake ≥ 3.15** (recommended for build)  
- **Visual Studio 2022 / g++ / Clang**  
- Optional:  
  - [Eigen](https://eigen.tuxfamily.org/) for matrix operations  
  - [Matplotlib C++](https://github.com/lava/matplotlib-cpp) for visualization  

---

## 🚀 Build and Run
### 1️⃣ Compile
Using CMake:
```bash
mkdir build
cd build
cmake ..
cmake --build .
```
Or directly with Visual Studio:

Open the .sln file

Set main.cpp as startup project

Press F5 to build and run

2️⃣ Run Simulation

Edit your input file (e.g. input/config.txt) to define parameters:
```ini
PipeLength = 2500
MassFlowRate = 0.2
Tinlet = 90
SoilConductivity = 2.0
TimeStep = 60
```
Then execute:
```bash
./enhance.exe
```
Outputs will be written to /output/ folder.

📈 Example Result

Outlet water temperature evolution

Soil temperature distribution

Energy efficiency and heat storage capacity comparison

(Plots can be generated using Python or MATLAB for post-processing.)

📚 Citation


🧑‍💻 Author

Yang Bi
Postdoctoral Researcher, Shanghai Jiao Tong University
Email: yang.bi@sjtu.edu.cn

📜 License

MIT License © 2025 Bi Yang
Free for research and educational use.


