# Quantum Energy Levels in Central Potentials

This repository contains a single MATLAB script that models the quantum energy levels of central potential wells using the finite difference method.  

Currently, the code implements and compares:

- **Harmonic Oscillator**  
- **Infinite Square Well**  

Finite square well (FSW) and spin–orbit coupling (SOC) versions have been removed to keep the project focused and minimal.

## How it works
- Constructs the kinetic energy operator using finite differences.  
- Builds the potential energy matrix for each potential.  
- Solves for the lowest eigenvalues and eigenfunctions using MATLAB's sparse eigensolver (`eigs`).  
- Plots the resulting energy levels for comparison.

## Requirements
- MATLAB R2022a or later (earlier versions may work).  
- No toolboxes beyond base MATLAB are required.  

## Usage
Clone the repository and run the script:

```matlab
quantum_central_potentials.m
```

This will compute and display the lowest energy levels for the harmonic oscillator and the infinite square well.

## Author
Developed by **Yunki Yau**  
GitHub: [yunkiyau](https://github.com/yunkiyau)  

---
