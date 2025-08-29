# Quantum Energy Levels in Central Potentials

This repository contains a single MATLAB script that models the first few quantum energy levels of central potential wells using the finite difference method. I did this as part of an assignment for my senior Advanced Quantum Physics course. The assignment also required creation of Finite square well models with and without spin–orbit coupling (SOC) which I have not included so that future students do not have access to the solution if a similar assignent is given again. The code computes the energy levels of:

- **Harmonic Oscillator**  
- **Infinite Square Well**  

These examples should give students a good idea on how other models can be computed. Please reach out to me at yunki.yau@gmail.com if you need more assistance and I would be happy to help. 

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
yunki.yau@gmail.com, yyau2516@uni.sydney.edu.au 
GitHub: [yunkiyau](https://github.com/yunkiyau)  

---
