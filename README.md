# Decoherence-of-Majorana-Qubits-From-1-f-Noise
Covariance Matrix Method Code for A. Alase, M. C. Goffage, M. C. Cassidy, S. N. Coppersmith
# Implementation of QPP_Library for: Decoherence in Majoraba Qubits by 1/f Noise

This repository contains the code used to generate all numerical results presented in our paper "A. Alase, M. C. Goffage, M. C. Cassidy, S. N. Coppersmith, Decoherence in Majorana Qubits by 1/f Noise (2025)". This code uses the covariance matrix method to calculate the time evolution of the tetron qubit, composed of two Kitaev chains, in the presence of a two-level-fluctuator (TLF) and returns the probability of exciting a quasiparticle pair in a single Kiteav chain of the qubit. See Appendix B of our earlier paper, "M. C. Goffage, A. Alase, M. C. Cassidy, S. N. Coppersmith, Leakage at zero temperature from changes in chemical potential in Majorana qubits, arXiv:2504.17485 (2025)", for further details on the method. 

Note that since our numerical results are statistical averages over multiple two-level-fluctuator (TLF) noise realisations, we expect minor differences in the plots between consecutive runs. 

## How to Run the Code
Locally: download the zip file titled "download_for_local_run" and unzip the contents in a local directory. Open the "code" directory, then open "Run_all_Figures.m" in MATLAB. Select the "editor tab" in Matlab and then select "run" to generate all figures in the paper (if prompted to Change Folder to the current directory by the Matlab editor, choose Change Folder, which is the highlighted option). The key figure, Fig. 2b, is generated first and takes under 2 hours. You may choose "stop run" after this time to view the key figure in the "results" directory. Generating the remaining figures may take more than 24 hours.


## Structure of the Code
This code runs our full covariance matrix code and generates all figures in our paper which present numerical results. In this order, the code generates the figures: Fig. 2c, Fig. 2b, Fig. S1, Fig. 3c, Fig. 3d, Fig. S2, and then Fig. S4. 

## Environment and Dependencies

Requires MATLAB2024a or more recent. 

### Output
All figures are outputted into the folder "results" and are titled according to their corresponding figure and panel in the main text or Supplementary Information of our paper. Additionally, the corresponding MATLAB workspaces generated for each figure are saved in the results folder. 

## Troubleshooting

If you run into any issues or have any questions, please contact m.goffage@unsw.edu.au. 
