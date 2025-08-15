# Decoherence-of-Majorana-Qubits-From-1-f-Noise
Code for all numerical calculations in "A. Alase, M. C. Goffage, M. C. Cassidy, S. N. Coppersmith (2025), Decoherence of Majorana Qubits from 1/f Noise, arxiv preprint: 	arXiv:2506.22394."
# Implementation of QPP_Library for: Decoherence in Majoraba Qubits by 1/f Noise

This repository contains the code used to generate all numerical results presented in our paper "A. Alase, M. C. Goffage, M. C. Cassidy, S. N. Coppersmith, Decoherence in Majorana Qubits by 1/f Noise (2025)". This code uses the covariance matrix method to calculate the time evolution of the tetron qubit, composed of two Kitaev chains, in the presence of a two-level-fluctuator (TLF) and returns the probability of exciting a quasiparticle pair in a single Kiteav chain of the qubit. See Methods and Sec. 3 of Supplementary Information of our paper for further details on our calculations. For a more in-depth overview of our numerical package, called "QPP_library", please see Appendix B of our earlier paper, "M. C. Goffage, A. Alase, M. C. Cassidy, S. N. Coppersmith, Leakage at zero temperature from changes in chemical potential in Majorana qubits, arXiv:2504.17485 (2025)", which employs the same library as the submitted paper. 

Note that since our numerical results are statistical averages over multiple two-level-fluctuator (TLF) noise realisations, we expect minor differences in the plots between consecutive runs. 

## How to Run the Code
1. Click the green **"Code"** drop-down menu at the top of this page.  
2. Select **"Download ZIP"**.  
3. Unzip the contents into a local directory.  
4. Open **`Run_Code_Demo.m`** or **`Run_all_Figures.m`** in MATLAB.  
5. In MATLAB, go to the **Editor** tab and click **Run**.  
6. If prompted to *Change Folder* to the current directory, choose **Change Folder** (highlighted option).  

## Code Demo - Program Details
**`Run_Code_Demo.m` generates the yellow line in Fig. 2b, which is $P_{QPP}$ versus time for a 3 micron nanowire. This is a key numerical result in the paper. 
Run_code_demo run time (tested on a Macbook bro): 4 minutes

## Run All Figures - Program Details
**`Run_all_Figures.m`** generates all numerical results and corresponding figures presented in the paper. 

Run_all_figures run time: >72 hours. 

You may also run individual figures by opening the "code" directory, then opening the figure run script of interest in Matlab (for example Run_Figure2b.m). You may then run the code for just that figure in MATLAB by selecting "run". Note that Run_Figurecd.m requires you to first execute Run_FigureS1.m and Run_FigureS2.m requires you to first executre Run_Figure2c.m. 


## Structure of the Code
This code runs our full covariance matrix code and generates all figures in our paper which present numerical results. Run_all_Figures.m generates the figures, in order: Fig. 2c, Fig. 2b, Fig. S1, Fig. 3c, Fig. 3d, Fig. S2, and then Fig. S4. 

## Environment and Dependencies

Requires MATLAB2024a or more recent. 

### Output
All figures are outputted into the folder "results" and are titled according to their corresponding figure and panel in the main text or Supplementary Information of our paper. Additionally, the corresponding MATLAB workspaces generated for each figure are saved in the results folder. 

## Troubleshooting

If you run into any issues or have any questions, please contact m.goffage@unsw.edu.au. 
