# Sel-Parity-recon-algorithm-paper
MATLAB implementation of two selective parity reconstruction algorithms from the paper Non-CPMG diffusion-weighted turbo spin-echo imaging using the selective parity approach (Aidin Arbabi, Vitaliy Khlebnikov, José P. Marques, and David G Norris). The best way to getting started is to read "readme.docs" file for each of the algorithm.

Folder “matlabStuff” contains the following dependencies for selpar recon:
(1) FID-A-master: scripts for reading Siemens raw “*.dat” files.
Available from https://github.com/CIC-methods/FID-A
(2) SPIRiT_v0.3: scripts for estimating the non-acquired k-space lines with SPIRIT-cg (SPIRIT conjugate gradient).
Available from https://people.eecs.berkeley.edu/~mlustig/Software.html
An alternative is to use SPIRIT-POCS implemented in-house.
(3) selpar_recon: scripts for selpar recon. 
(4) NIfTI read/write: scripts for reading/writing “*.nii” files, used to generate bias field.
