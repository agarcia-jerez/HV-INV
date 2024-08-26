![imagen](https://github.com/user-attachments/assets/8885d9f1-3044-4796-ac4c-a0b6be895db4)

HV-Inv is a computer code for forward calculation and inversion of H/V spectral ratios of ambient noise (HVSRN) based on the diffuse field assumption (DFA). It takes advantage of the recently stated connection between the HVSRN and the elastodynamic Green’s function which arises from the ambient noise interferometry theory. The software supports joint inversion of HVSRN and dispersion curves by using several local and global algorithms: Monte Carlo sampling, simulated annealing, downhill simplex and interior-point.

It has been written in Matlab R2015a

The exe folder contains compiled fortran code for the forward computation of H/V under the diffuse field approach.
The source code for that vile is available at https://github.com/agarcia-jerez/HV-DFA

Note that the definition of H/V used in this software is:

![imagen](https://github.com/user-attachments/assets/b7823b40-af5e-41a8-b6c5-d7068b6dbcf6)

where Ex represents the directional energy density of the component X (North, East or Vertical), that is, adding the horizontal components, not averaging. 

The official website of the HV-Inv project is https://w3.ual.es/GruposInv/hv-inv/
