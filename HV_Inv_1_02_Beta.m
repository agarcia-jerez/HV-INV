

%==========================================================================
%                            HV-Inv 1.0.2 Beta (6/24/2016)
%
%HV-Inv is a computer code for forward calculation and inversion of H/V 
%spectral ratios of ambient noise (HVSRN) based on the diffuse field 
%assumption (DFA). 
%
%It takes advantage of the connection between the HVSRN and the elasto-
%dynamic Green's function that arises from the ambient noise interferometry
%theory.
%
%The software support joint inversion of HVSRN and dispersion curves by 
%using several local and global algorithms: Montecarlo sampling, simulated 
%annealing, downhill simplex and interior-point.Forward calculations are 
%performed by HVf (exe folder). Type HVf -h in an OS terminal for help.
%
%Copyright (C) 2014-2016 by José Piña-Flores and Antonio Garcia-Jerez.
%
%Project webpage: http://www.ual.es/GruposInv/hv-inv/
%
%This program is free software: you can redistribute it and/or modify it
%under the terms of the GNU General Public License version 3 as published 
%by the Free Software Foundation.
% 
%This program is distributed in the hope that it will be useful,but WITHOUT
%ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
%FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
%more details.
%
%Refrences
%
%García-Jerez A., Piña-Flores J., Sánchez-Sesma F.J., Luzón F., Perton M.
%(2016) A computer code for forward calculation and inversion of the H/V 
%spectral ratio under the diffuse field assumption, Computers & Geosciences
%,in press.
%
%Piña-Flores J., Perton M., García-Jerez A., Carmona E., Luzón F., Molina-Villegas J.C., 
%Sánchez-Sesma F.J. (2016). The inversion of spectral ratio H/V in a layered system
%using the Diffuse Field Assumption (DFA), Geophysical Journal International, Submitted.
%
%Piña-Flores, J. 2015. Cálculo e inversión del cociente H/V a partir de 
%ruido ambiental. Unpublished M.Sc. Thesis, Universidad Nacional Autónoma 
%de México, México DF, 76 pp. In Spanish.
%
%Sánchez-Sesma, F.J., Rodríguez, M., Iturrarán-Viveros, U., Luzón, F., 
%Campillo, M., Margerin, L., García-Jerez, A., Suarez, M., Santoyo, M.A., 
%Rodríguez-Castellanos, A. (2011) A theory for microtremor H/V spectral 
%ratio: application for a layered medium, Geophysical Journal International
%186, 221-225.
%
%García-Jerez A., Luzón F., Sánchez-Sesma F. J., Lunedei E.,  Albarello D.,
%Santoyo M. A., Almendros J. (2013) Diffuse elastic wavefield within a 
%simple crustal model. Some consequences for low and high frequencies. 
%Journal of Geophysical Research 118(10), 5577-5595.
%
%Funding 
%
%Spanish Ministry of Economy and Competitiveness under grants CGL2014-59908
%and CGL2010-16250 and the European Union with FEDER.
%DGAPA-UNAM under Project IN104712. 
%AXA Research Fund under the 2012 AXA Project: "the use of historical 
%seismic records as realizations of a diffuse field for tomography of 
%alluvial basins: application for the valley of Mexico City".
%==========================================================================
%%
% MATLAB® Version: 8.5.0.197613 (R2015a) & 9.0.0 341360 (R2016a)
%  License Numbers: 125381, 40274058 , 40462396 
%==========================================================================

%% H/V input file format (See /ExampleHV/HV3layer.txt)

% Frequency(Hz) Amplitude-HV HV_std_(optional)
% 0.1 1.3666
% 0.113522 1.36972
% 0.128873 1.37324
% 0.1463 1.37715
% 0.166083 1.38157
% 0.188541 1.38658
% 0.214036 1.3923
% ..


%% DC input file format  DC (See /ExampleDC/DC_Ray_Fun_3layer.txt)
% Frequency(Hz) Velocity(m/s) Velocity_std_(optional) 

% 0.5	1374.474092
% 0.543574	1372.327564
% 0.590946	1369.988752
% 0.642446	1367.442827
% 0.698434	1364.671653
% 0.759301	1361.653919
% 0.825473	1358.367134
% 0.897411	1354.789588
% 0.975619	1350.893143
% 1.06064	1346.654439

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean memory
clear 
close all
clc

%% Script principal for execute GUI HV-INV
addpath(genpath('./etc'));
addpath(genpath('./bin/Functions/Func_HV'));
addpath(genpath('./bin/Functions/Func_Inv/F_GUI'));
addpath(genpath('./bin/Functions/Func_Inv/F_Plot'));
addpath(genpath('./bin/Functions/Func_Inv/F_Kernel/F_Fow'));
addpath(genpath('./bin/Functions/Func_Inv/F_Kernel/F_Inv'));
addpath(genpath('./bin/Functions/Func_Inv/F_Kernel/F_Mod'));

%% Execute HV-Inv 1.0.2 Beta GUI 
HVTI;



