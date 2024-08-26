%     Copyright (C) 2014,2017 José Piña-Flores, Antonio García-Jerez.
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License version 3 as
%     published by the Free Software Foundation.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GUI for forward calculation of H/V. It calls compiled FORTRAN code HVf.exe
function [HVT,f]=HVPD(m,fmin1,fmax1,n,NR,NL,dk,log,apsv)
%INPUT
%m; Model information
%fmin1 Minimum frequency
%fmax1 Maximun frequency
%n Number of samples
%NR Number of Rayleigh modes 
%NL Number of Love modes
%dk Number of integration points for body waves
%log Sampling Log or lineal
%apsv Parameter of body waves
 % OUTPUT
 %HVT Values amplitude H/V
 %f Vector for frequencies
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Choose \ or / to build paths and executable version depending on OS
if ispc
    barra='\';
    folder='exe_Win';
elseif ismac
    barra='/';
    folder='exe_Mac';
elseif isunix
    barra='/';
    folder='exe_Linux';
end
%% Write model file
mod=['etc',barra,'mod.txt'];
dlmwrite(mod,m(1,1));
dlmwrite(mod, m(2:end,:),'-append','delimiter','\t');
% Create list of frequencies
if log==1
    ff=(logspace(log10(fmin1),log10(fmax1),n).');
else
    ff=(linspace(fmin1,fmax1,n)).';
end
% Write frequency list to file
frec=['etc',barra,'fre.txt'];
dlmwrite(frec, ff,'delimiter','\t','precision',6);
arch=[' -f ',mod,];
%% Please Wait message during forward calulation

%% Run HVf.exe
[~,HVC]=system(['exe',barra,folder,barra,'HVf',' -nmr ',num2str(NR),' -nml ',num2str(NL),' -nks ',num2str(dk),' -apsv ',num2str(apsv),' -ash ','0.001',' -prec 1E-40 ',' -hv ',' -ff ',frec,arch]);
%% Read H/V curves from HVf output
OUT = sscanf(HVC,'%f',[2 inf] )';
[HVT,f]=deal(OUT(:,2),OUT(:,1));
%% Error Message
if isnan(sum(HVT)) || sum(HVT)==0;
    close (mss);
    errordlg('Imposible Calculate H/V','Error H/V');
    [HVT,f]=deal(0);
    return;
end