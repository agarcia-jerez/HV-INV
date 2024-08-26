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
%% GUI for forward calculation of Dispersion curves. It calls compiled FORTRAN code HVf.exe
function [DCT,f]=DCPD(m,fmin1,fmax1,n,Mod,Ray,VPhase,log)
%INPUT
%m: Model information
%fmin1: Minimum frequency
%fmax1: Maximun frequency
%n: Number of samples
%Mod: Number of modes 
%Ray: 0 for Love, 1 for Rayleigh
% VPhase: 0 for Group velocity, 1 for phase velocity 
%log: Logatithmic (1) or lineal (0) sampling 

% OUTPUT
%DCT Velocities
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
%% Run HVf.exe
if Ray==1
    if VPhase==1
        [~,DCC]=system(['exe',barra,folder,barra,'HVf',' -nmr ',num2str(Mod+1),' -nml ',num2str(0),' -prec ',num2str(1E-10) ,' -ph ',' -ff ',frec,arch]);
    else
        [~,DCC]=system(['exe',barra,folder,barra,'HVf',' -nmr ',num2str(Mod+1),' -nml ',num2str(0),' -prec ',num2str(1E-10) ,' -gr ',' -ff ',frec,arch]);
    end
else
    if VPhase==1
        [~,DCC]=system(['exe',barra,folder,barra,'HVf',' -nml ',num2str(Mod+1),' -nmr ',num2str(0),' -prec ',num2str(1E-10) ,' -ph ',' -ff ',frec,arch]);
    else
        [~,DCC]=system(['exe',barra,folder,barra,'HVf',' -nml ',num2str(Mod+1),' -nmr ',num2str(0),' -prec ',num2str(1E-10) ,' -gr ',' -ff ',frec,arch]);
    end
end
%% Read DC curves from HVf output
OUT = sscanf(DCC,'%f',inf);
DCT=OUT(3:end);
f=[];
for i=1:OUT(2);
    f=[f; ff];
end
% Error Message
if isnan(sum(DCT)) || sum(DCT)==0;
    close (mss);
    errordlg('Imposible Calculate DC','Error DC');
    [DCT,f]=deal(0);
    return;
end