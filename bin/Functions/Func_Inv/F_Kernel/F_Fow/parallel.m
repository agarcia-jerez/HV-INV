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

%% Enable parallel processing
%
% EP is the number of workers the user wants to activate
%
% Use EP=1 to ask for the number of available workers
% The answer will be stored in FGP. The workers will not be activated yet
%
% Use EP=1000 to close the workers

% function FGP=parallel(EP)
% if EP==1
%     %% Find out the number of available workers
%     if exist('distcomp','dir')
%         distcomp.feature( 'LocalUseMpiexec', false );
%         c = parcluster;
%         FGP=c.NumWorkers;
%     else
%         FGP=1;
%     end
% elseif EP>=2 && EP<512
%     %% Open EP workers
%     distcomp.feature( 'LocalUseMpiexec', false );
%     set(parcluster(), 'NumWorkers',EP)
%     c=parcluster;
%     parpool(c)
% elseif EP==1000;
%     %% Close workers
%     delete(gcp('nocreate'))
% end
function FGP=parallel(EP)
% EP==1 find out
% FGP == 0 means that the kernels have not been open
% If EP==1, the number of workers (FGP) is found out
% Identifica si existe nucleos en procesadores
if EP==1
    if exist('distcomp','dir')
        distcomp.feature( 'LocalUseMpiexec', false );
        c = parcluster();
        FGP=num2str(c.NumWorkers);
    else
        FGP=1;
    end
elseif EP>=2 && EP<512
    %% Disposición de los nucleos para procesos en paralelo
    distcomp.feature( 'LocalUseMpiexec', false );
    set(parcluster(), 'NumWorkers',EP)
    c=parcluster;
    parpool(c);
    FGP=num2str(c.NumWorkers);
elseif EP==1000;
    %% Close kernels
    delete(gcp('nocreate'))
    FGP=0;
else
    FGP=0;
end