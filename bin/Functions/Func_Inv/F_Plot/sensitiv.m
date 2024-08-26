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

%% Plot the dependence of the misfit on individual parameters during inversion
function sensitiv()
% INPUT
[a,L2]=GetV('MODALL','L2ALL');
%a all models during inversion
%L2 all misfit during inversion


mms=waitbar(0,{'Please wait'});
set(mms, 'CloseRequestFcn','');
% Assign variables
[~, b2]=sort(L2,'descend'); % Misfit
[cmp,mm,L2in,aale,pp] =deal(parula(length(b2)),length(a(:,1,1)),L2(b2),a(:,:,b2),1:4:(4*length(a(:,1,1))));
% Loop for the list of models
fig=figure(10);
set(fig,'name','Parameters space','numbertitle','off');
for i=1:mm
    %% Plot Vs
    figure(10);
    A1=subplot(mm,4,pp(i)+1);
    scatter(A1,aale(i,3,:), L2in, 20, cmp, 'filled');
    [A1.YScale,A1.Box,A1.YLabel.String,A1.XLabel.String,A1.YLabel.Interpreter,A1.XLabel.Interpreter]=...
     deal('log','on','$Misfit$','$Velocity$ $V_s$ [${m} \over {s}$]','latex','latex');
    %% Plot Vp
    A2=subplot(mm,4,pp(i)+2);
    scatter(A2,aale(i,2,:),L2in,20, cmp, 'filled')
      [A2.YScale,A2.Box,A2.YLabel.String,A2.XLabel.String,A2.YLabel.Interpreter,A2.XLabel.Interpreter]=...
     deal('log','on','$Misfit$','$Velocity$ $V_p$ [${m} \over {s}$]','latex','latex');
    %% Plot Density
    A3=subplot(mm,4,pp(i)+3);
    scatter(A3,aale(i,4,:),L2in,20, cmp, 'filled')
    [A3.YScale,A3.Box,A3.YLabel.String,A3.XLabel.String,A3.YLabel.Interpreter,A3.XLabel.Interpreter]=...
    deal('log','on','$Misfit$','$Density$ [${Kg} \over {m^3}$]','latex','latex');
    %% Condition of halfspace
    if i~=mm
        %% Plot Thickness
        A4=subplot(mm,4,pp(i));
        scatter(A4,aale(i,1,:),L2in,20, cmp, 'filled')
        [A4.YScale,A4.Box,A4.YLabel.String,A4.XLabel.String,A4.YLabel.Interpreter,A4.XLabel.Interpreter]=...
        deal('log','on','$Misfit$','$Thickness$ [$m$]','latex','latex');
    end
    waitbar(i/mm);   
end
delete(mms)

