%
%     Copyright (C) 2014,2017 Antonio García-Jerez, José Piña-Flores.
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

%% Esta función realiza el acondicionamiento del modelo para su graficación

function [BTA,ALFA,RO,D]=MODELO(aopt,dll)
%INPUT
% aopt information of model
% dll halfsapce depth (only for representation)
%OUTPUT
% BTA, vector for vs values
% ALFA, vector for vp values
% RO, Vector for density values
% D vector for thickness values
% Plots the velocity profile in the GUI Elastic parameters
if ~isnan(aopt)
    %% Inicialización de variables
    d=aopt(:,1)';
    dw=[0 d];
    NC=length(d);
    %pREPARANDO VARIABLES
    [BTA,ALFA,RO,D]=deal(zeros(1,2*NC));
    [dd,m,nn]=deal(zeros(1,NC+1),0,1);
    %Asignando puntos dobles para el formato
    for ii=1:NC
        dd(ii+1)=dw(ii+1)+dd(ii);
    end
    %% Built up profile vectors
    for ij=1:NC*2
        m=m+1;
        [ALFA(ij),BTA(ij),RO(ij),D(ij)]=deal(aopt(nn,2),aopt(nn,3),aopt(nn,4),dd(nn+1));
        if m==2
            [nn,m]=deal(nn+1,0);
        end
    end
    %% Set a "thickness" for the halfspace 
    D=[0 D(1:end-1)];
    D(end)=dll*1.25;
else
    [BTA,ALFA,RO,D]=deal([]);
end