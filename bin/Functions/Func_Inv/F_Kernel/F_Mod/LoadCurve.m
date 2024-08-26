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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pass]=LoadCurve(handles,flag1)
% Function for loading experimental data in *.txt format (H/V and velocity)
% and plotting these experimental curves
%INPUT
%handles elements for control GUI
% flag1 == 1 for H/V
% flag1 == 2 for Dispersion curve
% flag1 == 3 for Joint Inversion
%OUTPUT
% pass ID for agree operations
%% Initialization of variables
% C.HVFile={};% Name file HV
% C.DCFile={};% Name file dispersion curve
% C.HVFobs=[];% Frequencies of HV
% C.HVAobs=[];% Amplitudes of HV
% C.HVDS=[];% Standard deviation HV
% C.DCFobs=[];% Frequency of dispersion curve
% C.DCVobs=[];% Velocity of dispersion curve
% C.DCDS=[];% Standard deviation of dispersion curve
[C.HVFobs,C.HVAobs,C.HVDS,C.CovInv,C.Cov,C.det,C.DCFobs,C.DCVobs,C.DCDS]=deal([]);
[C.HVFile,C.DCFile]=deal({});
[C.nmHV,C.nmDC]=deal(1);
%Create variables
% pmy=0; % Weight curves for samples
% ww=0; % Weigth introduced by user
ZOMBIE=GetV('ZOMBIE');
LoadV('pmy',0,'ww',0,'stdgDC',3,'stdgHV',3,'flag',flag1);% Carga variables al workspace
[Pl,pass]=deal(struct,true);%PI= identificador de graficas; pass=1 No hay error de lectura en el archivo
%% Open file with HV curve
if flag1 ==1 || flag1==3
    if ZOMBIE.ON
        FileName='zombieHV.txt';
        Path='';
    else
        [FileName, Path]=uigetfile({'*.txt'},'Open HV Curve');
    end
    if FileName~=0
        fid=fopen([Path FileName]);
        HVC=(textscan(fid,' %f %f %f ','delimiter','\t','commentstyle','#'));% Read file  delimited simple tabulador
        fclose(fid);
        if length(HVC{1})<=1
            fid=fopen([Path FileName]);
            HVC=(textscan(fid,' %f %f %f ','delimiter',',','commentstyle','#'));% Read file delimited with simple comma
            fclose(fid);
        end
        if length(HVC{1})<=1
            fid=fopen([Path FileName]);
            HVC=(textscan(fid,' %f %f %f ','delimiter',' ','commentstyle','#'));% Read file delimited with simple space
            fclose(fid);
        end
        
        if isempty(HVC{1}) || length(HVC{1})<=1  % The file is empty or headers uncomment
            warndlg('The file is empty or the headers are not commented with #. Check the file. ','Warning ')
            pass=false;% Cancela la operación
        else
            if  ~isnan (HVC{3}(1,1)) % Si la std no esta vacia
                if length(HVC{1})==length(HVC{2}) && length(HVC{1})==length(HVC{3}) % Vectores del mismo tamaño
                    pass=true;
                else
                    pass=false;% Cuando el archivo no contiene los mismos elementos de frecuencia y amplitud
                    warndlg('Number of frequencies and amplitudes does not  match. Check the file ','Warning')
                end
            else
                if length(HVC{1})==length(HVC{2})
                    pass=true;
                else
                    pass=false;% Cuando el archivo no contiene los mismos elementos de frecuencia y amplitud
                    warndlg('Number of frequencies and amplitudes does not  match. Check the file ','Warning')
                end
            end
        end
        if pass
            if ~isnan(HVC{3}(1,1))
                HVC=cell2mat(HVC);                
                HVC=sortrows(HVC,1);% Sort to get increasing frequencies
                C.HVFobs=HVC(:,1);% frequencies
                C.HVAobs=HVC(:,2);% Amplitudes                
                C.HVDS=HVC(:,3);% Assign std values from file
                stdgHV=2;
            else
                HVC=[HVC{1} HVC{2}];%asigna valores de amplitud y frecuencia.
                HVC=sortrows(HVC,1);% Sort to get increasing frequencies
                C.HVFobs=HVC(:,1);% frequencies
                C.HVAobs=HVC(:,2);% Amplitudes
                %% Getting standard deviation
                if ~ZOMBIE.ON
                    AA=STDG; % call GUI to get the way to compute assumed std
                    uiwait(AA);% waiting GUI
                    [stdg,Vstd]=GetV('stdg','Vstd');%Loading settings made at interface
                    % stdg is the method (1= constant deviation Vstd; 2 = Vstd % of the measurement (frequency-dependent); 3 =(std==1)
                else
                    Vstd=10;stdg=2;% 10% of standard deviation for Zombie mode
                end
                % Computing assumed std 
                if stdg==1
                    C.HVDS=zeros(length(HVC(:,1)),1)+Vstd;
                    stdgHV=1;
                elseif stdg==2
                    C.HVDS=HVC(:,2)*(Vstd/100);
                    stdgHV=2;
                else
                    C.HVDS=ones(length(HVC(:,1)),1);
                    stdgHV=3;
                end
            end           
            C.nmHV=length(HVC(:,1));% number of samples
            pass=true;
            LoadV('stdgHV',stdgHV);
            C.HVFile={FileName};
            % Load Covariance Matrix
            %if exist([Path C.HVFile{1}(1:end-4) 'cov.txt'],'file')
            if exist([Path C.HVFile{1}(1:end-3) 'cov'],'file')            
                fid=fopen([Path C.HVFile{1}(1:end-3) 'cov'],'rt');
                C.Cov=fscanf(fid,'%f',[C.nmHV+1 C.nmHV]);
                C.Cov=C.Cov(2:end,:);% The first row is the list of frequencies
                fclose(fid);
                C.det=det(C.Cov);
                C.CovInv=inv(C.Cov);
            end
        end
    else
        pass=0;
    end
end
%% Open file for dispersion curve
if pass
    if flag1==2 || flag1==3
        if ~ZOMBIE.ON
            %% Getting type of mode dispersion curve
            AA=CDE;
            cancel=GetV('cancel');
            if cancel==1
                delete(AA);
            else
                uiwait(AA)
            end
            cancel=GetV('cancel');
        else
            % Settings for Zombie mode. One dispersion curve only.
            if ZOMBIE.PDC
                %MODE=0 Fundamental mode
                %POL=1 Rayleigh
                %VEL=1 Phase
                LoadV('Path','','FileName','zombiePDC.txt','VEL',1,'POL',1,'MODE',0);            
            else
                %MODE=0 Fundamental mode
                %POL=1 Rayleigh
                %VEL=2 Group
                LoadV('Path','','FileName','zombieGDC.txt','VEL',2,'POL',1,'MODE',0);                            
            end
            cancel=0;
        end
        if cancel==1
            delete(AA);
            pass=0;
        else
            [Path, FileName]=GetV('Path','FileName');
            fid=fopen([Path FileName]);
            CDC=(textscan(fid,' %f %f %f ','delimiter','\t','commentstyle','#'));% Read file  delimited simple tabulador
            fclose(fid);
            if length(CDC{1})<=1
                fid=fopen([Path FileName]);
                CDC=(textscan(fid,' %f %f %f ','delimiter',',','commentstyle','#'));% Read file delimited with simple comma
                fclose(fid);
            end
            if length(CDC{1})<=1
                fid=fopen([Path FileName]);
                CDC=(textscan(fid,' %f %f %f ','delimiter',' ','commentstyle','#'));% Read file delimited with simple space
                fclose(fid);
            end           
            if isempty(CDC{1}) || length(CDC{1})<=1  % The file is empty or headers uncommented
                warndlg('The file is empty or the headers are not commented with #. Check the file. ','Warning ');
                pass=false;
            else
                if  ~isnan (CDC{3}(1,1))% Si la std no esta vacia
                    if length(CDC{1})==length(CDC{2}) && length(CDC{1})==length(CDC{3})
                        pass=true;
                    elseif isnan (CDC{3}(1,1))% cuando std esta vacio
                        if  length(CDC{1})==length(CDC{2})
                            pass=true;
                        else
                            pass=false;% Cuando no estan comentados los headers
                            warndlg('The file is empty or the headers are not commented with #. Check the file. ','Warning ');
                        end
                    else
                        pass=false;% Cuando el archivo no contiene los mismos elementos de frecuencia, velocidad y std
                        warndlg('Number of frequencies and velocities does not  match. Check the file ','Warning')
                    end
                else
                    if length(CDC{1})==length(CDC{2})
                        pass=true;
                    else
                        pass=false;% Cuando el archivo no contiene los mismos elementos de frecuencia velocidad
                        warndlg('Number of frequencies and velocities does not  match. Check the file ','Warning');
                    end
                end
            end
            %% Condicion si el vector contiene el mismo numero de elementos
            if pass
                if ~isnan(CDC{3}(1,1))
                    CDC=cell2mat(CDC);
                    CDC=sortrows(CDC,1);
                    C.DCFobs=CDC(:,1);
                    C.DCVobs=CDC(:,2);
                    C.DCDS=CDC(:,3);
                    stdgDC=2;
                else
                    CDC=[CDC{1} CDC{2}];
                    CDC=sortrows(CDC,1);
                    C.DCFobs=CDC(:,1);
                    C.DCVobs=CDC(:,2);                    
                    %% Getting standard deviation
                    if ~ZOMBIE.ON
                        AA=STDG;
                        uiwait(AA);
                        [stdg,Vstd]=GetV('stdg','Vstd');
                    else
                        Vstd=10;stdg=2;% 10% of standard deviation
                        ww=0.5;% Peso 0.5 a cada curva si inversion conjunta
                        pmy=1;% Normalización por miestras si inversión conjunta
                        LoadV('ww',ww,'pmy',pmy);
                    end
                    if stdg==1% Cuando el usuario ingresa std constante
                        C.DCDS=zeros(length(CDC(:,1)),1)+Vstd;
                        stdgDC=1;
                    elseif stdg==2%Cuando std es proporcional a la amplitud
                        C.DCDS=CDC(:,2)*(Vstd/100);
                        stdgDC=2;
                    else% Cuando el usuario no ingresa std
                        C.DCDS=ones(length(CDC(:,1)),1);
                        stdgDC=3;
                    end
                end
                C.nmDC=length(CDC(:,1));                
                C.DCFile={FileName};
                LoadV('stdgDC',stdgDC);
            end
        end
    end
end

%% Plot HV & dispersion curve in the GUI
if pass
    delete(findobj(handles.PROFILE,'type','axes'))
    delete(findobj(handles.TARGETS,'type','axes'))
    if flag1==1 % Grafica para H/V
        Pl.HV=subplot(1,1,1,'Parent',handles.TARGETS);
        [Pl.HV.NextPlot,Pl.HV.XScale,Pl.HV.FontUnits]=deal('add', 'log','normalized');
    elseif flag1==2%Grafica para la CD
        Pl.DC=subplot(1,1,1,'Parent',handles.TARGETS);
        [Pl.DC.NextPlot,Pl.DC.XScale,Pl.DC.FontUnits]=deal('add', 'log','normalized');
    elseif flag1==3% Grafica para el H/V y CD
        Pl.HV=subplot(1,2,1,'Parent',handles.TARGETS);
        Pl.DC=subplot(1,2,2,'Parent',handles.TARGETS);
        [Pl.DC.NextPlot,Pl.DC.XScale,Pl.DC.FontUnits,...
            Pl.HV.NextPlot,Pl.HV.XScale,Pl.HV.FontUnits]=deal('add', 'log','normalized','add', 'log','normalized');
    end
    % Plot HV curve
    if flag1==1 || flag1==3
        if stdgHV~=3
            p1=errorbar(Pl.HV,C.HVFobs,C.HVAobs,C.HVDS,'-k');
        else
            p1=plot(Pl.HV,C.HVFobs,C.HVAobs,'-k');
        end
        % Control de las graficas (apariencia)
        [Pl.HV.XLim,Pl.HV.YLim,p1.LineWidth,p1.Parent.XScale]=deal([min(C.HVFobs) max(C.HVFobs)],...
            [min(C.HVAobs)-min(C.HVAobs)*0.5 max(C.HVAobs)+max(C.HVAobs)*0.3],2,'Log');
        l1=legend (Pl.HV,C.HVFile);
        [Pl.HV.XLabel.Interpreter,Pl.HV.YLabel.Interpreter,Pl.HV.Title.Interpreter]=deal('latex');
        [Pl.HV.Title.FontSize,l1.FontSize,Pl.DC.Box,Pl.DC.NextPlot,Pl.DC.XScale,l1.Interpreter,...
            Pl.HV.Title.String,Pl.HV.XLabel.String,Pl.HV.YLabel.String]=...
            deal(12,8,'on','add','log','none','$H/V$ $ratio$','$Frequency$ $[Hz]$','$Amplitude$ [${H} \over {V}$]');
        %% GUI menu for plotting HV curve
        c = uicontextmenu;
        uimenu(c,'Label','XScale: Log','Callback',@HV);
        uimenu(c,'Label','XScale: Linear','Callback',@HV);
        uimenu(c,'Label','YScale: Log','Callback',@HV);
        uimenu(c,'Label','YScale: Linear','Callback',@HV);
        set(Pl.HV,'UIContextMenu',c)
    end
    %% Plot Dispersion curve
    if flag1==2 || flag1==3
        if  stdgDC~=3
            p4=errorbar(Pl.DC,C.DCFobs,C.DCVobs,C.DCDS,'-k');
        else
            p4=plot(Pl.DC,C.DCFobs,C.DCVobs,'-k');
        end
        [Pl.DC.XLim,Pl.DC.YLim,p4.LineWidth,p4.Parent.XScale]=deal([min(C.DCFobs) max(C.DCFobs)],...
            [min(C.DCVobs)-min(C.DCVobs)*0.5 max(C.DCVobs)+max(C.DCVobs)*0.3],2,'Log');
        if GetV('VEL')==1
            Pl.DC.YLabel.String='$Phase$ $velocity$ [${m} \over {s}$]';
        else
            Pl.DC.YLabel.String='$Group$ $velocity$ [${m} \over {s}$]';
        end
        if GetV('POL')==1
            TitCurv='$Rayleigh$ $wave$';
        elseif GetV('POL')==2
            TitCurv='$Love$ $wave$';
        end
        l1=legend (Pl.DC,C.DCFile);
        [Pl.DC.XLabel.Interpreter,Pl.DC.YLabel.Interpreter,Pl.DC.Title.Interpreter]=deal('latex');
        [Pl.DC.Title.FontSize,l1.FontSize,Pl.DC.Box,Pl.DC.NextPlot,Pl.DC.XScale,l1.Interpreter,...
            Pl.DC.Title.String,Pl.DC.XLabel.String]=deal(12,8,'on','add','log','none',TitCurv,'$Frequency$ $[Hz]$');
        %% GUI menu of Visualization of dispersion curve
        c = uicontextmenu;
        uimenu(c,'Label','XScale: Log','Callback',@DC);
        uimenu(c,'Label','XScale: Linear','Callback',@DC);
        uimenu(c,'Label','YScale: Log','Callback',@DC);
        uimenu(c,'Label','YScale: Linear','Callback',@DC);
        Pl.DC.UIContextMenu=c;
    end
    LoadV('Pl',Pl,'C',C);
end

% GUI Function for plotting of HV curve
function HV(source,~)
Pl=GetV('Pl');
switch source.Label
    case 'XScale: Log'
        Pl.HV.XScale= 'log';
    case 'XScale: Linear'
        Pl.HV.XScale= 'linear';
    case 'YScale: Log'
        Pl.HV.YScale='log';
    case 'YScale: Linear'
        Pl.HV.YScale='linear';
end
% GUI Function for plotting Dispersion curve
function DC(source,~)
Pl=GetV('Pl');
switch source.Label
    case 'XScale: Log'
        Pl.DC.XScale='log';
    case 'XScale: Linear'
        Pl.DC.XScale= 'linear';
    case 'YScale: Log'
        Pl.DC.YScale='log';
    case 'YScale: Linear'
        Pl.DC.YScale='linear';
end
