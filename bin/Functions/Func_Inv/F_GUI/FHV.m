%     Copyright (C) 2014,2017 José Piña-Flores, Antonio García-Jerez.
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License version 3 as
%     published by the Free Software Foundation.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more deta
function varargout = FHV(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FHV_OpeningFcn, ...
    'gui_OutputFcn',  @FHV_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%#ok<*DEFNU>

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Variables for GUI and  BASE workspace      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FHV_OpeningFcn(hObject, ~, handles, varargin)
global ground_par keep_poi LockAxis change1thickness keep_hv %parametros globales

% GUI default Variables
handles.output = hObject;% Control de los elementos de la GUI
guidata(hObject, handles);% Control de los elementos de la GUI
%Panel HV
[handles.uitable1.Enable,handles.fmin.Enable,handles.fmax.Enable,handles.L.Enable,...
    handles.R.Enable,handles.bwi.Enable,handles.save.Enable,handles.samples.Enable,...
    handles.pushbutton1.Enable,handles.LMHV.Enable,handles.LMDC.Enable,handles.SM.Enable]=deal('off');% Control de los elementos de la GUI
set(handles.uipanel12.Children(:),'Enable','off');% Control de los elementos de la GUI

%Panel DC
set(handles.uipanel2.Children(:),'Enable','off');% Control de los elementos de la GUI
set(handles.Pol.Children(:),'Enable','off');% Control de los elementos de la GUI
set(handles.Velocity.Children(:),'Enable','off');% Control de los elementos de la GUI
[handles.fminDC.Enable,handles.fmaxDC.Enable,handles.ModDC.Enable,...
    handles.SamplesDC.Enable,handles.SaveDC.Enable]=deal('off');% Control de los elementos de la GUI

% Choice HV or DC
[handles.EnableHV.Enable,handles.EnableDC.Enable,handles.Update.Enable]=deal('off');

% Default settings for calling the forward H/V program:
% Frequency,samples,Modes SW and BW
handles.fmin.String=0.2;
handles.fmax.String=10;
handles.samples.String=100;
handles.bwi.String=500;
handles.L.String=10;
handles.R.String=10;

% Default settings for calling the forward H/V program:
% Frequency,samples,Modes for dispersion curve
handles.fminDC.String=0.2;
handles.fmaxDC.String=10;
handles.SamplesDC.String=100;
handles.ModDC.String=0;
% Global variables (Elastic parameter)
[ground_par,change1thickness,LockAxis,keep_poi,keep_hv]=deal(3,0,0,false,false);
% Identificator of graphical panels in GUI
pl=subplot(1,1,1,'Parent',handles.uipanel11);
pl.FontUnits='normalized';
PosSize=pl.FontSize;
PosPlot=pl.Position;

% HV panel
spHV=subplot(2,1,1,'Parent',handles.uipanel11);% subplot for HVs
spHV.NextPlot='add';
ExHVid=errorbar(spHV,nan,nan,nan,'-r','LineWidth',5);% Experimental HV
set(get(ExHVid,'Children'),{'LineWidth'},{5; 0.5})
ThHVBM=plot(spHV,nan,nan,'-.k','LineWidth',3);% Theoretical HV of loaded (back) model
ThHVid=plot(spHV,nan,nan,'-b','LineWidth',3);% Theoretical HV ID
SubHVPlot=spHV.Position;
spHV.Title.String='H/V curve';
% Activate Zoom and Pan tool into HV panel with mouse
spHV.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));

% Profile Panel
spMDL=subplot(1,1,1,'Parent',handles.uipanel14);% subplot for the ground model
spMDL.NextPlot='add';
auxMDLid=plot(spMDL,nan,nan,'.-b','LineWidth',3);% Auxiliary profile
% Activate Zoom, Pan tool, forward computation into model panel with mouse
spMDL.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));

% DC Panel
spDC=subplot(2,1,2,'Parent',handles.uipanel11);% subplot for HVs
spDC.NextPlot='add';
ExDCid=errorbar(spDC,nan,nan,nan,'-r','LineWidth',5);% Experimental HV
set(get(ExDCid,'Children'),{'LineWidth'},{5; 0.5})
ThDCBM=plot(spDC,nan,nan,'-.k','LineWidth',3);% Theoretical HV of loaded (back) model
ThDCid=plot(spDC,nan,nan,'-b','LineWidth',3);% Theoretical HV ID
[spDC.XLabel.String,spDC.YLabel.String]=deal...
    ('Frequency [Hz]','Phase Velocity(${m} \over {s}$)');% Labels
[spDC.XLabel.Interpreter,spDC.YLabel.Interpreter]=deal('latex');
[spDC.FontUnits]='normalized';
SubDCPlot=spDC.Position;
SubDCfont=spDC.FontSize;
% Activate Zoom and Pan tool into DC panel
spDC.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));

% Zooms
ZoomAction=zoom;
ZoomAction.UIContextMenu=uicontextmenu;
ZoomAction.UIContextMenu.Callback=@(hObject,eventdata)FHV('DesactiveTool',hObject,eventdata,guidata(hObject));
ZoomAction.ActionPreCallback=@MarkZoomedMdlAxis;

% Menu context for graphic
cF = uicontextmenu; %contexmenu for HV
cD = uicontextmenu; %contexmenu for DC
cP = uicontextmenu;% contexmenu for model
[ExHVLegend,ThHVBMLegend,ClearHVDataOption,LastLoadedModelFile]=deal(' '); % legends
[ExDCLegend,ThDCBMLegend,ClearDCDataOption]=deal(' '); % legends
[spHV.XLabel.String,spHV.YLabel.String,spMDL.XLabel.String,spMDL.YLabel.String]=deal...
    ('Frequency [Hz]','Amplitude (${H} \over {V}$)','Velocity Vs [m/s]','Depth [m]');% Labels
[spHV.XLabel.Interpreter,spHV.YLabel.Interpreter,spMDL.XLabel.Interpreter,spMDL.YLabel.Interpreter]=deal('latex');
[spHV.FontUnits,spMDL.FontUnits]=deal('normalized');
SubHVfont=spHV.FontSize;
% Load default variables
Lims=[0.1; 10; 0.5; 10]; %Axes limits HV
Limsdc=[0.1; 10; 200; 2000]; % Axes Limits DC;
%%Default settings
LoadV('fminFHV',0.2,'fmaxFHV',10,'nFHV',100,'nksBW',500,'LFHV',10,'RFHV',10,'IsLogF',1,'apsvFHV',0.001,'ShowAuxMDL',0,...
    'bmodel',nan,'limddFHV',nan,'HVFHV',nan,'FFHV',nan,'HVFHVD',nan,'FFHVD',nan,'spHV',spHV,'spDC',spDC,'spMDL',spMDL,...
    'ThHVid',ThHVid,'ThDCid',ThDCid,'ThDCBM',ThDCBM,'ThHVBM',ThHVBM,...
    'auxMDLid',auxMDLid,'ExHVid',ExHVid,'ExDCid',ExDCid,'ExHVLegend',ExHVLegend,'ThHVBMLegend',ThHVBMLegend,'cF',cF,...
    'ClearHVDataOption',ClearHVDataOption,'ZoomAction',ZoomAction,'ZoomOn',0,'Lims',Lims,'Limsdc',Limsdc,...
    'ClearDCDataOption',ClearDCDataOption,'ExDCLegend',ExDCLegend,'ThDCBMLegend',ThDCBMLegend,'LastLoadedModelFile',LastLoadedModelFile);
% Variables DC
LoadV('fminFDC',0.2,'fmaxFDC',10,'nFDC',100,'ModDC',0,'IsLogFDC',1,'VPhase',1,'RayDC',1);
LoadV('PosPlot',PosPlot,'SubHVPlot',SubHVPlot,'SubDCPlot',SubDCPlot,'PosSize',PosSize,...
    'SubHVfont',SubHVfont,'SubDCfont',SubDCfont)
LoadV('DCbm',0,'HVbm',0);
[spHV.UIContextMenu,spMDL.UIContextMenu,spDC.UIContextMenu]=deal(cF,cP,cD);
% Legend for Menu context for HV
uimenu(cF,'Label','XScale: Log','Callback',@setlinestyle);
uimenu(cF,'Label','XScale: Linear','Callback',@setlinestyle);
uimenu(cF,'Label','YScale: Log','Callback',@setlinestyle);
uimenu(cF,'Label','YScale: Linear','Callback',@setlinestyle);
uimenu(cF,'Label','Axes limits','Callback',@setlinestyle);
% Legend for Menu context for DC
uimenu(cD,'Label','XScale: Log','Callback',@setlinestyleDC);
uimenu(cD,'Label','XScale: Linear','Callback',@setlinestyleDC);
uimenu(cD,'Label','YScale: Log','Callback',@setlinestyleDC);
uimenu(cD,'Label','YScale: Linear','Callback',@setlinestyleDC);
uimenu(cD,'Label','Axes limits','Callback',@setlinestyleDC);
% % Legend for Menu context for model
uimenu(cP,'Label','Lock X Axis','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Lock Y Axis','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Keep Poisson ratio','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Change one thickness','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Move one interface','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Vp Profile','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Vs Profile','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Density Profile','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Keep HV ratio','Callback',@HandleMouseMenu);
LoadV('cP',cP,'cD',cD);
% Apariencia del menu context
[cP.Children(3).Checked,cP.Children(5).Checked,cP.Children(4).Separator,...
    cP.Children(6).Separator,cP.Children(1).Separator,cF.Children(3).Separator,cF.Children(1).Separator,...
    cD.Children(3).Separator,cD.Children(1).Separator]=deal('on');
% Visible graphics
spHV.FontSize=PosSize;
spHV.Position=PosPlot;
[spMDL.Visible,spHV.Visible,spDC.Visible]=deal('off');
% Graphic control of the menu GUI
Open = findall(handles.figure1, 'tag', 'Standard.FileOpen');
delete(Open)
Open = findall(handles.figure1, 'tag', 'Exploration.Rotate');
delete(Open)
Open = findall(handles.figure1, 'tag', 'Exploration.Brushing');
delete(Open)
Open = findall(handles.figure1, 'tag',  'DataManager.Linking');
delete(Open)
Open = findall(handles.figure1, 'tag',  'Annotation.InsertLegend');
delete(Open)
Open = findall(handles.figure1, 'tag',  'Annotation.InsertColorbar');
delete(Open)
Open = findall(handles.figure1, 'tag',  'Standard.SaveFigure');
delete(Open)
Open = findall(handles.figure1, 'tag', 'Standard.NewFigure');
delete(Open)
Open = findall(handles.figure1, 'tag', 'Standard.EditPlot');
delete(Open)
movegui(handles.figure1,'center')
% Verscroll=findjobj(handles.uitable1);
% VScbarFHV=Verscroll.get;
% LoadV('VScbarFHV',VScbarFHV);

% auxMDLid.Parent.Parent.Parent.WindowButtonUpFcn=@(hObject,eventdata)FHV('ClickAction2',hObject,eventdata,guidata(hObject));

%% Update curve automatically
function ClickAction2(a,~,handles)
if strcmp(a.CurrentObject.Tag,'impoly vertex');
    CallHV(handles,1);
end

%% Update
function Update_Callback(~,~,handles)
if handles.Update.Value
    handles.figure1.WindowButtonUpFcn=@(hObject,eventdata)FHV('ClickAction2',hObject,eventdata,guidata(hObject));
else
    handles.figure1.WindowButtonUpFcn='';
end

%% GUI Function edit minimum frequency by user
function fmin_Callback(hObject, ~, ~)
fminFHV = str2double(hObject.String); % Obtiene el valor de la frecuencia minima introducido por el usuario
if isnan(fminFHV)
    hObject.String=0.2;
    errordlg('Input must be a number','Error');
elseif fminFHV<=0
    hObject.String=0.2;
    errordlg('The number of minimun frequency must be greater than 0','Error');
end
LoadV('fminFHV',fminFHV);% Asignacion del valor al workspace

%% GUI Function edit maximum frequency by user
function fmax_Callback(hObject, ~, ~)
fminFHV=GetV('fminFHV');%Get fmin variable from workspace
fmaxFHV = str2double(hObject.String);%Get fmax variable from gui
if isnan(fmaxFHV)
    [hObject.String,fmaxFHV]=deal(10);
    errordlg('Input must be a number','Error');
elseif fmaxFHV<=fminFHV
    [hObject.String,fmaxFHV]=deal(10);
    errordlg('The number of maximun frequency must be greater than minimum frequency','Error');
end
LoadV('fmaxFHV',fmaxFHV);% Asignacion del valor de frecuencia maxima al workspace

%% GUI Function edit number of samples
function samples_Callback(hObject, ~, ~)
nF = str2double(hObject.String);% Número de muestras
%Condiciones
if isnan(nF)
    [hObject.String,nF]=deal(100);
    errordlg('Input must be a number','Error');
elseif nF<2 || not(mod(nF,1))==0;
    [hObject.String,nF]=deal(100);
    errordlg('The number of samples must be greater than 2 and integer','Error');
end
LoadV('nFHV',nF);% Load variable

%% GUI sampling options
function uipanel12_SelectionChangeFcn(hObject,~, handles)
if hObject==handles.Log%get variable
    IsLogF=1;
else
    IsLogF=0;
end
LoadV('IsLogF',IsLogF);%Load Variable

%% GUI Function edit Body wave integrals by user
function bwi_Callback(hObject, ~, ~)
nksBW = str2double(hObject.String);% Get variable
%Condiciones
if isnan(nksBW)
    [hObject.String,nksBW]=deal(50);
    errordlg('Input must be a number','Error');
elseif nksBW<0 || not(mod(nksBW,1))==0;
    [hObject.String,nksBW]=deal(50);
    errordlg('The number of integration must be positive and integer','Error');
end
LoadV('nksBW',nksBW);% Load Variable

%% GUI Function edit number of Love modes by user
function L_Callback(hObject, ~, ~)
LF = str2double(hObject.String);% Get Variable
if isnan(LF)
    [hObject.String,LF]=deal(10);
    errordlg('Input must be a number','Error');
elseif LF<0  || not(mod(LF,1))==0;
    [hObject.String,LF]=deal(10);
    errordlg('The number of modes must be greater than 0 and integer','Error');
end
LoadV('LFHV',LF);% Load Variable

%% GUI Function edit number of Rayleigh modes by user
function R_Callback(hObject, ~, ~)
RF = str2double(hObject.String);%Get Variable
%Condiciones
if isnan(RF)
    [hObject.String,RF]=deal(10);
    errordlg('Input must be a number','Error');
elseif RF<0  || not(mod(RF,1))==0;
    [hObject.String,RF]=deal(10);
    errordlg('The number of modes must be greater than 0 and integer','Error');
end
LoadV('RFHV',RF);%Load Variable

%% GUI Function save file of HV data
function save_Callback(~, ~, ~)
[HV,FFHV]=GetV( 'HVFHV', 'FFHV');%Load data from theorical HV
aop=[FFHV HV];
[file,path] = uiputfile('HV.txt','Save file name');%Path
if file~=0
    tex=[path,file];
    fileID = fopen(tex,'w');
    fprintf(fileID,'#Spectral Ratio H/V');
    fprintf(fileID,'\n');
    fprintf(fileID,['#Frequency [Hz]' ' # Amplitude H/V ']);
    fprintf(fileID,'\n');
%    fclose (fileID);
%    dlmwrite(tex, aop, '-append','delimiter', ' ','precision',6);% Save data
    fprintf(fileID,'%.6f %.6f\n',aop');
    fclose (fileID);
end

%% GUI Function save file of DC data
function SaveDC_Callback(~, ~, ~)
[DC,FFHV,nFDC,RayDC,VPhase]=GetV( 'HVFDC', 'FFDC','nFDC','RayDC','VPhase');%Load data from theorical HV
% Number of modes
NModes=length(FFHV)/nFDC;
% Polarization
if RayDC
    pol='# Rayleigh wave dispersion curve';
else
    pol='# Love wave dispersion curve';
end
% Velocity
if VPhase
    Vel='# Frequency [Hz]   # Phase Velocity [m/s]';
else
    Vel='# Frequency [Hz]   # Group Velocity [m/s]';
end
[file,path] = uiputfile('DC.txt','Save file name');%Path
if file~=0
    tex=[path,file];
    fileID = fopen(tex,'w');
    for i=1:NModes
        if i==1            
            fprintf(fileID,pol);
            fprintf(fileID,'\n');
            fprintf(fileID,Vel);
            fprintf(fileID,'\n');
            fprintf(fileID,'#Fundamental Mode');
            fprintf(fileID,'\n');
            %fclose (fileID);
            %dlmwrite(tex, [FFHV(1:nFDC) 1./DC(1:nFDC)], '-append','delimiter', ' ','precision',6);% Save data
            fprintf(fileID,'%.6f %.6f\n',[FFHV(1:nFDC)';1./DC(1:nFDC)']);            
        else
            fprintf(fileID,['# ' num2str(i-1) ' Higher mode']);
            fprintf(fileID,'\n');
            id=find(DC((nFDC*(i-1))+1:nFDC*i)~=0)+nFDC*(i-1);
            %dlmwrite(tex, [FFHV(id) 1./DC(id)], '-append','delimiter', ' ','precision',6);% Save data
            if ~isempty(id),fprintf(fileID,'%.6f %.6f\n',[FFHV(id)';1./DC(id)']);end                        
        end
    end
    fclose (fileID);
end

%% GUI Stabilization factor (slight damping for body wave integrals)
function RFact_Callback(~, ~, ~)
NUM=(inputdlg('Regularization factor','Waves Parameters',1,{'0.01'}));
if ~isempty(NUM)
    apsvFHV=str2double(NUM{1});
else
    apsvFHV=0.01;
end
LoadV('apsvFHV',apsvFHV);% Load Variable

%% GUI Function elastic parameters to default
function edit1_Callback(hObject, ~, handles)
global dataFHV
NCF = str2double(hObject.String);
if isnan(NCF) || isinf(NCF);
    [hObject.String,NCF]=deal(2);
    errordlg('Input must be a number','Error');
elseif NCF<1 || not(mod(NCF,1))==0;
    hObject.String= 2;
    NCF = 2;
    errordlg('The number of layers should be greater or equal to 1 and integer','Error');
end
% Control of GUI
[handles.uipanel14.Children.Visible,handles.uitable1.Enable,handles.SM.Enable,...
    handles.EnableHV.Enable,handles.EnableDC.Enable,handles.LDCD.Enable,...
    handles.Update.Enable,handles.LHVD.Enable]=deal('on');
% [handles.pushbutton1.Enable,handles.fmin.Enable,handles.fmax.Enable,...
%     handles.samples.Enable,handles.bwi.Enable,handles.L.Enable,handles.R.Enable,handles.SM.Enable,...
%     handles.uipanel12.Children(:).Enable]=deal('On');
% handles.EnableHV.Value=1;
handles.LMHV.Value=0;
handles.LMHV.Enable= 'Off';
handles.LMDC.Value=0;
handles.LMDC.Enable= 'Off';
% Values for the table (default)
numeleF=zeros(NCF,5);
numeleF(1:end-1,1)=10;% Thickness
numeleF(1:end-1,2)=500;numeleF(end,2)=1000;% Vp
numeleF(1:end-1,3)=200;numeleF(end,3)=500;% Vs
numeleF(1:end-1,4)=1500;numeleF(end,4)=2500;%Density
numeleF(:,5)=(0.5-(numeleF(:,3)./numeleF(:,2)).^2)./(1-(numeleF(:,3)./numeleF(:,2)).^2);%Poisson
dataFHV=numeleF;% to the global variable
handles.uitable1.Data=num2cell(numeleF);% to the table
handles.uitable1.ColumnEditable=true(1,5);
% Cleaning backmodel
[ThHVBM,ThDCBM,auxMDLid]=GetV('ThHVBM','ThDCBM','auxMDLid');
[ThHVBM.XData,ThHVBM.YData,ThDCBM.XData,ThDCBM.YData]=deal(nan);
if length(auxMDLid.XData)~=1% =1-> properties of impoly (g) do not exist; ~=1 g does exist
    g=GetV('g'); % Load local copy of g from workspace
    delete((g.Parent.Children(1)));% delete impoly without line
end
% Call to visual model edition function
VisualFit_v8(auxMDLid,handles);
g=GetV('g');
for ii=1:2*NCF
    delete(g.Children(ii).UIContextMenu.Children)
end;% borra el menu del botón derecho
delete(g.Children(end).UIContextMenu.Children)
[g.Parent.Children(end).XData,g.Parent.Children(end).YData]=deal(nan(size(g.Parent.Children(end).YData)));
% [g.Parent.Children(end).Marker,g.Parent.Children(end).LineStyle]=deal('none');
LIM=sum(dataFHV(:,1))*1.25;
handles.uipanel14.Children.YLim=[0 sum(dataFHV(:,1))*1.25 ];
LoadV('bmodel',nan,'limddFHV',nan,'NCFHV',NCF,'ZoomOn',0,'LIM',LIM)% Load Variables

%% GUI Function edit elastic parameters by user (Table of parameters)
function uitable1_CellEditCallback(hObject, ~, handles)
global ground_par keep_poi dataFHV pol_backup keep_hv
[dataF, dataFa] = deal(cell2mat(hObject.Data));% Obtaining values from table
% dataFa is the table contents just after user edition. Vs,Vp and Poisson can be mutually incompatible
% dataFa is preserved in this function for comparissons
% dataF is a copy of dataFa for work (to be fited). It is edited to impose Vs-Vp-Poisson compatibility
% dataFVH is the official model (shared variable) it is updated to dataF at the end of the function
[auxMDLid,ShowAuxMDL,g,bmodel,h]=GetV('auxMDLid','ShowAuxMDL','g','bmodel','h');% get model
% Obtain vertexes por ploting the model
[B,A,R,E]=MODELO(bmodel,GetV('limddFHV'));% beta alpha rho thickness
%LockAxisBack=LockAxis;%LockAxis=evalin('base','LockAxis');

changed=find(dataFa(:,4)~=dataFHV(:,4),1);
if ~isempty(changed)
    % User changed density
    if keep_hv
         % Updating all densities
         dataF(:,4)=rround((dataFa(changed,4)/dataFHV(changed,4))*dataFHV(:,4),2);   
    end
else
    changed=find(dataFa(1:end-1,1)~=dataFHV(1:end-1,1),1);
    if ~isempty(changed)
        % User changed thickness
        if keep_hv
            % Updating all thichnesses and velocities
            dataF(:,1:3)=rround((dataFa(changed,1)/dataFHV(changed,1))*dataFHV(:,1:3),2);            
        end    
    else
        changed=find(dataFa(:,5)~=dataFHV(:,5),1);
        if ~isempty(changed),
            % Poisson's ratio changed, update Vp
            if keep_hv
                % We could restore the Poisson ratio and abort here
                % Otherwise, this setting is ignored
                dataF(changed,5)=dataFHV(changed,5); % Restore
            end            
            % Compute new Vp
            dataF(changed,2)=rround(real(dataF(changed,3).*(sqrt((2*(1-dataF(changed,5)))./...
                                                           (1-(2.*dataF(changed,5))))) ),2);
        else % A velocity changed
            if keep_hv %We want to preserve h/v 
                changed=find(dataFa(:,3)~=dataFHV(:,3),1);
                if ~isempty(changed),
                    dataF(:,1:3)=rround((dataFa(changed,3)/dataFHV(changed,3))*dataFHV(:,1:3),2);
                else
                    changed=find(dataFa(:,2)~=dataFHV(:,2),1);
                    if ~isempty(changed),
                        dataF(:,1:3)=rround((dataFa(changed,2)/dataFHV(changed,2))*dataFHV(:,1:3),2);
                    end
                end
            elseif keep_poi % Keep Poisson's ratio           
                changed=find(dataFa(:,3)~=dataFHV(:,3),1);
                if ~isempty(changed),
                    % Vs changed, update Vp
                    dataF(changed,2)=rround(real(dataF(changed,3).*(sqrt((2*(1-dataF(changed,5)))./...
                                                                  (1-(2.*dataF(changed,5))))) ),2);
                else
                    changed=find(dataFa(:,2)~=dataFHV(:,2),1);                        
                    if ~isempty(changed)
                        % Vp changed, update Vs                            
                        dataF(changed,3)=rround(real(dataF(changed,2).*(sqrt((1-(2.*dataF(changed,5)))./...
                                                                      (2*(1-dataF(changed,5))))) ),2);
                    end
                end
            else % We want to preserve the other velocity
                changed=find(dataFa(:,3)~=dataFHV(:,3),1);
                if isempty(changed)
                    changed=find(dataFa(:,2)~=dataFHV(:,2),1);
                end
                dataF(changed,5)=rround((0.5-(dataF(changed,3)./dataF(changed,2)).^2)./...
                                        (1-(dataF(changed,3)./dataF(changed,2)).^2),3);
            end
        end
    end
end
dataF(end,1)=0;
handles.uitable1.Data=num2cell(dataF); % Upload dataF to Table

% Change the graphic model to show the variable changed in the table
if ~isempty(find(dataFHV(:,2)~=dataFa(:,2),1))
    %Since the user changed Vp, select Vp in the graphic model
    %Asignando valores a variables
    [ground_par,auxMDLid.Parent.XLabel.String,auxMDLid.Parent.Parent.Title,auxMDLid.Parent.UIContextMenu.Children(4).Checked,...
        auxMDLid.Parent.UIContextMenu.Children(3).Checked,auxMDLid.Parent.UIContextMenu.Children(2).Checked]=...
        deal(2,'Velocity Vp [m/s]','Vp Profile','on','off','off');
    if ShowAuxMDL
        [g.Parent.Children(end).XData,g.Parent.Children(end).YData]=deal(A,E);
    end
elseif sum(dataFHV(:,3)~=dataFa(:,3))
    %The user changed Vs, select Vs in the graphic model
    ground_par=3;
    auxMDLid.Parent.XLabel.String='Velocity Vs [m/s]';
    auxMDLid.Parent.Parent.Title='Vs Profile';
    auxMDLid.Parent.UIContextMenu.Children(4).Checked='off';
    auxMDLid.Parent.UIContextMenu.Children(3).Checked='on';
    auxMDLid.Parent.UIContextMenu.Children(2).Checked='off';
    if ShowAuxMDL
        g.Parent.Children(end).XData=B;
        g.Parent.Children(end).YData=E;
    end
elseif sum(dataFHV(:,4)~=dataFa(:,4))
    %The user changed density, select density in the graphic model
    [ground_par,auxMDLid.Parent.XLabel.String,auxMDLid.Parent.Parent.Title,auxMDLid.Parent.UIContextMenu.Children(4).Checked,...
        auxMDLid.Parent.UIContextMenu.Children(3).Checked,auxMDLid.Parent.UIContextMenu.Children(2).Checked]=...
        deal(4,'Density [kg/m3]','Density Profile','off','off','on');
    if ShowAuxMDL
        [g.Parent.Children(end).XData,g.Parent.Children(end).YData]=deal(R,E);
    end
end
% Update dataF to dataFHV (the shared variable):
dataFHV=dataF;
%Update blue profile to dataF (i.e. dataFHV):
[B,A,R,E]=MODELO(dataF(:,1:4),GetV('limddFHV'));% Change dataF to a format to draw
E(end)=sum(dataF(:,1))*1.25;
ProfilePolys=[E;A;B;R].';
pol_backup=[ProfilePolys(:,ground_par) ProfilePolys(:,1)];
setConstrainedPosition(h,pol_backup);
%LockAxis=LockAxisBack;% Restore LockAxis status

%% GUI Function for start computing HV curve Button
function pushbutton1_Callback(~, ~, handles)
CallHV(handles,1);

function y=rround(x,n)
if ismac
    y=round(x,n);
else
    y=roundn(x,-n);
end

%% Function for calculate HV and activate Zoom Pan tools by mouse click into Profile
function ClickAction(~,event,handles)
persistent click
if event.Button==1 %Action with the primary button
    if isempty(click),% One click
        click = 1;
        pause(0.5); %Add a delay to distinguish single click from a double click
        if click == 1
            click = [];
        end
    else % Double click
        CallHV(handles,1)% Call Calculate Theoretical HV 
        click = [];
    end
elseif event.Button==2 % Action with the scroll button
    if isempty(click)% One click
        click = 1;
        pause(0.5); %Add a delay to distinguish single click from a double click
        if click == 1
            % Assignin Uicontextmenu Pan Tool
            ZoomAction=GetV('ZoomAction');
            ZoomAction.Direction='in';
            ZoomAction.Enable='on';
            LoadV('ZoomOn',1);
            click = [];
        end
    else % Double click
        if  strcmp(event.Source.Parent.Title,'HV ratio')
            axis(handles.uipanel11.Children(end),'tight');% HV curve axis
            axis(handles.uipanel11.Children(end),'manual');
        else %% Model axis
            axis(handles.uipanel14.Children,'tight');
            LoadV('ZoomOn',0);
            axis(handles.uipanel14.Children,'manual');
        end
        click = [];
    end
end

%% Desactivation Pan and Zoom tool
function DesactiveTool(~,~,~)%Action with the secundary button
zoom off ; %Desactiva el zoom
[ThHVid,ExHVid,ThHVBM,auxMDLid]=GetV('ThHVid','ExHVid','ThHVBM','auxMDLid');
auxMDLid.Parent.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));
ThHVid.Parent.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));
ExHVid.Parent.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));
ThHVBM.Parent.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));

function HandleMouseMenu(source,~)
global keep_hv ground_par keep_poi LockAxis change1thickness dataFHV pol_backup
[h,spMDL,auxMDLid,ShowAuxMDL,bmodel]=GetV('h','spMDL','auxMDLid','ShowAuxMDL','bmodel');
LockAxisBack=LockAxis;
% % Formato para pintar back model
[BTA,ALFA,DEN,E]=MODELO(bmodel,GetV('limddFHV'));
switch get(source,'Label')
    case 'Vs Profile' % Caso Vs
        [source.Parent.Children(2).Checked,source.Parent.Children(4).Checked]=deal('off');
        [source.Parent.Children(3).Checked,spMDL.XLabel.String,spMDL.Parent.Title,ground_par]...
            =deal('on','Velocity Vs [m/s]','Vs Profile',3);
        if ShowAuxMDL % Existencia del back model
            [auxMDLid.XData,auxMDLid.YData(1:end-1)]=deal(BTA,E(1:end-1));
        end
        
    case 'Vp Profile' % Caso Vp
        [source.Parent.Children(4).Checked,spMDL.XLabel.String,spMDL.Parent.Title,ground_par]...
            =deal('on','Velocity Vp [m/s]','Vp Profile',2);% Asignando variables
        [source.Parent.Children(2).Checked,source.Parent.Children(3).Checked]=deal('off');% Apariencia de la GUI
        if ShowAuxMDL% Existencia del back model
            [auxMDLid.XData,auxMDLid.YData(1:end-1)]=deal(ALFA,E(1:end-1));
        end
        axis auto;
    case 'Density Profile' % Caso de la densidad
        [source.Parent.Children(2).Checked,spMDL.XLabel.String,spMDL.Parent.Title,ground_par]...
            =deal('on','Density [kg/m3]','Density Profile',4);
        [source.Parent.Children(3).Checked,source.Parent.Children(4).Checked]=deal('off');
        if ShowAuxMDL% Existencia del back model
            [auxMDLid.XData,auxMDLid.YData(1:end-1)]=deal(DEN,E(1:end-1));
        end
        axis auto;
    case 'Change one thickness' % Caso del espesor
        [source.Parent.Children(6).Checked,source.Parent.Children(5).Checked,change1thickness]...
            =deal('on','off',1);
    case 'Move one interface'% Caso de la interface
        [source.Parent.Children(5).Checked,source.Parent.Children(6).Checked,change1thickness]...
            =deal('on','off',0);
    case 'Keep Poisson ratio' % Caso Del coediciente de Poisson
        if strcmp(source.Parent.Children(7).Checked,'on')
            [source.Parent.Children(7).Checked,keep_poi]=deal('off',false);
        else
            [source.Parent.Children(7).Checked,keep_poi]=deal('on',true);
        end
    case 'Lock X Axis'
        if LockAxisBack==1,
            [source.Parent.Children(9).Checked,LockAxisBack]=deal('off',0);
        elseif LockAxisBack==2,
            [source.Parent.Children(9).Checked,source.Parent.Children(8).Checked, LockAxisBack]...
                =deal('on','off',1);
        else
            [source.Parent.Children(9).Checked, LockAxisBack]=deal('on',1);
        end
    case 'Lock Y Axis'
        if LockAxisBack==2,
            [source.Parent.Children(8).Checked, LockAxisBack]=deal('off',0);
        elseif LockAxisBack==1,
            [source.Parent.Children(8).Checked,source.Parent.Children(9).Checked,LockAxisBack]...
                =deal('on','off',2);
        else
            [source.Parent.Children(8).Checked, LockAxisBack]=deal('on',2);
        end
    case 'Keep HV ratio'
        if strcmp(source.Parent.Children(1).Checked,'off')
            source.Parent.Children(1).Checked='on';
            keep_hv=true;
        else
            source.Parent.Children(1).Checked='off';
            keep_hv=false;
        end
end
% % Formato para pintar back model
[B,A,R,E]=MODELO(dataFHV(:,1:4),GetV('limddFHV'));
% Actualizando grafica del perfil
E(end)=sum(dataFHV(:,1))*1.25;
ProfilePolys=[E;A;B;R].';
pol_backup =[ProfilePolys(:,ground_par) ProfilePolys(:,1) ];
setConstrainedPosition(h,pol_backup);

% Actualizando limites de la gráfica del perfil
[bmodel,limddFHV]=GetV('bmodel','limddFHV');
if length(bmodel)==1,% bmodel no existe
    source.Parent.Parent.CurrentAxes.XLim=[min(dataFHV(:,ground_par))-min(dataFHV(:,ground_par))*0.05...
        max(dataFHV(:,ground_par))+max(dataFHV(:,ground_par))*0.05  ];
    source.Parent.Parent.CurrentAxes.YLim=[0 ,max(limddFHV*1.25, E(end)) ];
else
    source.Parent.Parent.CurrentAxes.XLim=[min(min(bmodel(:,ground_par))-min(bmodel(:,ground_par))*0.05 , min(dataFHV(:,ground_par))-min(dataFHV(:,ground_par))*0.05)...
        max(max(bmodel(:,ground_par))+max(bmodel(:,ground_par))*0.05 , max(dataFHV(:,ground_par))+max(dataFHV(:,ground_par))*0.05)];
    source.Parent.Parent.CurrentAxes.YLim=[0 ,max(limddFHV*1.25, E(end)) ];
end
source.Parent.Parent.CurrentAxes.Children(end).YData(length(E))= max(limddFHV*1.25, E(end));
LoadV('ZoomOn',0);%Zoom desactivado
LockAxis=LockAxisBack;


%% Import model from *.txt
function LM_Callback(~, ~, handles)
global dataFHV
%Leyendo archivo
[FileName, Path]=uigetfile({'*.txt'},'Load model');
if FileName~=0
    fid=fopen([Path FileName]);
    datos=cell2mat(textscan(fid,' %f %f %f %f ','commentstyle','#'));
    fclose(fid);
    NCF=datos(1,1);
    if isnan(NCF)
        errordlg('Input must be a number','Error');
        return
    elseif NCF<1 || not(mod(NCF,1))==0;
        errordlg('The number of layers should be greater or equal to 1 and integer','Error');
        return
    else
        % Ordenando valores del archivo
        datos(:,5)=(0.5-((datos(:,3))./(datos(:,2))).^2)./(1-((datos(:,3))./(datos(:,2))).^2);
        datos(end,1)=0;%Half space
        [handles.edit1.String,handles.uitable1.Data]=deal(NCF,num2cell(datos(2:end,:)));
        dataFHV=datos(2:end,:);
        % Apariencia de la GUI
        [handles.uitable1.Enable,handles.EnableHV.Enable,handles.EnableDC.Enable,...
            handles.LMHV.Enable,handles.LMDC.Enable,handles.uipanel14.Children.Visible,...
            handles.Update.Enable,handles.SM.Enable]=deal('On');
        handles.LMHV.Value=0;
        handles.LMDC.Value=0;
        LoadV('bmodel',datos(2:end,1:4),'NCFHV',NCF,'LastLoadedModelFile',FileName);% Cargando valores al workspace        
    end
    auxMDLid=GetV('auxMDLid');
    if length(auxMDLid.XData)~=1,% El modelo auxiliar estaba inicializado, destruyo poligono
        g=GetV('g');
        delete(g.Parent.Children(1));
    end
    %LLamando a rutina para perfil interactivo
    VisualFit_v8(auxMDLid,handles);
    g=evalin('base','g');
    % borra el menu del botón derecho
    for ii=1:2*NCF
        delete(g.Children(ii).UIContextMenu.Children)
    end
    % Propiedades de linea para el perfil del backmodel
    delete(g.Children(end).UIContextMenu.Children)
    g.Parent.Children(end).LineStyle='-';
    g.Parent.Children(end).Color=[0 0 0];
    LoadV('limddFHV',sum(datos(2:end,1)),'ShowAuxMDL',1);% Cargando variables
end

%% H/V for Loaded model
function LMHV_Callback(~, ~,handles)
[handles.LMHV.Enable,handles.uipanel11.Children(end).Visible]=deal('on');
[ThHVBM,FileName]=GetV('ThHVBM','LastLoadedModelFile');
if handles.LMHV.Value
    [handles.pushbutton1.Enable,handles.save.Enable]=deal('off');
    handles.EnableHV.Value=1;
    LoadV('HVbm',1,'ThHVBMLegend',FileName);
    CallHV(handles,2);% Compute H/V for Background Model (case 2)
    control.Value=1;
    EnableHV_Callback(control,[],handles)
    handles.pushbutton1.Enable='on';
    MakeLegendHV;
else
    [ThHVBM.XData,ThHVBM.YData]=deal(nan);
    MakeLegendHV;
end
LoadV('HVbm',0);

%% DC for Loaded model
function LMDC_Callback(~, ~,handles)
handles.LMDC.Enable='on';
[ThDCBM,FileName]=GetV('ThDCBM','LastLoadedModelFile');
if handles.LMDC.Value
    [handles.pushbutton1.Enable,handles.SaveDC.Enable]=deal('off');
    handles.EnableDC.Value=1;
    control.Value=1;
    EnableDC_Callback(control,[],handles)
    LoadV('DCbm',1,'ThDCBMLegend',FileName);
    CallHV(handles,2)% Compute DC for Background Model (case 2)
    handles.pushbutton1.Enable='on';
    MakeLegendDC;
else
    [ThDCBM.XData,ThDCBM.YData]=deal(nan);
    MakeLegendDC;
end
LoadV('DCbm',0);

%% Export model to *.txt
function SM_Callback(~, ~,~ )
global dataFHV
[file,path] = uiputfile('Model.txt','Save Model');
if file~=0
    tex=[path,file];
    fileID = fopen(tex,'w');
    fprintf(fileID,'#Number of layers included the half-space');
    fprintf(fileID,'\n');
    fprintf(fileID,['#Thickness [m] ' ' #Vp [m/s] ' ' #Vs [m/s] '  ' #Density [kg/m3] ']);
    fprintf(fileID,'\n');
    fprintf(fileID,num2str(size(dataFHV,1)));
    fprintf(fileID,'\n');
    %fclose (fileID);
    %dlmwrite(tex, dataFHV(:,1:4) ,'-append', 'delimiter', ' ','precision','%8.2f');
    fprintf(fileID,'%8.2f %8.2f %8.2f %8.2f\n',dataFHV(:,1:4)');
    fclose (fileID);
end
%% Open HV data curve
function LHVD_Callback(~, ~, handles)
[HVFile, Path]=uigetfile({'*.txt'},'Open HV Curve (Data)');
if HVFile~=0
    fid=fopen([Path HVFile]);
    HVD=(textscan(fid,' %f %f %f ','delimiter','\t','commentstyle','#'));% Read file  delimited simple tabulador
    fclose(fid);
    if length(HVD{1})<=1
        fid=fopen([Path HVFile]);
        HVD=(textscan(fid,' %f %f %f ','delimiter',',','commentstyle','#'));% Read file delimited with simple comma
        fclose(fid);
    end
    if length(HVD{1})<=1
        fid=fopen([Path HVFile]);
        HVD=(textscan(fid,' %f %f %f ','delimiter',' ','commentstyle','#'));% Read file delimited with simple space
        fclose(fid);
    end
    if isempty(HVD{1}) || length(HVD{1})<=1  % The file is empty or headers uncomment
        warndlg('the file is empty or the headers are not comment with ´´#´. Check the file ','!! Warning !!')
        pass=false;% Cancela la operación
    else
        if  ~isnan (HVD{3}(1,1)) % Si la std no esta vacia
            if length(HVD{1})==length(HVD{2}) && length(HVD{1})==length(HVD{3}) % Vectores del mismo tamaño
                pass=true;
            elseif isnan (HVD{3}(1,1)) % cuando std esta vacio
                if  length(HVD{1})==length(HVD{2})
                    pass=true;
                else
                    pass=false;
                    warndlg('the file is empty or the headers are not comment with ´´#´. Check the file ','!! Warning !!')
                end
            else
                pass=false;
                warndlg('The data are not match. Check the file ','!! Warning !!')
            end
        else
            if length(HVD{1})==length(HVD{2})
                pass=true;
            else
                pass=false;
                warndlg('The data are not match. Check the file ','!! Warning !!')
            end
        end
    end
    if pass
        ClearHVDataOption=['Clear ' HVFile];
        [ExHVid,cF]=GetV('ExHVid','cF');
        if length(cF.Children)==6, % Borrar del menu el rotulo de borrado de un posible antiguo modelo
            delete(cF.Children(1));
        end
        uimenu(cF,'Label',ClearHVDataOption,'Callback',@setlinestyle);
        [ExHVid.XData,ExHVid.YData,ExHVid.LData,ExHVid.UData]=deal(nan);
        [ExHVid.XData,ExHVid.YData]=deal(HVD{1},HVD{2});
        if length (HVD{3}) == length (HVD{2})
            [ExHVid.LData,ExHVid.UData]=deal(HVD{3});
        end
        LoadV('ExHVLegend',HVFile,'ClearHVDataOption',ClearHVDataOption,'HVFHVD',HVD{2},'FFHVD',HVD{1})% Load variables workspace
        handles.uipanel11.Children(end).Visible='on';%GUI Changes appearance
        MakeLegendHV;% Leyenda al eje de los HV
        control.Value=1;
        EnableHV_Callback(control,[],handles)
        handles.EnableHV.Value=1;
        if handles.EnableDC.Value
            control.Value=1;
            EnableDC_Callback(control,[],handles)
        end
    end
end

%% Open DC data curve
function LDCD_Callback(~, ~, handles)
[DCFile, Path]=uigetfile({'*.txt'},'Open DC Curve (Data)');
if DCFile~=0
    fid=fopen([Path DCFile]);
    DCD=(textscan(fid,' %f %f %f ','delimiter','\t','commentstyle','#'));% Read file  delimited simple tabulador
    fclose(fid);
    if length(DCD{1})<=1
        fid=fopen([Path HVFile]);
        DCD=(textscan(fid,' %f %f %f ','delimiter',',','commentstyle','#'));% Read file delimited with simple comma
        fclose(fid);
    end
    if length(DCD{1})<=1
        fid=fopen([Path HVFile]);
        DCD=(textscan(fid,' %f %f %f ','delimiter',' ','commentstyle','#'));% Read file delimited with simple space
        fclose(fid);
    end
    if isempty(DCD{1}) || length(DCD{1})<=1  % The file is empty or headers uncomment
        warndlg('the file is empty or the headers are not comment with ´´#´. Check the file ','!! Warning !!')
        pass=false;% Cancela la operación
    else
        if  ~isnan (DCD{3}(1,1)) % Si la std no esta vacia
            if length(DCD{1})==length(DCD{2}) && length(DCD{1})==length(DCD{3}) % Vectores del mismo tamaño
                pass=true;
            elseif isnan (DCD{3}(1,1)) % cuando std esta vacio
                if  length(DCD{1})==length(DCD{2})
                    pass=true;
                else
                    pass=false;
                    warndlg('the file is empty or the headers are not comment with ´´#´. Check the file ','!! Warning !!')
                end
            else
                pass=false;
                warndlg('The data are not match. Check the file ','!! Warning !!')
            end
        else
            if length(DCD{1})==length(DCD{2})
                pass=true;
            else
                pass=false;
                warndlg('The data are not match. Check the file ','!! Warning !!')
            end
        end
    end
    if pass
        %% Select modes
        n=1;
        if length(DCD{1})==length(DCD{3}) 
            dc=zeros(1,3);
            for i=1:length(DCD{1})-1
                if DCD{1}(i)<DCD{1}(i+1)
                    dc(n,:)=[DCD{1}(i) DCD{2}(i) DCD{3}(i)];
                    n=n+1;
                else
                    dc(n,:)=[DCD{1}(i) DCD{2}(i) DCD{3}(i)];
                    n=n+1;
                    dc(n,:)=[nan nan nan];
                    n=n+1;
                end
            end
            dc(n,:)=[DCD{1}(i+1) DCD{2}(i+1) DCD{3}(i+1)];
            DCD{1}=dc(:,1);
            DCD{2}=dc(:,2);
            DCD{3}=dc(:,3);
        else
            dc=zeros(1,2);
            for i=1:length(DCD{1})-1
            if DCD{1}(i)<DCD{1}(i+1)
                dc(n,:)=[DCD{1}(i) DCD{2}(i) ];
                n=n+1;
            else
                dc(n,:)=[DCD{1}(i) DCD{2}(i)];
                n=n+1;
                dc(n,:)=[nan nan];
                n=n+1;
            end
            end
            DCD{1}=dc(:,1);
            DCD{2}=dc(:,2);
        end
        
        ClearDCDataOption=['Clear ' DCFile];
        [ExDCid,cD]=GetV('ExDCid','cD');
        if length(cD.Children)==6, % Borrar del menu el rotulo de borrado de un posible antiguo modelo
            delete(cD.Children(1));
        end
        uimenu(cD,'Label',ClearDCDataOption,'Callback',@setlinestyleDC);
        [ExDCid.XData,ExDCid.YData,ExDCid.LData,ExDCid.UData]=deal(nan);
        [ExDCid.XData,ExDCid.YData]=deal(DCD{1},DCD{2});
        if length (DCD{3}) == length (DCD{2})
            [ExDCid.LData,ExDCid.UData]=deal(DCD{3});
        end
        LoadV('ExDCLegend',DCFile,'ClearDCDataOption',ClearDCDataOption,'HVFDCD',DCD{2},'FFDCD',DCD{1})% Load variables workspace
        handles.uipanel11.Children(1).Visible='on';%GUI Changes appearance
        MakeLegendDC;% Leyenda al eje de los HV
        control.Value=1;
        EnableDC_Callback(control,[],handles)
        handles.EnableDC.Value=1;
        if handles.EnableHV.Value
            control.Value=1;
            EnableHV_Callback(control,[],handles)
        else
            [PosPlot,PosSize,spDC]=GetV('PosPlot','PosSize','spDC');
            spDC.FontSize=PosSize;
            spDC.Position=PosPlot;
        end
    end
end

%% Calculate HV Curve from HVf.exe
function CallHV(handles,idHV)
global dataFHV
mss=msgbox({'Please wait'});
delete(findobj(mss,'string','OK'));
delete(findobj(mss,'style','frame'));
DCbm=GetV('DCbm');
HVbm=GetV('HVbm');
if handles.EnableHV.Value && ~DCbm
    pas=1;% id control
    if idHV==1 % Id theorical HV
        id='ThHVid';
    else
        id='ThHVBM'; % Id back model
        data =GetV('bmodel'); % back model
    end
    %Get variables
    [HVid,apsvFHV,NC,fmin,fmax,R,L,nksBW,n,IsLogF]=...
        GetV(id,'apsvFHV','NCFHV', 'fminFHV','fmaxFHV','RFHV', 'LFHV','nksBW','nFHV','IsLogF');
    [handles.save.Enable,handles.pushbutton1.Enable,handles.uipanel11.Children(end).Visible]=deal('off', 'off','on');
    % Condicionantes Ramas surface waves and body
    if R==0 && L~=0 && nksBW==0
        errordlg('Impossible to obtain the ratio H/V','Error')
    elseif R==0 && L==0 && nksBW==0
        errordlg('Impossible to obtain the ratio H/V','Error')
    else
        % Creacion de modelo *.txt
        if idHV==1
            par=([[NC; dataFHV(:,1)]  [0 ;dataFHV(:,2) ] [0 ;dataFHV(:,3) ] [0 ;dataFHV(:,4) ]]); % Modelo Test HV
        else
            par=([[NC;data(:,1)] [0;data(:,2)] [0;data(:,3)] [0;data(:,4)]]);% backmodel HV
        end
        % Function of execute HVf.exe
        [HVFHV,FFHV]=HVPD(par,fmin,fmax,n,R,L,nksBW,IsLogF,apsvFHV);
        % Condicionantes
        if isnan(sum(HVFHV)) || isinf(sum(HVFHV))  ;
            imal=isnan(HVFHV) | isinf(HVFHV);
            if length(imal)>length(HVFHV)*0.1
                errordlg('Impossible to obtain the H/V ratio','Error')
                pas=0;% id control
            else
                HVFHV(imal)=interp1(A.Fobs(~imal),HVFHV(~imal),A.Fobs(imal),'linear','spline'); % Interpolación
                pas=1; % id control
            end
        end
        % GUI Graphic to curve
        if pas==1% id control
            HVid.XData=FFHV;
            HVid.YData=HVFHV;            
            MakeLegendHV;% Update Legend
            if idHV==1% solo guarda el HV test
                LoadV('HVFHV',HVFHV,'FFHV',FFHV); % Save data HV (frequency, Amplitude)
            end
            % graphical Settings
            if IsLogF
                HVid.Parent.XScale='log';
            else
                HVid.Parent.XScale='linear';
            end
            %HVid.Parent.XLim=[FFHV(1) FFHV(end)]; %set limits
        end
    end
    [handles.save.Enable,handles.pushbutton1.Enable]=deal('on');
elseif  ~handles.EnableHV.Value
    [ThHVid,ThHVBM]=GetV('ThHVid','ThHVBM');
    [ThHVid.XData,ThHVid.YData,ThHVBM.XData,ThHVBM.YData]=deal(nan);
    MakeLegendHV;
end

%% DIspersion curve
if handles.EnableDC.Value && ~HVbm
    pas=1;
    if idHV==1 % Id theoretical DC
        idDC='ThDCid';
    else
        idDC='ThDCBM'; % Id back model
        data =GetV('bmodel'); % back model
    end
    %Get variables
    [DCid,NC,fminDC,fmaxDC,ModDC,n,IsLogFDC,RayDC,VPhase]=...
        GetV(idDC,'NCFHV', 'fminFDC','fmaxFDC','ModDC','nFDC','IsLogFDC','RayDC','VPhase');
    [handles.SaveDC.Enable,handles.pushbutton1.Enable,]=deal('off', 'off');
    % Condicionantes Ramas surface waves and body
    % Creacion de modelo *.txt
    if idHV==1
        par=([[NC; dataFHV(:,1)]  [0 ;dataFHV(:,2) ] [0 ;dataFHV(:,3) ] [0 ;dataFHV(:,4) ]]); % Modelo Test HV
    else
        par=([[NC;data(:,1)] [0;data(:,2)] [0;data(:,3)] [0;data(:,4)]]);% backmodel HV
    end
    % Function of execute HVf.exe
    [HVFDC,FFDC]=DCPD(par,fminDC,fmaxDC,n,ModDC,RayDC,VPhase,IsLogFDC);
%    MakeLegendDC;% Update Legend
    % Condicionantes
    if isnan(sum(HVFDC)) || isinf(sum(HVFDC))  ;
        imal=isnan(HVFDC) | isinf(HVFDC);
        if length(imal)>length(HVFDC)*0.1
            errordlg('Impossible to obtain the dispersion curve','Error')
            pas=0;% id control
        end
    end
    if isempty(HVFDC)
        errordlg('Impossible to obtain the dispersion curve','Error')
        pas=0;% id control
    end
    % GUI Graphic to curve
    if pas==1% id control
        DCid.XData=reshape([reshape(FFDC,n,ModDC+1);nan(1,ModDC+1)],(n+1)*(ModDC+1),1);
        DCid.YData=reshape([reshape(1./HVFDC,n,ModDC+1);nan(1,ModDC+1)],(n+1)*(ModDC+1),1);                
        MakeLegendDC;% Update Legend
        if idHV==1% solo guarda el HV test
            LoadV('HVFDC',HVFDC,'FFDC',FFDC); % Save data HV (frequency, Amplitude)
        end
        % graphical Settings
        if IsLogFDC
            DCid.Parent.XScale='log';
        else
            DCid.Parent.XScale='linear';
        end
    end
    % Legend of dispersion curve axis
    if handles.VPhase.Value
        DCid.Parent.YLabel.String='Phase Velocity(${m} \over {s}$)';
    else
        DCid.Parent.YLabel.String='Group Velocity(${m} \over {s}$)';
    end
    if handles.RayDC.Value
        DCid.Parent.Title.String='Rayleigh wave';
    else
        DCid.Parent.Title.String='Love wave';
    end
    [handles.SaveDC.Enable,handles.pushbutton1.Enable]=deal('on');
elseif  ~handles.EnableDC.Value
    % Clear data
    [ThDCid,ThDCBM,]=GetV('ThDCid','ThDCBM');
    [ThDCid.XData,ThDCid.YData,ThDCBM.XData,ThDCBM.YData]=deal(nan);
    MakeLegendDC;
end
close (mss);

if ~DCbm && ~HVbm
    if ~handles.EnableDC.Value
        handles.LMDC.Value=0;
        %         LMDC_Callback([], [],handles) %
    end
    if ~handles.EnableHV.Value
        handles.LMHV.Value=0;
        %         LMHV_Callback([], [],handles) %
    end
end

%% Function for Legend for HV axes
function MakeLegendHV()
[spHV,ExHVLegend,ThHVBMLegend]=GetV('spHV','ExHVLegend','ThHVBMLegend');
valid_oids=[];
valid_labels={};
if any(~isnan(spHV.Children(2).XData))
    [valid_oids,valid_labels{end+1}]=deal([valid_oids,spHV.Children(2)],ThHVBMLegend);
end
if ~isnan(spHV.Children(3).XData(1))
    [valid_oids,valid_labels{end+1}]=deal([valid_oids,spHV.Children(3)],ExHVLegend);
end
if any(~isnan(spHV.Children(1).XData))
    [valid_oids,valid_labels{end+1}]=deal([valid_oids,spHV.Children(1)],'H/V test model');
end
if ~isempty(valid_oids)
    legend(spHV,valid_oids,valid_labels,'interpreter','none');
    legend(spHV,'show');
else
    legend(spHV,'hide');
end

%% Function for Legend for DC axes
function MakeLegendDC()
[spDC,ExDCLegend,ThDCBMLegend]=GetV('spDC','ExDCLegend','ThDCBMLegend');
valid_oids=[];
valid_labels={};
if any(~isnan(spDC.Children(2).XData))
    [valid_oids,valid_labels{end+1}]=deal([valid_oids,spDC.Children(2)],ThDCBMLegend);
end
if ~isnan(spDC.Children(3).XData(1))
    [valid_oids,valid_labels{end+1}]=deal([valid_oids,spDC.Children(3)],ExDCLegend);
end
if any(~isnan(spDC.Children(1).XData))
    [valid_oids,valid_labels{end+1}]=deal([valid_oids,spDC.Children(1)],'DC test model');
end
if ~isempty(valid_oids)
    legend(spDC,valid_oids,valid_labels,'interpreter','none');
    legend(spDC,'show');
else
    legend(spDC,'hide');
end

%% GUI Function for Visualization of HV subplot
function setlinestyle(source,~)
[spHV,ClearHVDataOption,cF,Lims]=GetV('spHV','ClearHVDataOption','cF','Lims');
switch get(source,'Label')
    case 'XScale: Log'
        spHV.XScale='log';
    case 'XScale: Linear'
        spHV.XScale='linear';
    case 'YScale: Log'
        spHV.YScale='log';
    case 'YScale: Linear'
        spHV.YScale='linear';
    case ClearHVDataOption
        [ExHVid]=GetV('ExHVid');
        [ExHVid.XData,ExHVid.YData,ExHVid.LData,ExHVid.UData]=deal(nan);
        delete(cF.Children(1));
        MakeLegendHV;
        %         if ~isempty(find(~isnan(FFHV),1)),xlim(spHV,[min(FFHV) max(FFHV)]);end
        %         if ~isempty(find(~isnan(HVFHV),1)),ylim(spHV,[min(HVFHV)-min(HVFHV)*0.5 max(HVFHV)+max(HVFHV)*0.1]);end
        LoadV('HVFHVD',nan,'FFHVD',nan,'ExHVLegend','');
    case 'Axes limits'
        NUM=str2double(inputdlg({'X min','X max','Y min' ,'Y max'},...
            'Limits',1,{num2str(Lims(1)),num2str(Lims(2)),num2str(Lims(3)),num2str(Lims(4))}));
        if ~isempty (NUM)
            Lims=NUM;
            spHV.XLim=[Lims(1) Lims(2)];
            spHV.YLim=[Lims(3) Lims(4)];
            LoadV('Lims',Lims)
        else
            spHV.YLimMode='Auto';
            spHV.XLimMode='Auto';
        end
end

%% GUI Function for Visualization of DC subplot
function setlinestyleDC(source,~)
[spDC,ClearDCDataOption,cD,Limsdc]=GetV('spDC','ClearDCDataOption','cD','Limsdc');
switch get(source,'Label')
    case 'XScale: Log'
        spDC.XScale='log';
    case 'XScale: Linear'
        spDC.XScale='linear';
    case 'YScale: Log'
        spDC.YScale='log';
    case 'YScale: Linear'
        spDC.YScale='linear';
    case ClearDCDataOption
        [ExDCid]=GetV('ExDCid');
        [ExDCid.XData,ExDCid.YData,ExDCid.LData,ExDCid.UData]=deal(nan);
        delete(cD.Children(1));
        MakeLegendDC;
        %         if ~isempty(find(~isnan(FFHV),1)),xlim(spHV,[min(FFHV) max(FFHV)]);end
        %         if ~isempty(find(~isnan(HVFHV),1)),ylim(spHV,[min(HVFHV)-min(HVFHV)*0.5 max(HVFHV)+max(HVFHV)*0.1]);end
        LoadV('HVFDCD',nan,'FFDCD',nan,'ExDCLegend','');
    case 'Axes limits'
        NUM=str2double(inputdlg({'X min','X max','Y min' ,'Y max'},...
            'Limits',1,{num2str(Limsdc(1)),num2str(Limsdc(2)),num2str(Limsdc(3)),num2str(Limsdc(4))}));
        if ~isempty (NUM)
            Limsdc=NUM;
            spDC.XLim=[Limsdc(1) Limsdc(2)];
            spDC.YLim=[Limsdc(3) Limsdc(4)];
            LoadV('Limsdc',Limsdc)
        else
            spDC.YLimMode='Auto';
            spDC.XLimMode='Auto';
        end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dispersion curve
%% GUI Function edit minimum frequency by user
function fminDC_Callback(hObject, ~, ~)
fminFDC = str2double(hObject.String); % Obtiene el valor de la frecuencia minima introducido por el usuario
if isnan(fminFDC)
    hObject.String=0.2;
    errordlg('Input must be a number','Error');
elseif fminFDC<=0
    hObject.String=0.2;
    errordlg('The number of minimun frequency must be greater than 0','Error');
end
LoadV('fminFDC',fminFDC);% Asignacion del valor al workspace

%% GUI Function edit maximum frequency by user
function fmaxDC_Callback(hObject, ~, ~)
fminFDC=GetV('fminFDC');%Get fmin variable from workspace
fmaxFDC = str2double(hObject.String);%Get fmax variable from gui
if isnan(fmaxFDC)
    [hObject.String,fmaxFDC]=deal(10);
    errordlg('Input must be a number','Error');
elseif fmaxFDC<=fminFDC
    [hObject.String,fmaxFDC]=deal(10);
    errordlg('The number of maximun frequency must be greater than minimum frequency','Error');
end
LoadV('fmaxFDC',fmaxFDC);% Asignacion del valor de frecuencia maxima al workspace

%% GUI Function edit number of samples
function samplesDC_Callback(hObject, ~, ~)
nFDC = str2double(hObject.String);% Número de muestras
%Condiciones
if isnan(nFDC)
    [hObject.String,nFDC]=deal(100);
    errordlg('Input must be a number','Error');
elseif nFDC<2 || not(mod(nFDC,1))==0;
    [hObject.String,nFDC]=deal(100);
    errordlg('The number of samples must be greater than 2 and integer','Error');
end
LoadV('nFDC',nFDC);% Load variable

%% GUI sampling options
function uipanel2_SelectionChangeFcn(hObject,~, handles)
if hObject==handles.LogDC%get variable
    IsLogFDC=1;
else
    IsLogFDC=0;
end
LoadV('IsLogFDC',IsLogFDC);%Load Variable
% # Modes DC
function ModDC_Callback(hObject, ~, ~)
ModDC = str2double(hObject.String);%Get Variable
%Condiciones
if isnan(ModDC)
    [hObject.String,ModDC]=deal(0);
    errordlg('Input must be a number','Error');
elseif ModDC<0  || not(mod(ModDC,1))==0;
    [hObject.String,ModDC]=deal(0);
    errordlg('The number of modes must be greater than 0 and integer','Error');
end
LoadV('ModDC',ModDC);%Load Variable

function Pol_SelectionChangeFcn(hObject,~, handles)
if hObject==handles.RayDC%get variable
    RayDC=1;
else
    RayDC=0;
end
LoadV('RayDC',RayDC);%Load Variable

function Velocity_SelectionChangeFcn(hObject,~, handles)
if hObject==handles.VPhase%get variable
    VPhase=1;
else
    VPhase=0;
end
LoadV('VPhase',VPhase);%Load Variable
%%

function EnableHV_Callback(control, ~,handles) %
if control.Value==1
    %Panel HV
    spHV=GetV('spHV');
    [handles.fmin.Enable,handles.fmax.Enable,handles.L.Enable,...
        handles.R.Enable,handles.bwi.Enable,handles.save.Enable,handles.samples.Enable,...
        handles.pushbutton1.Enable]=deal('on');% Control de los elementos de la GUI
    set(handles.uipanel12.Children(:),'Enable','on');% Control de los elementos de la GUI
    handles.LHVD.Enable='on';
    if ~isnan(sum(spHV.Children(1).XData))
        handles.save.Enable='on';
    else
        handles.save.Enable='off';
    end
else
    [handles.fmin.Enable,handles.fmax.Enable,handles.L.Enable,...
        handles.R.Enable,handles.bwi.Enable,handles.save.Enable,handles.samples.Enable]=deal('off');% Control de los elementos de la GUI
    set(handles.uipanel12.Children(:),'Enable','off');% Control de los elementos de la GUI
    [PosPlot,spHV,PosSize,spDC]=GetV('PosPlot','spHV','PosSize','spDC');
    spDC.FontSize=PosSize;
    spDC.Position=PosPlot;
    spHV.Visible='off';
    handles.LHVD.Enable='off';
    spDC.Visible='on';
    handles.LDCD.Enable='on';
    if handles.EnableDC.Value==0
        handles.EnableDC.Value=1;
        set(handles.uipanel2.Children(:),'Enable','on');% Control de los elementos de la GUI
        set(handles.Pol.Children(:),'Enable','on');% Control de los elementos de la GUI
        set(handles.Velocity.Children(:),'Enable','on');% Control de los elementos de la GUI
        [handles.fminDC.Enable,handles.fmaxDC.Enable,handles.ModDC.Enable,...
            handles.SamplesDC.Enable]=deal('on');% Control de los elementos de la GUI
        [spDC.Children(1).Visible,spDC.Children(2).Visible,spDC.Children(3).Visible]=deal('on');
        [PosPlot,spHV,PosSize,spDC]=GetV('PosPlot','spHV','PosSize','spDC');
        spDC.FontSize=PosSize;
        spDC.Position=PosPlot;
        spHV.Visible='off';
        handles.LHVD.Enable='off';
        spDC.Visible='on';
        handles.LDCD.Enable='on';
        if ~isnan(sum(spDC.Children(1).XData))
            handles.SaveDC.Enable='on';
        else
            handles.SaveDC.Enable='off';
        end
    end
end
if handles.EnableHV.Value==1 && handles.EnableDC.Value==0
    [PosPlot,spHV,PosSize,spDC]=GetV('PosPlot','spHV','PosSize','spDC');
    spHV.FontSize=PosSize;
    spHV.Position=PosPlot;
    spDC.Visible='off';
    spHV.Visible='on';
    handles.LHVD.Enable='on';
    handles.LDCD.Enable='off';
end
if handles.EnableHV.Value==1 && handles.EnableDC.Value==1
    [SubHVPlot,SubDCPlot,spHV,SubHVfont,SubDCfont,spDC]=...
        GetV('SubHVPlot','SubDCPlot','spHV','SubHVfont','SubDCfont','spDC');
    spHV.FontSize=SubHVfont;
    spHV.Position=SubHVPlot;
    spDC.FontSize=SubDCfont;
    spDC.Position=SubDCPlot;
    spDC.Visible='on';
    spHV.Visible='on';
    handles.LHVD.Enable='on';
    handles.LDCD.Enable='on';
end

function EnableDC_Callback(control, ~,handles) %
if control.Value==1
    %Panel DC
    set(handles.uipanel2.Children(:),'Enable','on');% Control de los elementos de la GUI
    set(handles.Pol.Children(:),'Enable','on');% Control de los elementos de la GUI
    set(handles.Velocity.Children(:),'Enable','on');% Control de los elementos de la GUI
    [handles.fminDC.Enable,handles.fmaxDC.Enable,handles.ModDC.Enable,...
        handles.SamplesDC.Enable]=deal('on');% Control de los elementos de la GUI
    spDC=GetV('spDC');
    [spDC.Children(1).Visible,spDC.Children(2).Visible,spDC.Children(3).Visible]=deal('on');
    %MakeLegendDC;
    handles.LDCD.Enable='on';
    if ~isnan(sum(spDC.Children(1).XData))
        handles.SaveDC.Enable='on';
    else
        handles.SaveDC.Enable='off';
    end
else
    %Panel DC
    set(handles.uipanel2.Children(:),'Enable','off');% Control de los elementos de la GUI
    set(handles.Pol.Children(:),'Enable','off');% Control de los elementos de la GUI
    set(handles.Velocity.Children(:),'Enable','off');% Control de los elementos de la GUI
    [handles.fminDC.Enable,handles.fmaxDC.Enable,handles.ModDC.Enable,...
        handles.SamplesDC.Enable,handles.SaveDC.Enable]=deal('off');% Control de los elementos de la GUI
    [PosPlot,spHV,PosSize,spDC]=GetV('PosPlot','spHV','PosSize','spDC');
    spHV.FontSize=PosSize;
    spHV.Position=PosPlot;
    spDC.Visible='off';
    [spDC.Children(1).Visible,spDC.Children(2).Visible,spDC.Children(3).Visible]=deal('off');
    spHV.Visible='on';
    handles.LHVD.Enable='on';
    handles.LDCD.Enable='off';
    legend(spDC,'hide');
    if handles.EnableHV.Value==0
        handles.EnableHV.Value=1;
        [handles.fmin.Enable,handles.fmax.Enable,handles.L.Enable,...
            handles.R.Enable,handles.bwi.Enable,handles.samples.Enable,...
            handles.pushbutton1.Enable]=deal('on');% Control de los elementos de la GUI
        set(handles.uipanel12.Children(:),'Enable','on');% Control de los elementos de la GUI
        [PosPlot,spHV,PosSize,spDC]=GetV('PosPlot','spHV','PosSize','spDC');
        spHV.FontSize=PosSize;
        spHV.Position=PosPlot;
        spDC.Visible='off';
        [spDC.Children(1).Visible,spDC.Children(2).Visible,spDC.Children(3).Visible]=deal('off');
        spHV.Visible='on';
        handles.LHVD.Enable='on';
        handles.LDCD.Enable='off';
        if ~isnan(sum(spHV.Children(1).XData))
            handles.save.Enable='on';
        else
            handles.save.Enable='off';
        end
    end
end
if handles.EnableDC.Value==1 && handles.EnableHV.Value==0
    [PosPlot,spHV,PosSize,spDC]=GetV('PosPlot','spHV','PosSize','spDC');
    spDC.FontSize=PosSize;
    spDC.Position=PosPlot;
    spHV.Visible='off';
    spDC.Visible='on';
    handles.LHVD.Enable='off';
    handles.LDCD.Enable='on';
    [spDC.Children(1).Visible,spDC.Children(2).Visible,spDC.Children(3).Visible]=deal('on');
end
if handles.EnableHV.Value==1 && handles.EnableDC.Value==1
    [SubHVPlot,SubDCPlot,spHV,SubHVfont,SubDCfont,spDC]=...
        GetV('SubHVPlot','SubDCPlot','spHV','SubHVfont','SubDCfont','spDC');
    spHV.FontSize=SubHVfont;
    spHV.Position=SubHVPlot;
    spDC.FontSize=SubDCfont;
    spDC.Position=SubDCPlot;
    spDC.Visible='on';
    [spDC.Children(1).Visible,spDC.Children(2).Visible,spDC.Children(3).Visible]=deal('on');
    spHV.Visible='on';
    handles.LHVD.Enable='on';
    handles.LDCD.Enable='on';
end

%% GUI close interface option
function figure1_CloseRequestFcn(~, ~, ~)
choice = questdlg('Quitting?','Are you Sure?','Yes','No','No');
switch choice
    case 'Yes'
        closereq
    otherwise
        return
end

function Exit_Callback(~, ~, ~)
choice = questdlg('Quitting?','Are you Sure?','Yes','No','No');
switch choice
    case 'Yes'
        closereq
    otherwise
        return
end

function Menu_Callback(~, ~, ~)
function edit1_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function varargout = FHV_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
function fmin_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function fmax_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function samples_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function bwi_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function L_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function R_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function Auxiliary(hObject)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MarkZoomedMdlAxis(hFig,~)
[spMDL,ZoomAction]=GetV('spMDL','ZoomAction');
%hFig.CurrentAxes
if hFig.CurrentAxes==spMDL 
    if strcmp(ZoomAction.Direction,'in')% && strcmp(ZoomAction.Enable,'on')
        LoadV('ZoomOn',1);
    else
        LoadV('ZoomOn',0);
    end
end
