%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        HV-Inv 2.5 Beta GUI
%
% Inversion of the H/V spectral ratio based on the diffuse wavefield theory
%
% Matlab functions deal with GUIs and inversion algorithms. HVTI is the
% main program of HV-Inv 2.5 Beta.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2014,2018 José Piña-Flores, Antonio García-Jerez.
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License version 3 as
%     published by the Free Software Foundation.
%
%     This program is distributed in the hope that it will be useful, but
%     WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function varargout = HVTI(varargin)
%% Function for HV-Inv GUI
% Type of OS
if ispc % Windows
    addpath(genpath('./bin/Functions/Func_Fig/Fig_Win'));
elseif ismac %Mac
    addpath(genpath('./bin/Functions/Func_Fig/Fig_Mac'));
elseif isunix %Linux
    addpath(genpath('./bin/Functions/Func_Fig/Fig_Linux'));
end
gui_Singleton = 1;
gui_State = struct('gui_Name',mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @HVTI_OpeningFcn, ...
    'gui_OutputFcn',  @HVTI_OutputFcn, ...
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
%% Assign default values starting the GUI
function HVTI_OpeningFcn(hObject, ~, handles,  ~)

%% Check for ZOMBIE BEHAVIOR
[ZOMBIE.PARA,ZOMBIE.MDL,ZOMBIE.HV,ZOMBIE.PDC,ZOMBIE.GDC,ZOMBIE.CONF,ZOMBIE.ON,ZOMBIE.ASV]=deal(false);
if exist('zombie.para','file'),ZOMBIE.PARA=true;end
if exist('zombieCONF.m','file'),ZOMBIE.CONF=true;end
if exist('zombieHV.txt','file'),ZOMBIE.HV=true;end
if exist('zombiePDC.txt','file'),ZOMBIE.PDC=true;end
if exist('zombieGDC.txt','file'),ZOMBIE.GDC=true;end
if exist('zombieMDL.txt','file'),ZOMBIE.MDL=true;end
if (ZOMBIE.HV||ZOMBIE.PDC||ZOMBIE.GDC)&&(ZOMBIE.PARA||ZOMBIE.MDL),
    ZOMBIE.ON=true;
    ZOMBIE.ASV=true; % Save tested models and forward calculations in a file, in addition to from the best model.
                     % The file names are results_Profiles.txt,results_HVcurves.txt,results_DCcurves.txt
                     % The best model goes to zombieBEST.txt
                     % This setting can be disabled by means of a ZOMBIE.ASV=false statement in zombieCONF.m
end

if ~ZOMBIE.ON
    img = imread('LogoHVinv.png');
    %% Display logo Java
    jimg = im2java(img);
    frame = javax.swing.JFrame;
    frame.setUndecorated(true);
    icon = javax.swing.ImageIcon(jimg);
    label = javax.swing.JLabel(icon);
    frame.getContentPane.add(label);
    frame.pack;
    imgSize = size(img);
    frame.setSize(imgSize(2),imgSize(1));
    screenSize = get(0,'ScreenSize');  % Get the screen size from the root object
    frame.setLocation((screenSize(3)-imgSize(2))/2,...  % Center on the screen
        (screenSize(4)-imgSize(1))/2);
    frame.show;
    pause(0.5)
    frame.hide;
end

if ismac
    setenv ('DYLD_LIBRARY_PATH','/usr/local/bin:/opt/local/lib:')
end
% GUI Variables
handles.output = hObject;
guidata(hObject, handles);
%% Kernels control
FGP=parallel(1);
if FGP==1
    handles.PARALL.Enable= 'off';
    open=2;% parallel computing not available
else
    LoadV('NUMW',FGP);
    open=0;% Parallel computing is currently dissabled. (=1 means enabled)
end

%% Default settings at workspace
% Identificators of graphical panels in GUI
Perr=[];% Global graphic for misfits
fgd1=[];% Legend of misfit graphic
F1=[];% Graphic for current misfit
F2=[];% Graphic for current misfit

% Default settings for calling the forward H/V program:
kb_by_dk=1000;% Minimum points integration
mkb_by_dk=2000;% Maximum points integration
ramasL=5;% Max. number of Love modes to be computed
ramasR=5;% Max. number of Rayleigh modes to be computed
apsv=0.005;% Stabilization factor (slight damping for body wave integrals); 0 = no dampiing

% General characteristics of allowed models
NC=0;% Number of layers
ZLVs=0;% Low S-wave velocity zones (Vs can decrease as depth increases for some layers); 0 = forbidden; 1 = allowed
ZLVp=1;% Low P-wave velocity zones (Vp can decrease as depth increases for some layers); 0 = forbidden; 1 = allowed
HHS=1;% Hard halfspace (largest Vs) required; 1 = required; 0 = not required
ISDEPTH=false;% Use depths instead of thicknesses
             % Individual model can be either in thickness or depth inside:
             % Optimize,ini1_plus

% Settings for inversion algorithm
Tinv=1;%ID for inversion algorithm. %1=Monte Carlo Sampling (MCS), 2=Simulated Annealing (SA),
       %3=Modified SA (MSA), 4=Simplex ,5=Interior point (IP)
ini1=[];% Initial model
INPO=50;%Default for the number of models in the random initial population:

%% Variables for inversion
% Iterations,% Number of iterations % Number of last iterations,Search radius (for MCS, SA and MSA)
A.Niter=0;% Number of iterations (SA - MSA)
A.Nsubiter=0;% Number temperatures (MSA)
A.Niterf=100;% Number of iterations (MCS) or final temperature interations (SA - MSA)
A.NN=10; % Perturbation range (%) (MCS - SA - MSA)

% Additional settings for SA - MSA inversion algorithm
RT=0.9;%Reduction factor for Temperature
T0=1;% Method for initial Temperature for SA and MSA; T0=1 => use a method based on a probability of increasing the misfit of the initial model.
% If T0=1, the probability of acceptance (in the begining) of a model which multiplies the misfit L2 of the initial model by EA is AP,
% that is, exp(-(EA*L2)/T)=AP or T=-A.EA*A.L2/log(A.AP)
EA=0.1;%Relative misfit increment
AP=0.5;% Acceptance probability

% Overwrite (some) default settings from file zombieCONF.m
if ZOMBIE.CONF
    fid=fopen('zombieCONF.m','rt');    
    while ~feof(fid)
        linea=fgetl(fid);
        eval(linea);
    end
    fclose(fid);
    if Tinv==1,% Monte Carlo Sampling (MC)
        if exist('InvParam','var'),% otherwise we assume default settings 
            A.Niter=0;
            A.Niterf=InvParam(1);
            A.Nsubiter=0;
            A.NN=InvParam(2);            
        end
    elseif Tinv==2,% Simulated Annealing (SA)
        A.Niter=InvParam(1);
        A.Nsubiter=max(InvParam(2)-1,1);        
        A.Niterf=InvParam(3);
        A.NN=InvParam(4);
    elseif Tinv==3,% Modified Simulated Annealing (MSA)
        A.Niter=InvParam(1);
        A.Nsubiter=InvParam(2);
        A.Niterf=InvParam(3);
        A.NN=InvParam(4); 
    elseif Tinv==4 || Tinv==5,% Simplex or Interior Point
        A.Niter=InvParam(1);% Number of iters
        A.Nsubiter=InvParam(2);%MaxFunEvals
        if length(InvParam)>2,A.Niterf=InvParam(3);else A.Niterf=1e-34;end% TolFun
        if length(InvParam)>3,A.NN=InvParam(4);else A.NN=1e-34;end%TolX
    end
end

% If we have changed to DEPTHS (in zombieCONF.m) update header of table.
if ISDEPTH
    LoadV('ISDEPTH',false);% Flag to thickness mode
    ThickOrDepth_Callback([],[], handles);% Press button to set depth mode again, with table header updating
end

% Upload values to workspace
LoadV('Tinv',Tinv,'F2',F2,'F1',F1,'fgd1',fgd1,'Perr',Perr,'ini1',ini1,'A',A,'open',open,'HHS',HHS,'INPO',INPO,...
    'NC',NC,'kb_by_dk',kb_by_dk,'mkb_by_dk',mkb_by_dk,'ramasL',ramasL,'ramasR',ramasR,'ZLVs',ZLVs,'ZLVp',ZLVp,...
    'ISDEPTH',ISDEPTH,'T0',T0,'EA',EA,'RT',RT,'AP',AP,'apsv',apsv,'ZOMBIE',ZOMBIE);

%% Graphic control of the GUI
cla(handles.ploterror,'reset')
[handles.START.Visible,handles.edit1.Enable,handles.PS.Enable,handles.Models.Enable,...
    handles.Pon.Checked,handles.LP.Enable,handles.SPara.Enable,handles.WPara.Enable,...
    handles.ploterror.Visible,handles.LZVS.Checked,handles.MSA.Checked,handles.SA.Checked,...
    handles.MW.Checked,handles.LS.Checked,handles.LF.Checked,handles.LV.Enable,handles.Inversion.Enable,...
    handles.T0.Enable,handles.Redu.Enable,handles.Inversiont.Enable,handles.ThickOrDepth.Enable]=deal('off');
[handles.HHS.Checked,handles.T0misfit.Checked,handles.Poff.Checked,handles.LZVP.Checked]=deal('on');
%     [handles.T0.Enable,handles.Redu.Enable,handles.Inversiont.Enable]=deal('off');

%% GUI Menu Control
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

% Load data in zombie mode
if ZOMBIE.ON
    fprintf(1,'Loading data in zombie mode\n');
    if (ZOMBIE.PDC || ZOMBIE.GDC) && ZOMBIE.HV
        THVDC_Callback([],[], handles);        
    elseif (ZOMBIE.PDC || ZOMBIE.GDC)
        TDC_Callback([], [], handles);
    else
        THV_Callback([],[], handles);
    end
end

% Load .para file (table of ranges) in zombie mode
if ZOMBIE.ON&&ZOMBIE.PARA
    fprintf(1,'Loading parameter ranges in zombie mode\n');
    LP_Callback([],[],handles);
end

% Load initial model (zombieMDL.txt) in zombie mode
if ZOMBIE.ON&&ZOMBIE.MDL
    fprintf(1,'Loading inital model in zombie mode\n');
    LoadV('handles',handles);
    initialmodel;
end

% Start inversion in zombie mode, save best model and close
if ZOMBIE.ON
    figure(handles.figure1);
    START_Callback([],[],handles);
    smodel_Callback([],[],[]);
    if ZOMBIE.ASV
        MODELS;
    end
    % Close program
    if ispc
        delete('etc\*.txt');
    elseif ismac || isunix
        delete('etc/*.txt');
    end
    closereq
end

%% Control of the table of model parameters ranges(GUI)

function edit1_Callback (hObject, ~, handles)
%% Set number of layers (halfspace included)
NC = str2double(hObject.String);
if isnan(NC) || isempty(NC) || isinf(NC)
    [hObject.String,NC]=deal( 2);
    errordlg('Input must be a number','Error');
elseif NC<=1 || not(mod(NC,1))==0
    [hObject.String,NC]=deal( 2);
    errordlg('The number of layers should be greater or equal to 2 and integer','Error');
end
[ini1,INPO,ISDEPTH]=GetV('ini1','INPO','ISDEPTH');
% Check if the number of layers matches that of the initial model
if NC==size(ini1,1) || INPO~=0 % No change in number of layers || no initial model that needs matching with NC
    NCok=1;% OK to put (or reset defaults for) NC layers
else
    choice = questdlg(sprintf('The number of layers does not match the initial model.\n\n Initial model will be deleted'),'Warning','Ok','Cancel','Ok');
    % Handle response
    switch choice
        case 'Ok'
            NCok=1;% OK to put (or reset defaults for) NC layers
        case 'Cancel'
            NCok=0;
            NC=GetV('NC');
            hObject.String=NC;
    end
end
%% Control of fields in the GUI (table of model parameters)
if NCok
    [handles.START.Visible,handles.TableParameters.Enable,...
        handles.SPara.Enable,handles.LV.Enable,handles.T0.Enable,...
        handles.Redu.Enable,handles.Inversiont.Enable,handles.ThickOrDepth.Enable]=deal('on');
    %% Initial ranges for the model parameters
    var=zeros(NC,15);% contains the table
    %default parameters H Vp Vs Density Poisson
    if ISDEPTH
        var(1:end-1,1)=25+25*(0:NC-2)';%Zmin
        var(1:end-1,2)=var(1:end-1,1); %Zmax      
    else        
        var(1:end-1,1)=10;%Hmin
        var(1:end-1,2)=50;%Hmax
    end
    var(:,4)=400;%Vpmin
    var(:,5)=6000; %Vpmax
    var(:,7)=200;% Vsmin
    var(:,8)=3400;%Vsmax
    var(:,10:11)=2000;% min and max density    
    var(:,13)=0.25;%Poisson min
    var(:,14)=0.40;%Poisson max 
    var(:,[1 2 4 5 7 8 10 11])=round(var(:,[1 2 4 5 7 8 10 11]),2);% Some roundings to prevent from invisible figures in the tables
    var(:,3:3:12)=round(var(:,2:3:11)-var(:,1:3:10),2);% Diferences (max-min) rounded
    var(:,13:14)=round(var(:,13:14),4);    
    var(:,15)=round(var(:,14)-var(:,13),4);
    limdd=sum(var(:,2)); %max. depth down to the halfspace
    %% Asignation of variables to workspace
    LoadV('NC',round(NC,0),'var',var,'limdd',limdd,'INPO',50,'ini1',[])
    %Edicion GUI
    handles.TableParameters.Data=num2cell(var(:,[1 2 4 5 7 8 10 11 13 14]));
    % handles.TableParameters.ColumnEditable',true(1,10));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% "Target" Structure %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Target is HV
function THV_Callback(~, ~, handles)
% Reading input file (HV only)
pass=LoadCurve(handles,1);
handles.WPara.Enable= 'on';
LoadV('flag',1)% flag=1 for inversion of H/V only
if pass
    [handles.edit1.Enable,handles.START.Enable,handles.LP.Enable,handles.ThickOrDepth.Enable]=deal('on');
    nm=GetV('C.nmHV');% obteniendo numero de muestras de HV cargado
    if nm>100
        warndlg(sprintf(['The effective number of samples for HV curve is greater than 100 '...
            '(' num2str(nm) ').\n\n',...
            'This is likely to slow down the inversion process.\n\n',...
            'A usual number of samples is 50. Resample the curve if necessary.']),'!! Warning !!')
    end
    handles.Inversion.Enable='on';
    %     [handles.T0.Enable,handles.Redu.Enable,handles.Inversiont.Enable]=deal('off');
end

%% Target DC
function TDC_Callback(~, ~, handles)
% Reading input file (Dispersion curve only)
pass=LoadCurve(handles,2);
handles.WPara.Enable='off';
LoadV('flag',2)% flag=2 for inversion of DC only
if pass
    [handles.edit1.Enable,handles.START.Enable,handles.LP.Enable,handles.ThickOrDepth.Enable]=deal('on');    
    nm=GetV('C.nmDC');
    if nm>100
        warndlg(sprintf(['The effective number of samples for dispersion curve is greater than 100 '...
            '(' num2str(nm) ').\n\n',...
            'This is likely to slow down the inversion process.\n\n',...
            'A usual number of samples is 50. Resample the curve if necessary.']),'!! Warning !!')
    end
    handles.Inversion.Enable='on';
    %     [handles.T0.Enable,handles.Redu.Enable,handles.Inversiont.Enable]=deal('off');
end

%% Target HV+DC
function THVDC_Callback(~,~, handles)
% Reading input file (H/V and dispersion curve)
pass=LoadCurve(handles,3);
LoadV('flag',3) % flag=3 for joint inversion
handles.WPara.Enable='on';
if pass
    [handles.edit1.Enable,handles.START.Enable,handles.LP.Enable,handles.ThickOrDepth.Enable]=deal('on');
    [nmD,nmH]=GetV('C.nmDC','C.nmHV');
    if nmH>100 && nmD>100
        warndlg(sprintf(['The effective number of samples for dispersion curve and DC curve ', '(' num2str(nmD) ') & HV curve ' '(' num2str(nmH) ')' ...
            ' are greater than 100.\n\n',...
            'This is likely to slow down the inversion process.\n\n',...
            'A usual number of samples is 50. Resample all curves if necessary.']),'!! Warning !!')
    elseif nmD>100
        warndlg(sprintf(['The effective number of samples for dispersion curve is greater than 100 '...
            '(' num2str(nmD) ').\n\n',...
            'This is likely to slow down the inversion process.\n\n',...
            'A usual number of samples is 50. Resample the curve if necessary.']),'!! Warning !!')
    elseif nmH>100
        warndlg(sprintf(['The effective number of samples for HV curve is greater than 100 '...
            '(' num2str(nmH) ').\n\n',...
            'This is likely to slow down the inversion process.\n\n',...
            'A usual number of samples is 50. Resample the curve if necessary.']),'!! Warning !!')
    end
    handles.Inversion.Enable='on';
    %     [handles.T0.Enable,handles.Redu.Enable,handles.Inversiont.Enable]=deal('off');
end

%% Foward HV
% New GUI for forward problem
function FHV_Callback(~, ~, ~)
FHV;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Start Inversion%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function START_Callback(~, ~, handles)
[var,varoriginal,ZLVs,ZLVp,NC,HHS,ISDEPTH] =GetV ('var','var','ZLVs','ZLVp','NC','HHS','ISDEPTH');
id=0;%error type in parameters table. 3= Abort; 4=Depths need update
fow=0;%id for updating parameters
[cap1,cap2,cap4]=deal(cell(0));
if ~any(any(var(:,end-2:end-1)>=0.5)), % 0.5 not reached for poisson
    if ~any(any(isnan(var))),% No NaN value
        if any(var(:,2)<varoriginal(:,1))
            if ISDEPTH
                errordlg(sprintf('Some minimum depths are larger than the maximum\n\n Please, change parameters'),'Bad Condition');
            else
                errordlg(sprintf('Some minimum thicknesses are larger than the maximum\n\n Please, change parameters'),'Bad Condition');
            end
            id=3;
        end
        if any(var(:,5)<varoriginal(:,4))
            errordlg(sprintf('Some minimum Vp values are larger than the maximum\n\n Please, change parameters'),'Bad Condition');
            id=3;
        end
        if any(var(:,8)<varoriginal(:,7))
            errordlg(sprintf('Some minimum Vs values are larger than the maximum\n\n Please, change parameters'),'Bad Condition');
            id=3;
        end
        if any(var(:,11)<varoriginal(:,10))
            errordlg(sprintf('Some minimum density values are larger than the maximum\n\n Please, change parameters'),'Bad Condition');
            id=3;
        end
        if any(var(:,14)<varoriginal(:,13))
            errordlg(sprintf('Some minimum Poisson ratio values are larger than the maximum\n\n Please, change parameters'),'Bad Condition');
            id=3;
        end
        if id~=3
            if ISDEPTH
                % Z must increase as the layer index increases
                % Replacing var(kl,2) (max Z) with the minimum among the
                % maximum value of Z in the kl-th layer and the maximun
                % values of the shalower layers
                for kl=NC-2:-1:1
                    var(kl,2)=min(var(kl,2),var(kl+1,2));
                end
                % Replacing var(kl,1) (min Z) with the maximum among the
                % minimum values of Z in the kl-th layer and the minimum
                % values of the shallower layers
                for kl=2:NC-1
                    var(kl,1)=max(var(kl,1),var(kl-1,1));
                end
                var(1:NC-1,3)=var(1:NC-1,2)-var(1:NC-1,1);
                % Check for errors
                if find(var(:,3)<0)
                    listerr=find(var(:,3)<0);
                    for kl=1:length(listerr)
                        cap4{kl}=[' ' num2str(listerr(kl)) ' '];
                    end
                    id=4;
                end
            end
            if ZLVp                
                % Low velocities are allowed in Vp. In this case we perform random model selection in Vs.
                % We find the effective the Vs limits from the Vp range and Poisson ratio range  
%%%                if var(:,14)<0.4999
                    var(:,7)=round_UD(max(var(:,7),var(:,4).*sqrt((1-(2*var(:,14)))./(2*(1-var(:,14))))),2,'U');% We increase minimum Vs if is smaller than the result from maximum poi and minimum Vp
%%%                end
                var(:,8)=round_UD(min(var(:,8),var(:,5).*sqrt((1-(2*var(:,13)))./(2*(1-var(:,13))))),2,'D');% We decrease maximum Vs if is greater than the result from minimum poi and maximum Vp
                var(:,9)=round(var(:,8)-var(:,7),2);
            end
            if ZLVs==0
                % Vs must increase as depth increases
                % Replacing var(kl,8) (max Vs) with the minimum among the
                % maximum value of Vs in the kl-th layer and the maximun
                % values of the shalower layers
                for kl=NC-1:-1:1
                    var(kl,8)=min(var(kl,8),var(kl+1,8));
                end
                % Replacing var(kl,7) (min Vs) with the maximum among the
                % minimum values of Vs in the kl-th layer and the minimum
                % values of the shallower layers
                for kl=2:NC
                    var(kl,7)=max(var(kl,7),var(kl-1,7));
                end
                var(:,9)=round(var(:,8)-var(:,7),2);           
                % Check
                if find(var(:,9)<0)
                    listerr=find(var(:,9)<0);            
                    for kl=1:length(listerr)
                        cap1{kl}=[' ' num2str(listerr(kl)) ' '];
                    end
                    id=1;                    
                end
            else
                if find(var(:,9)<0)
                    listerr=find(var(:,9)<0);
                    for kl=1:length(listerr)
                        cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                    end
                    id=2;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ZLVp==0, % Forbidden Low Velocity Zones in Vp.
                var(:,4)=round_UD(max(var(:,4),var(:,7).*sqrt((2*(1-var(:,13)))./(1-(2*var(:,13))))),2,'U'); % We increase minimum Vp if is smaller than the result from minumum poi and minimum Vs
                var(:,5)=round_UD(min(var(:,5),var(:,8).*sqrt((2*(1-var(:,14)))./(1-(2*var(:,14))))),2,'D'); % We decrease maximum Vp if is greater than the result from maximum poi and maximum Vs
                var(:,6)=round(var(:,5)-var(:,4),2);
                % Replacing var(kl,5) (max Vp) with the minimum among the
                % maximum value of Vp in the kl-th layer and the maximun
                % value in the shalower layers
                for kl=NC-1:-1:1
                    var(kl,5)=min(var(kl,5),var(kl+1,5));
                end
                % Replacing var(kl,4) (min Vp) with the maximum among the
                % minimum values of Vp in the kl-th layer and the minimum
                % value in the shallower layers
                for kl=2:NC
                    var(kl,4)=max(var(kl,4),var(kl-1,4));
                end
                var(:,6)=round(var(:,5)-var(:,4),2);
                if find(var(:,6)<0)
                    listerr=find(var(:,6)<0);
                    for kl=1:length(listerr)
                        cap1{kl}=[' ' num2str(listerr(kl)) ' '];
                    end
                    id=1;
                end
            else
                if find(var(:,6)<0)
                    listerr=find(var(:,6)<0);
                    for kl=1:length(listerr)
                        cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                    end
                    id=2;
                end
            end
        end
        if id==4
            errordlg(['Depth of bottom interface could not be generated for layer(s) ',cap4,' You are coding thickness instead of depth.']);
        end
        if ZLVp~=0
            if id==1
                errordlg(['Vs values could not be generated for layer(s) ',cap1,' Low velocity zones are forbidden.']);
            elseif id==2
                errordlg(['Range of Vp, Vs and Poisson ratio are incompatible for layer(s) ',cap2]);
            end
        elseif ZLVp==0
            if id==1
                errordlg(['Vp values could not be generated for layer(s) ',cap1,' Low velocity zones are forbidden.']);
            elseif id==2
                errordlg(['Range of Vp, Vs and Poisson ratio are incompatible for layer(s) ',cap2]);
            end
        end
        % Ask for update varoriginal -> var
        [~,col]=find(var(:,:)~=varoriginal(:,:));
        if ~isempty(col)  && id==0
%            if ZLVp~=0
%                dia=sprintf('Some values in Vs intervals are incompatible with other constraints!\n\n Upgrade ranges?');
%            else
%                dia=sprintf('Some values in Vp intervals are incompatible with other constraints!\n\n Upgrade ranges?');
%            end
            texto='';
            if any(col==1|col==2)
                texto=[texto,'Some values in Z intervals are incompatible with other constraints!\n'];
            end
            if any(col==4|col==5)
                texto=[texto,'Some values in Vp intervals are incompatible with other constraints!\n'];
            end
            if any(col==7|col==8)
                texto=[texto,'Some values in Vs intervals are incompatible with other constraints!\n'];
            end
            if ~isempty(texto)
                dia=sprintf([texto,'\n Update ranges?']);
                choice = questdlg(dia, ...
                    'Bad condition', ...
                    'Update ranges','Cancel','Update ranges');
                % Handle response
                switch choice
                    case 'Update ranges'
                        handles.TableParameters.Data=num2cell(var(:,[1 2 4 5 7 8 10 11 13 14]));
                        %fow=1;
                        LoadV('var',var);
                    case 'Cancel'
                        warndlg('Change parameters','!! Warning !!')
                        fow=1;
                end
            else
%                 kkk=find(var(:,:)~=varoriginal(:,:))
%                 var(kkk)-varoriginal(kkk)
                handles.TableParameters.Data=num2cell(var(:,[1 2 4 5 7 8 10 11 13 14]));
                LoadV('var',var);                
            end
        end
        [ini1]=GetV('ini1');% Comprobacion del modelo inicial versus parametros de la tabla
        if fow==0 && ~isempty(ini1)
            if ISDEPTH,IniDepths=round(cumsum(ini1(1:end-1,1)),2);end            
            if ~ISDEPTH && any(var(:,1)>ini1(:,1)|var(:,2)<ini1(:,1))
                listerr=find(var(:,1)>ini1(:,1)|var(:,2)<ini1(:,1));
                for kl=1:length(listerr)
                    cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                end
                errordlg(['Thickness of initial model is out of range for layer(s) ',cap2]);
                fow=1;
            elseif ISDEPTH && any(var(1:end-1,1)>IniDepths | var(1:end-1,2)<IniDepths)                
                listerr=find(var(1:end-1,1)>IniDepths|var(1:end-1,2)<IniDepths);
                for kl=1:length(listerr)
                    cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                end
                errordlg(['Depth of initial model is out of range for layer(s) ',cap2]);
                fow=1;                
            elseif any(var(:,4)>ini1(:,2)|var(:,5)<ini1(:,2))
                listerr=find(var(:,4)>ini1(:,2)|var(:,5)<ini1(:,2));
                for kl=1:length(listerr)
                    cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                end
                errordlg(['Vp of initial model is out of range for layer(s) ',cap2]);
                fow=1;
            elseif any(var(:,7)>ini1(:,3)|var(:,8)<ini1(:,3))
                listerr=find(var(:,7)>ini1(:,3)|var(:,8)<ini1(:,3));
                for kl=1:length(listerr)
                    cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                end
                errordlg(['Vs of initial model is out of range for layer(s) ',cap2]);
                fow=1;
            elseif any(var(:,10)>ini1(:,4)|var(:,11)<ini1(:,4))
                listerr=find(var(:,10)>ini1(:,4)|var(:,11)<ini1(:,4));
                for kl=1:length(listerr)
                    cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                end
                errordlg(['Density of initial model is out of range for layer(s) ',cap2]);
                fow=1;
            elseif ZLVp==0 && any(diff(ini1(:,2))<0),
                errordlg(sprintf('Low velocity zones in Vp are forbidden.\n Incompatible initial model.'));
                fow=1;
            elseif ZLVs==0 && any(diff(ini1(:,3))<0),
                errordlg(sprintf('Low velocity zones in Vs are forbidden.\n Incompatible initial model.'));
                fow=1;
            elseif HHS~=0 && any(ini1(end,3)-ini1(:,3)<0),
                errordlg(sprintf('Vs is not maximum for the halfspace in the initial model.\n Check model or settings.'));
                fow=1;
            else
                poi_ini1=(0.5-(ini1(:,3)./ini1(:,2)).^2)./(1-(ini1(:,3)./ini1(:,2)).^2);
                if any(var(:,13)>poi_ini1 | var(:,14)<poi_ini1)
                    listerr=find(var(:,13)>poi_ini1 | var(:,14)<poi_ini1);%Poisson
                    for kl=1:length(listerr)
                        cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                    end
                    errordlg(['Poisson`s ratio of initial model is out of range for layer(s) ',cap2]);
                    fow=1;
                end
            end
        end
        
        if fow==0 && id==0
            % Cleaning memory
            evalin('base','clear(''asens'')');  evalin('base','clear(''ER'')');  evalin('base','clear(''DATAERROR'')');
            evalin('base','clear(''HVALL'')');  evalin('base','clear(''aopt'')');  evalin('base','clear(''MODALL'')');
            evalin('base','clear(''CTALL'')');  evalin('base','clear(''ERROR'')');  evalin('base','clear(''L2ALL'')');
            evalin('base','clear(''L2sens'')');
            % GUI edition
            [handles.edit1.Enable,handles.TDC.Enable,handles.THV.Enable,handles.FHV.Enable,...
                handles.THVDC.Enable,handles.PS.Enable,handles.Models.Enable,handles.Inversion.Enable,...
                handles.LP.Enable,handles.SPara.Enable,handles.smodel.Visible,handles.START.Visible,...
                handles.Report.Enable,handles.Menu.Enable,handles.ThickOrDepth.Enable]=deal('off');
            [handles.ploterror.Visible,handles.STOP.Visible]=deal('on');
            %%Inicializacion de variables
            handles.STOP.UserData=0;
            
            %% Inversion of HV and/or DC
            HVDCINV(handles);
            
            %% GUI edition
            [handles.edit1.Enable,handles.TDC.Enable,handles.THV.Enable,handles.FHV.Enable,...
                handles.THVDC.Enable,handles.PS.Enable,handles.Models.Enable,handles.Inversion.Enable,...
                handles.LP.Enable,handles.SPara.Enable,handles.Menu.Enable,handles.START.Visible,handles.ThickOrDepth.Enable ]=deal('on');
            [handles.STOP.Visible,handles.MSA.Checked,handles.SA.Checked,...
                handles.MW.Checked,handles.LS.Checked,handles.LF.Checked]=deal('off');
            Tinv=GetV('Tinv');
            if Tinv==1
                handles.MW.Checked='on';
            elseif Tinv==2
                handles.SA.Checked='on';
            elseif Tinv==3
                handles.MSA.Checked='on';
            elseif Tinv==4
                handles.LS.Checked='on';
            elseif Tinv==5
                handles.LF.Checked='on';
            end
            if ispc
                delete('etc\*.txt');
            elseif ismac
                delete('etc/*.txt');
            elseif isunix
                delete('etc/*.txt');
            end
        end
    else
        errordlg(sprintf('Some values are NAN\n\n Please, change parameters'),'Bad Condition');
    end
else
    errordlg(sprintf('The Poisson ratio values are larger than 0.5\n\n Please, change parameters'),'Bad Condition');
end

%% Stop Inversion
function STOP_Callback(~, ~, handles)
% Funcion for stopping the program during inversion
choice = questdlg('Stop Inversion?','Finished Inversion','Yes','No','Yes');
switch choice
    case 'Yes'
        %Editing GUI
        handles.STOP.UserData=1;
        handles.STOP.Visible= 'off';
        [handles.smodel.Visible,handles.START.Visible]=deal( 'on');
        msgbox('Finished program. ','H/V Inversion');
    case 'No'
        handles.STOP.UserData=0;
    otherwise
        handles.STOP.UserData=0;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  "Model Parameters" Structure %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TableParameters_CellEditCallback(hObject, ~, handles)
% Rounds model written in the table and save to var
NC=GetV('NC');
var=zeros(NC,15);
var(:,[1 2 4 5 7 8 10 11 13 14]) =cell2mat( hObject.Data);
var(:,[1 2 4 5 7 8 10 11])=round(var(:,[1 2 4 5 7 8 10 11]),2);
var(:,13:14)=round(var(:,13:14),4);
var(:,3:3:12)=round(var(:,2:3:11)-var(:,1:3:10),2);
var(:,15)=round(var(:,14)-var(:,13),4);
limdd=sum(var(:,2));
%Assigning var to the Workspace
LoadV('var',var,'limdd',limdd);
handles.TableParameters.Data=num2cell(var(:,[1 2 4 5 7 8 10 11 13 14]));

function ThickOrDepth_Callback(~,~, handles)
% hObject,eventdata,guidata(hObject)
%% Change between depths and Thickness
ISDEPTH=GetV('ISDEPTH');
ch2=get(handles.uipanel26,'children');% Handle of the Panel
if ISDEPTH
    % Change to Thickness
    ch2(6).ColumnName{1}='Thickness min (m)';
    ch2(6).ColumnName{2}='Thickness max (m)';
    ch2(1).String='Use depths';
else
    % Change to Thickness
    ch2(6).ColumnName{1}=' Depth min  (m) ';
    ch2(6).ColumnName{2}=' Depth max  (m) ';
    ch2(1).String='Use thickness';
end
LoadV('ISDEPTH',~ISDEPTH);

function LP_Callback(hObject, ~, handles)
%% Load model Parameters from *.para
ZOMBIE=GetV('ZOMBIE');
if ZOMBIE.ON&&ZOMBIE.PARA
    FileName='zombie.para';
    Path='';
else
    [FileName, Path]=uigetfile({'*.para'},'Open Parameters');
end
if (FileName)~=0
    fid=fopen([Path FileName]);
    datos=textscan(fid,' %f %f %f %f %f %f %f %f %f %f ','commentstyle','#');
    fclose(fid);
    NC=datos{1}(1);
    if isnan(NC)
        %hObject.String=3;
        %NC = str2double(hObject.String');
%       NC=3;
%       handles.edit1.String=NC;        
        errordlg('Input must be a number','Error');
return        
    elseif NC<2 || not(mod(NC,1))==0;
        %hObject.String= 3;
        %NC = str2double(hObject.String);
%       [handles.TableParameters.Enable,handles.START.Visible]=deal('on');
%        NC=3;
%       handles.edit1.String=NC;        
        errordlg('The number of layers should be equal or greater than 2 and integer.','Error');
return
    else
        [handles.TableParameters.Enable,handles.START.Visible]=deal('on');
        handles.edit1.String=NC;
    end
    var=zeros(NC,15);
    for i=1:2:9% Minima to var
        for j=1:NC
            var(j,(3*i-1)/2)= datos{:,i}(j+1);
        end
    end
    for i=2:2:10% Maxima to var
        for j=1:NC
            var(j,(3*i-2)/2)= datos{:,i}(j+1);
        end
    end
    var(:,[1 2 4 5 7 8 10 11])=round(var(:,[1 2 4 5 7 8 10 11]),2);
    var(:,13:14)=round(var(:,13:14),4);
    var(:,3:3:12)=round(var(:,2:3:11)-var(:,1:3:10),2);
    var(:,15)=round(var(:,14)-var(:,13),4);    
    limdd=sum(var(:,2));
    handles.TableParameters.Data=num2cell(var(:,[1 2 4 5 7 8 10 11 13 14]));
    LoadV('limdd',limdd,'var',var,'NC',NC);
    [handles.START.Visible,handles.TableParameters.Enable,...
        handles.SPara.Enable,handles.LV.Enable,handles.T0.Enable,...
        handles.Redu.Enable,handles.Inversiont.Enable,handles.edit1.Enable,handles.ThickOrDepth.Enable]=deal('on');
end

function SPara_Callback(~, ~, ~)
%% Export model parameters to *.para
[var,NC]=GetV('var','NC');
[file,path] = uiputfile('Parameters.para','Save file');
if file~=0
    ISDEPTH=GetV('ISDEPTH');          
        tex=[path,file];
        fileID = fopen(tex,'w');
        fprintf(fileID,'#Number of layers included the half-space');
        fprintf(fileID,'\n');
        if ISDEPTH
        fprintf(fileID,['#Depth min [m] ' ' #Depth max [m] ' ' #Vp min [m/s] ' ' #Vp max [m/s] ' ...
            ' #Vs min [m/s] ' ' #Vs max [m/s] ' ' #Density min [kg/m3] ' ' #Density max [kg/m3] ' ...
            ' #Poisson min ' ' #Poisson max']);
        else
          fprintf(fileID,['#Thickness min [m] ' ' #Thickness max [m] ' ' #Vp min [m/s] ' ' #Vp max [m/s] ' ...
            ' #Vs min [m/s] ' ' #Vs max [m/s] ' ' #Density min [kg/m3] ' ' #Density max [kg/m3] ' ...
            ' #Poisson min ' ' #Poisson max ']);
        end
        fprintf(fileID,'\n');
%        fclose (fileID);
%        dlmwrite(tex,NC,'-append', 'delimiter', ' ','precision','%d');
%        dlmwrite(tex, var(:,[1 2 4 5 7 8 10 11 13 14]) ,'-append', 'delimiter', ' ','precision','%8.2f');
        fprintf(fileID,'%d\n',NC);
        fprintf(fileID,'%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n',...
                        var(:,[1 2 4 5 7 8 10 11 13 14])');
        fclose (fileID);
end

function NIP_Callback(~, ~, handles)
%% Ask for the nummer of models in the initial population
INPO=GetV('INPO');
if INPO==0
    defaultanswer={num2str(50)};
else
    defaultanswer={num2str(INPO)};
end
NUM=str2double((inputdlg({sprintf('Number of models in initial population(>0)\n\n(Introduce 0 to set an initial model)\n')},...
    ' ',1,defaultanswer)));
if  ~isempty(NUM)
    if ~isnan (NUM);
        LoadV('ini1',[],'INPO',NUM);
    end
else
    LoadV('INPO',INPO);
end
if NUM==0
    LoadV('handles',handles);
    initialmodel;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Funtions for the  Menu %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function PS_Callback(~, ~, ~)
%% Call function which assess sensitivity to parameters
sensitiv;% Uses MODALL and L2ALL

function sal_Callback(~, ~, ~)
%% "Exit" (quit program)
choice = questdlg('Quitting?','Are you Sure?','Yes','No','No');
switch choice
    case 'Yes'
        delete('etc\*.txt');
        closereq
    otherwise
        return
end

function Models_Callback(~, ~, ~)
%% New GUI in "Models" menu. Function to show best fitting models
MODELS

%% Two Functions in "Settings" Menu that enable kernels for parallel inversion
function Poff_Callback(~, ~, handles)
open=GetV('open');
if open==1
    handles.Poff.Checked='on';
    handles.Pon.Checked='off';
    mss=msgbox({'Parallelization Settings','Please Wait'},'H/V');
    delete(findobj(mss,'string','OK'));
    delete(findobj(mss,'style','frame'));
    [~]=parallel(1000);
    close(mss);
    LoadV('open',0);
    handles.textNumwork.Visible='off';
end

function Pon_Callback(~, ~, handles)
[open,pop]=GetV('open','NUMW');
if open==1
    mss=msgbox({'Parallelization Settings','Please Wait'},'H/V');
    delete(findobj(mss,'string','OK'));
    delete(findobj(mss,'style','frame'));
    delete(gcp('nocreate'))
    close(mss);
end
NUMW=str2double((inputdlg({'Enter number of Workers (Labs): (>=2 & <512)'},'',1,{num2str(pop)})));
if ~isnan(NUMW) || NUMW>=2
    handles.Pon.Checked='on';
    [handles.Poff.Checked,handles.edit1.Enable,handles.TDC.Enable,handles.THV.Enable,...
        handles.FHV.Enable,handles.THVDC.Enable,handles.PS.Enable,handles.Models.Enable,...
        handles.Inversion.Enable,handles.LP.Enable,handles.SPara.Enable,handles.START.Enable,...
        handles.STOP.Enable,handles.Menu.Enable]=deal('off');
    mss=msgbox({sprintf('The Parallelization Settings will take a few seconds.\n\nPlease Wait')},'HV-inv');
    delete(findobj(mss,'string','OK'));
    delete(findobj(mss,'style','frame'));
    set(mss, 'CloseRequestFcn','')
    parallel(NUMW);
    delete(mss);
    handles.textNumwork.String=['Connection with',' ',num2str(NUMW),' ', 'Labs.'];
    [handles.edit1.Enable,handles.textNumwork.Visible,handles.TDC.Enable,handles.THV.Enable,...
     handles.FHV.Enable,handles.THVDC.Enable,handles.PS.Enable,handles.Models.Enable,...
     handles.Inversion.Enable,handles.LP.Enable,handles.SPara.Enable,handles.START.Enable,...
     handles.STOP.Enable,handles.Menu.Enable]=deal('on');
    LoadV('NUMW',NUMW,'open',1);
end

function T0free_Callback(~, ~,  handles)
%% Function for user defined "Initial Temperature" menu
AP=str2double((inputdlg({'Temperature'},'T0',1,{'10'})));
LoadV('AP',AP,'T0',2)
handles.T0free.Checked='on';
handles.T0misfit.Checked='off';

function T0misfit_Callback(~, ~,handles)
%% Function for "Initial Temperature" menu (alternative input method)
NUM=str2double(inputdlg({'Relative misfit increment (>0)',' Probability of acceptance (<1)'},'T0',1,{'0.1','0.5'}));
if ~isempty(NUM)
    set(handles.T0free,'Checked','off')
    set(handles.T0misfit,'Checked','on')
    LoadV('T0',1,'EA',NUM(1),'AP',NUM(2));
end

function Redu_Callback(~, ~, ~)
%% Function for "Cooling Schedule" menu
RT=str2double(inputdlg({'Temperature ratio (Ti+1/Ti)'},'Cooling Schedule',1,{'0.9'}));
LoadV('RT',RT);

%%Functions for "Algorithm details" menu

function MW_Callback(~, ~, handles)
% Metropoli Walk (Monte Carlo Sampling)
NUM=str2double((inputdlg({'Number of iterations','Perturbation range (%)'},' ',1,{'100','10'})));
if ~isempty(NUM)
    A=GetV('A');
    handles.MW.Checked=('on');
    [handles.MSA.Checked,handles.SA.Checked,handles.LS.Checked,...
        handles.LF.Checked,handles.Redu.Enable,handles.T0.Enable]=deal('off');
    A.Nsubiter=0;
    A.Niter=0;
    A.Niterf=NUM(1);
    A.NN=NUM(2);
    LoadV('A',A,'Tinv',1);
end

function SA_Callback(~, ~, handles)
%% Function for Simulated Annealing (SA) menu
NUM=str2double(inputdlg({['M',char(225),'rkov`s chain length'],'Number of temperatures',...
    ['Last M',char(225),'rkov`s chain length'],'Perturbation range (%)'},' ',1,{'100','1','100','10'}));
if ~isempty(NUM)
    A=GetV('A');
    [handles.SA.Checked,handles.Redu.Enable,handles.T0.Enable]=deal('on');
    [handles.LS.Checked,handles.LF.Checked,handles.MW.Checked,handles.MSA.Checked]=deal('off');
    A.Nsubiter=max(NUM(2),1);  
    A.Niter=NUM(1);
    A.Niterf=NUM(3);
    A.NN=NUM(4);
    LoadV('A',A,'Tinv',2);
end


function MSA_Callback(~, ~, handles)
%% Function for Modified Simulated Annealing (MSA) menu
NUM=str2double((inputdlg({'Number of iterations','Number of Reheatings',...
    'Number of Last iterations','Perturbation range (%)'},' ',1,{'100','0','100','10'})));
if ~isempty(NUM)
    A=GetV('A');
    [handles.SA.Checked,handles.LS.Checked,handles.LF.Checked,handles.MW.Checked]=deal('off');
    [handles.Redu.Enable,handles.T0.Enable,handles.MSA.Checked]=deal('on');
    A.Nsubiter=NUM(2);
    A.Niter=NUM(1);
    A.Niterf=NUM(3);
    A.NN=NUM(4);
    LoadV('A',A,'Tinv',3);
end
open=GetV('open');
if open==0
    choice = questdlg(sprintf('The inversion method is designed to run in parallel.\n\nTo use this method in parallel?'));
    if strcmp('Yes',choice);
        pop=GetV('NUMW');
        NUMW=str2double(cell2mat(inputdlg({sprintf('Enter number of Workers (Labs): (>=2 & <512)')},...
            ' ',1,{num2str(pop)})));
        if isnan(NUMW) || NUMW<2
        else
            mss=msgbox({sprintf('The Parallelization Settings will take a few seconds.\n\nPlease Wait')},'HV-inv');
            delete(findobj(mss,'string','OK'));
            delete(findobj(mss,'style','frame'));
            set(mss, 'CloseRequestFcn','')
            [handles.Poff.Checked,handles.edit1.Enable,handles.TDC.Enable,handles.THV.Enable,...
                handles.FHV.Enable,handles.THVDC.Enable,handles.PS.Enable,handles.Models.Enable,...
                handles.Inversion.Enable,handles.LP.Enable,handles.SPara.Enable,handles.START.Enable,...
                handles.STOP.Enable,handles.Menu.Enable]=deal('off');
            parallel(NUMW);
            delete(mss);
            LoadV('NUMW',NUMW,'open',1);
            handles.textNumwork.String=['Connection with',' ',num2str(NUMW),' ', 'Labs.'];
            [handles.textNumwork.Visible,handles.edit1.Enable,handles.TDC.Enable,handles.THV.Enable,handles.FHV.Enable,...
                handles.THVDC.Enable,handles.PS.Enable,handles.Models.Enable,handles.Inversion.Enable,...
                handles.LP.Enable,handles.SPara.Enable,handles.START.Enable,handles.STOP.Enable,...
                handles.Menu.Enable,handles.Pon.Checked]=deal('on');
        end
    end
end

function LS_Callback(~, ~, handles)
%% Function for Simplex Downhill menu
NUM=str2double(inputdlg({'Maximum number of iterations allowed','Maximum number of evaluations allowed',...
    'Termination tolerance on the function value (TolFun)','Termination tolerance on parameters (TolX)'},...
    ' ',1,{'500','500','1.e-34','1.e-34'}));
if ~isempty(NUM)
    [handles.MSA.Checked,handles.SA.Checked,handles.MW.Checked,...
        handles.LF.Checked,handles.Redu.Enable,handles.T0.Enable]=deal('off');
    [A.Niter,A.Nsubiter,A.Niterf,A.NN,handles.LS.Checked]=deal(NUM(1),NUM(2),NUM(3),NUM(4),'on');
    LoadV('Tinv',4,'A',A);
end

%% Function for Interior point (IP) menu
function LF_Callback(~, ~, handles)
NUM=str2double(inputdlg({'Maximum number of iterations allowed','Maximum number of evaluations allowed',...
    'Termination tolerance on the function value (TolFun)','Termination tolerance on parameters (TolX)'},...
    ' ',1,{'500','500','1.e-34','1.e-34'}));
if ~isempty(NUM)
    [handles.MSA.Checked,handles.SA.Checked,handles.MW.Checked,...
        handles.LS.Checked,handles.Redu.Enable,handles.T0.Enable]=deal('off');
    [A.Niter,A.Nsubiter,A.Niterf,A.NN,handles.LF.Checked]=deal(NUM(1),NUM(2),NUM(3),NUM(4),'on');
    LoadV('Tinv',5,'A',A);
end

%% Function for the "Wave Parameters" menu
function WPara_Callback(~, ~, ~)
NUM=str2double(inputdlg({'Rayleigh waves modes','Love wave modes','Minimum number of integration points for BW' ,'Maximum Number of integration points for BW','Regularization factor'},...
    'Waves Parameters',1,{'5','5','1000','2000','0.001'}));
if ~isempty(NUM)
    [NR,NL,dk,mdk,apsv]=deal(NUM(1),NUM(2),NUM(3),NUM(4),NUM(5));
    if NR<0 || NL <0
        errordlg('The minimum number of waves modes should be greater to zero','Error');
    end
    if dk>mdk
        errordlg('The maximun number of integration should be greater or equal to minimum number of integration ','Error');
        dk=mdk;
    end
    LoadV('ramasR',NR,'ramasL',NL,'kb_by_dk',dk,'mkb_by_dk',mdk,'apsv',apsv)
end

%% Function for quitting the program
function figure1_CloseRequestFcn(~, ~, ~)
choice = questdlg('Quitting?','Are you Sure?','Yes','No','No');
switch choice
    case 'Yes'
        delete('etc\*.txt');
        closereq
    otherwise
        return
end

%% Three functions for the "L_V Zones" Menu
%% Function for Low S-wave velocity zones
function LZVS_ButtonDownFcn(~, ~, handles)
%ZLVp=GetV('ZLVp');
if strcmp(get(handles.LZVS,'Checked'),'off')
    handles.LZVS.Checked='on';
    ZLVs=1;
    LoadV('ZLVs',ZLVs);
else
    handles.LZVS.Checked='off';
    ZLVs=0;
    %     HHS=1; set(handles.HHS,'Checked','on')
    handles.LZVP.Checked='on';
    ZLVp=1;
    LoadV('ZLVs',ZLVs,'ZLVp',ZLVp);
end


%% Function for Low P-wave velocity zones
function LZVP_ButtonDownFcn(~, ~, handles)
%ZLVs=GetV('ZLVs');
if strcmp(get(handles.LZVP,'Checked'),'off')
    handles.LZVP.Checked='on';
    ZLVp=1;
    LoadV('ZLVp',ZLVp);
else
    ZLVp=0;
    handles.LZVP.Checked='off';
    ZLVs=1;
    handles.LZVS.Checked='on';
    % To be removed in future version:
    handles.HHS.Checked='off';
    HHS=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LoadV('ZLVs',ZLVs,'ZLVp',ZLVp,'HHS',HHS);
end

%% Function for Hard half space
function HHS_ButtonDownFcn(~, ~, handles)
if strcmp(get(handles.HHS,'Checked'),'off')
    handles.HHS.Checked='on';
    handles.LZVS.Checked='on';
    HHS=1;
    ZLVs=1;
    % To be removed in future version:
    ZLVp=1;
    handles.LZVP.Checked='on';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LoadV('HHS',HHS,'ZLVs',ZLVs,'ZLVp',ZLVp);
else
    handles.HHS.Checked='off';
    handles.LZVS.Checked='on';
    HHS=0;
    %ZLVs=1;
    %LoadV('HHS',HHS,'ZLVs',ZLVs);
    LoadV('HHS',HHS);
end

function varargout = HVTI_OutputFcn(~, ~, ~)
%% Function for "Help-About" menu
%varargout{1} = handles.output;
varargout{1} = {};
ZOMBIE=GetV('ZOMBIE');
if ~ZOMBIE.ON
    About_Callback([],[],[]);
end

function About_Callback(~, ~, ~)
msgbox(sprintf(['HVInv Project (v 2.5) 2014-2018\n\n',...
    'This program is free software: you can redistribute it and/or modify\n',...
    'it under the terms of the GNU General Public License version 3 as\n',...
    'published by the Free Software Foundation. For details see the copy\n',...
    'included in this pakage.\n\n',...
    'This program is distributed in the hope that it will be useful,\n',...
    'but WITHOUT ANY WARRANTY; without even the implied warranty of\n',...
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\nSee the ',...
    'GNU General Public License for more details.\n\n',...
    'You should have received a copy of the GNU General Public License\n',...
    'along with this program.  If not, see <http://www.gnu.org/licenses/>..\n\n',...
    'Contributors:\n\n',...
    'Antonio Garc',char(237),'a-Jerez (agarcia-jerez@ual.es)\n',...
    'Jos',char(233), ' Pi',char(241),'a-Flores (ead2009@hotmail.com)\n',...
    'Mathieu Perton (mathieu.perton@gmail.com)\n',...
    'Francisco J. S',char(225),'nchez-Sesma (sesma@unam.mx)\n',...
    'Francisco Luz',char(243),'n (fluzon@ual.es)\n',...
    'Marc Wathelet\n\n',...
    'Forward calculation based on original ideas of:\n\n',...
    'F. S',char(225),'nchez-Sesma, A. Garc',char(237),'a-Jerez and M. Perton\n\n',...
    'Support for GUI and inverse problem:\n\nJos',char(233), ' Pi',char(241),'a-Flores\n\n',...
    'Support for forward calculations:\n\nA. Garc',char(237),'a-Jerez\n\n',...
    'More information at http://www.ual.es/GruposInv/hv-inv/\n\n'...
    'MATLAB' char(174),',', 'Version: 8.5.0.197613 (R2015a) &  9.0.0 341360 (R2016a)\n'...
    'License Numbers: 125381, 40274058, 40462396 \n\n\n'...
    char(169),' Jos',char(233), ' Pi',char(241),'a-Flores & Antonio Garc',char(237),'a-Jerez, 2014-2018\n']),'About');

%% Open User manual
function Manual_Callback(~, ~, ~)
open('./Document/UserManual.html')

function smodel_Callback(~, ~, ~)
%% Save optimum model in *.txt format
[aopt,ZOMBIE]=GetV('aopt','ZOMBIE');
if(ZOMBIE.ON)
    file='zombieBEST.txt';
    path='';
else
    [file,path] = uiputfile('Model.txt','Save file name');
end
if ~isempty(file)
    tex=[path,file];
    fileID = fopen(tex,'w');
    fprintf(fileID,'#Number of layers, halfspace included');
    fprintf(fileID,'\n');
    fprintf(fileID,['#Thickness [m] ' ' #Vp [m/s] ' ' #Vs [m/s] '  ' #Density [kg/m3] ']);
    fprintf(fileID,'\n');
    fprintf(fileID,num2str(size(aopt,1)));
    fprintf(fileID,'\n');
    %fclose (fileID);
    %dlmwrite(tex, aopt ,'-append', 'delimiter', ' ','precision','%8.2f');
    fprintf(fileID,'%8.2f %8.2f %8.2f %8.2f\n',aopt');    
    fclose (fileID);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Additional GUI Functions %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NIP_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



