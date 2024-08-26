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


%% Graphical interface for selecting
%% the type of dispersion curve
%% for individual or joint inversion

%% CDE function (GUI)
function varargout = CDE(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CDE_OpeningFcn, ...
    'gui_OutputFcn',  @CDE_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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

%% Variables by default
function CDE_OpeningFcn(hObject, ~, handles, ~)
handles.output = hObject;
guidata(hObject, handles);
%% Number of curve mode
%MODE=0; % MODE==0 Fundamental mode
%% Weight percentage for joint inversion
% ww=0.5; % Porcentage Individual (Joint inversion)
% pmy=1;% (pmy=1 asigna el peso en función de las muestras, pmy=0 asigna peso segun el usuario)
% NFile=0; % Número de archivos (Inabilitado)
%% Assigning environment variables GUI
[handles.PESO.Enable,handles.MODE.Enable,handles.pushbutton1.Enable,...
    handles.RAY.Enable,handles.LOV.Enable,handles.PH.Enable, handles.GR.Enable,...
    handles.checkbox1.Enable]=deal('off');% Control de los elementos de la gui
[handles.MODE.String,handles.PESO.String]=deal(0,0.5);
LoadV('MODE',0,'pmy',1,'ww',0.5,'NFile',0,'cancel',0);% Asignando variables al workspace
[FileName, Path]=uigetfile({'*.txt'},'Open Dispersion Curve');% Datos de la curva de dispersion desde archivo
if FileName~=0 % Si no esta vacio en nombre del archivo, continua con la operación
    % Control de los elementos de la GUI
    [handles.pushbutton1.Enable,handles.RAY.Enable,handles.LOV.Enable,...
        handles.PH.Enable,handles.GR.Enable]=deal('on');
    handles.text3.String=FileName;% Muestra el nombre de archivo
    
    flag=GetV('flag');% Flag de inversion conjunta flag=3 ,flag=1 
    if flag==2 || flag==1;
        [handles.PESO.Enable,handles.checkbox1.Enable]=deal( 'off');%Deshabilita el valor del peso de las curvas en la inversion
    else
        handles.checkbox1.Enable='on';%habilita el valor del peso de las curvas en la inversion
    end
    LoadV('FileName',FileName,'Path',Path);% Carga ruta y nombre de archivo para su lectura posterior
    movegui(handles.figure1,'center')
else
    LoadV('cancel',1);% Cancelacion de operacion por el usuario
end
%% Module polarization (Rayleigh)
function POL_CreateFcn(~, ~, ~)
%% value to default
LoadV('POL',1);
function POL_SelectionChangeFcn(hObject, ~, handles)
%% Assigning the value of polarization of Rayleigh curve by user
if hObject==handles.RAY
    POL=1; % Rayleigh
else
    POL=2; % Love
end
LoadV('POL',POL)
%% %% Assigning the value of velocity of phase or group by user
function VEL_CreateFcn(~, ~, ~)
%% value to default
LoadV('VEL',1);
function VEL_SelectionChangeFcn(hObject, ~, handles)
%% Assigning the value of velocity
if hObject==handles.PH
    VEL=1; % Phase Velocity
else
    VEL=2;% Group velocity
end
LoadV('VEL',VEL);
%% Module number assignment mode (Deshabilitado)
function MODE_Callback(hObject, ~, ~)
MODE = str2double(hObject.String);
if isnan(MODE) % Si en nan 
    [hObject.String,MODE]=deal(0);
    errordlg('Input must be a number','Error');
elseif MODE <0 %si MODE es menor a cero
    [hObject.String,MODE]=deal(0);
    errordlg('The number of Mode Index should be greater to 0 ','Error')
end
LoadV('MODE',MODE);% MODE=0 Fundamental mode
%% Weight allocation module in joint inversion
function PESO_Callback(hObject, ~, ~)
ww = str2double(hObject.String);% Obtiene valor del factor de peso introducido por el usuario
if isnan(ww)
    [ww,hObject.String]=deal(0.5);
    errordlg('Input must be a number','Error');
elseif ww<0 || ww>1
    [ww,hObject.String]=deal(0.5);
    errordlg('The  weight should be between 0 and 1 ','Error')
end
LoadV('ww',ww);%Factor peso de las curvas en la inversion conjunta
%% Enabling the weight factor
function checkbox1_Callback(~, ~, handles)
if handles.checkbox1.Value==1;% Factor de peso en funcion de las muestras
    pmy=1;
    handles.PESO.Enable= 'off';% GUI
else % Factor de peso en funcion del usuario
    pmy=0;
    handles.PESO.Enable= 'on';%GUI
end
LoadV('pmy',pmy);
%% Load next curve (Disable)
function Next_Callback(~, ~, handles)
 if handles.Next.Value==1;
%     NFile=1;
% else
%     NFile=0;
 end
% LoadV('NFile',NFile);
%% Button of accept process
function pushbutton1_Callback(~, ~, ~)
hfig=get(0,'CurrentFigure');
delete(hfig);% Eliminacion de la gui secundaria
%%  Button of cancel process
function pushbutton2_Callback(~, ~, ~)
LoadV('cancel',1); % cancelacion de opreacion por el usuario
hfig=get(0,'CurrentFigure');
delete(hfig);% Eliminacion de la gui secundaria
%% Open data file (deshabilitado)
function pushbutton3_KeyPressFcn (~, ~,~)
% [FileName, Path]=uigetfile({'*.txt'},'Open Dispersion Curve');
% if FileName~=0
%     [handles.pushbutton1.Enable,handles.RAY.Enable,handles.LOV.Enable,...
%         handles.PH.Enable,handles.GR.Enable ]=deal('on');
%     handles.text3.String=FileName;
%     yes=GetV('yes');
%     if yes==2;
%         [handles.PESO.Enable,handles.checkbox1.Enable]=deal('off');
%     else
%         handles.checkbox1.Enable='On';
%     end
%     LoadV('FileName',FileName,'Path',Path);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary functions of gui
function pushbutton3_Callback(~, ~, ~)
function figure1_CloseRequestFcn(~, ~, ~)
function figure1_ResizeFcn(~, ~, ~)
function LOV_CreateFcn(~, ~, ~)
function MODE_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function varargout = CDE_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
function PESO_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end