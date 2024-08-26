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
%% Function for input standard deviation of curve by user
function varargout = STDG(varargin)
%% Function for STD GUI
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @STDG_OpeningFcn, ...
    'gui_OutputFcn',  @STDG_OutputFcn, ...
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
function STDG_OpeningFcn(hObject,~, handles, varargin)
%% Default settings
%Vstd=10; %Value Standar desviation default
%stdg=3;% Flag for SD: std=1 constant values; std=2 proportional values; std=3 'no' desviation (1);
handles.output = hObject;
guidata(hObject, handles);
movegui(handles.figure1,'center')
LoadV('Vstd',10,'stdg',3);
handles.Vstd.String=10;
%% GUI Function edit SD by user
function Vstd_Callback(hObject, ~, ~)
Vstd = str2double(hObject.String);
if isnan(Vstd)
    hObject.String= 10;
    Vstd =10;
    errordlg('Input must be a number','Error');
elseif Vstd <0
    hObject.String= 10;
    Vstd =10;
    errordlg('The number of Standard deviation should be greater to 0 ','Error')
end
LoadV('Vstd',Vstd);
%% GUI Function assigning SD by user
function Cstd_Callback(~, ~, ~)
LoadV('stdg',1);
close
function Pstd_Callback(~, ~, ~)
LoadV('stdg',2);
close
%% GUI Function cancel SD by user
function Cancel_Callback(~,~,~)
LoadV('Vstd',1,'std',3);
close
%% Additional GUI Functions
function varargout = STDG_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
function Vstd_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end