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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%}

%% Function  for selection inversion method by user
% If no method is selected (closing) we awwignTinv is preserved
%

function varargout = selectmethod(varargin)%#codegen
%% Function for selecthmetod GUI
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @selectmethod_OpeningFcn, ...
    'gui_OutputFcn',  @selectmethod_OutputFcn, ...
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
%% Function GUI
function selectmethod_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
movegui(handles.figure1,'center')
%Tinv=0;% Type Inversion
Tinv=GetV('Tinv');
LoadV('Tinv',-Tinv);
function varargout = selectmethod_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
%% Function for election Metropoli method
function pushbutton1_Callback(~, ~, ~)
LoadV('Tinv',1);
close
%% Function for election SA method
function pushbutton2_Callback(~, ~, ~)
LoadV('Tinv',2);
close
%% Function for election SAM method
function pushbutton3_Callback(~, ~, ~)
LoadV('Tinv',3);
close
%% Function for election Simplex Downhill method
function pushbutton4_Callback(~, ~, ~)
LoadV('Tinv',4);
close
%% Function for election Interior-point method
function pushbutton5_Callback(~, ~, ~)
LoadV('Tinv',5);
close