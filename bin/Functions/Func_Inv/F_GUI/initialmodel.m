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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function  for input for initial model by user
%% Function for initialmodel GUI

function varargout = initialmodel(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @initialmodel_OpeningFcn, ...
    'gui_OutputFcn',  @initialmodel_OutputFcn, ...
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

%% GUI Function elastic parameters by default
function initialmodel_OpeningFcn(hObject, ~, handles, varargin)
% Suggest initial model within the table of Model Parameter ranges if it
% exists or ask for model file otherwise
handles.output = hObject;
guidata(hObject, handles);
movegui(handles.figure1,'center');
[NC,ZOMBIE,ISDEPTH]=GetV('NC','ZOMBIE','ISDEPTH');
%% Elastic parameters
if NC~=0
    % Suggest initial model from the settings in the main table of ranges    
    var=GetV('var');
    iniTry=var(:,1:3:10)+var(:,3:3:12)./2;        
    if ISDEPTH 
        % Make a valid model in terms of depths
        for index=1:NC-1
            if index==NC-1
                overlap=0;
            else
                overlap=var(index,2)-var(index+1,1);
            end
            if overlap>0,
                if index>1
                    iniTry(index,1)=mean([max(iniTry(index-1),var(index,1)) var(index+1,1)+overlap/2]);
                else
                    iniTry(index,1)=mean([var(index,1) var(index+1,1)+overlap/2]);
                end
            else
                if index>1
                    iniTry(index,1)=mean([max(iniTry(index-1),var(index,1)) var(index,2)]);
                else
                    iniTry(index,1)=mean([var(index,1) var(index,2)]);
                end                
            end
        end
        % Change to thickness (standard for models)
        iniTry(:,1)=[iniTry(1,1);diff(iniTry(1:end-1,1));0];        
    end
    iniTry(end,2:3)=var(end,[4 7])+var(end,[6 9])*0.66;% Try to suggest a stiffer halfspace
    iniTry=round(iniTry,2);% Thickness, Vp, Vs % density with two decimals from the begining
    iniTry(:,5)=round((0.5-(iniTry(:,3)./iniTry(:,2)).^2)./(1-(iniTry(:,3)./iniTry(:,2)).^2),4);% Poisson's ratio, 4 decimal
    %Edit table of model parameters
    LoadV('iniTry',iniTry); % Save 
    handles.uitable2.Data=num2cell(iniTry);
else
    [handles.uitable2.Enable,handles.pushbutton1.Enable]=deal('off');
    handles.text2.String='Input the Initial model from file *.txt';
end
if ZOMBIE.ON
    pushbutton3_Callback([],[],handles);
end

function pushbutton3_Callback(~, ~,handles)
% Load initial model from file
ZOMBIE=GetV('ZOMBIE');
if ZOMBIE.ON
    FileName='zombieMDL.txt';
    Path='';
else
    [FileName, Path]=uigetfile({'*.txt'},'Open Parameters');
end
if FileName~=0
    fid=fopen([Path FileName]);
    datos=textscan(fid,'%f','commentstyle','#');
    fclose(fid);
    iniTry=round(cell2mat(datos),2);% Warning can be displayed is rounding is significative
    iniTry=(reshape(iniTry(2:end),4,[]))';
    iniTry(:,5)=round((0.5-(iniTry(:,3)./iniTry(:,2)).^2)./(1-(iniTry(:,3)./iniTry(:,2)).^2),4);
    iniTry(end,1)=0; 
    LoadV('iniTry',iniTry);
    handles.uitable2.Data=num2cell(iniTry);
    [handles.uitable2.Enable,handles.pushbutton1.Enable]=deal('on');
    handles.uitable2.ColumnEditable=true(1,5);
end
if ZOMBIE.ON
    pushbutton1_Callback([],[],[]); %Accept initial model and close
end

%% Accept initial model
function pushbutton1_Callback(~,~,~)
[iniTry,handles2,NC,ZOMBIE,ISDEPTH]=GetV('iniTry','handles','NC','ZOMBIE','ISDEPTH');
NC_IM=size(iniTry,1);
if NC>0 && NC~=NC_IM && ~ZOMBIE.ON, % Number of layers of the inital model does not meet model parametes table
    choice = questdlg('The number of layers doesn`t match the Model Parameters table. The table will be restarted.','Warning','Ok','Cancel','Ok');
    % Handle response
    switch choice
        case 'Cancel'
            close;
            evalin('base', 'clear iniTry');
            return
    end
end
% Permision granted:
LoadV('ini1',iniTry(:,1:4)); % Set initial model
if NC~=NC_IM, % Updating Model Parameters table is required
    [handles2.START.Visible,handles2.TableParameters.Enable,...
        handles2.SPara.Enable,handles2.LV.Enable,handles2.T0.Enable,...
        handles2.Redu.Enable,handles2.Inversiont.Enable,handles2.ThickOrDepth.Enable]=deal('On');
    %% Set variables
    handles2.edit1.String=NC_IM;
    var=zeros(NC_IM,15);
    if ISDEPTH, % Initial model is still interpreted as thickness. A range of depths is created
        % We propose moving +-30% of the thickness of the thinner layer separated by the interface
        var(:,1)=[cumsum(iniTry(1:end-1,1))-min(iniTry(1:end-1,1),[iniTry(2:end-1,1);inf])*0.3;0];
        var(:,2)=[cumsum(iniTry(1:end-1,1))+min(iniTry(1:end-1,1),[iniTry(2:end-1,1);inf])*0.3;0];
    else
        var(:,1)=iniTry(:,1)-iniTry(:,1).*0.5;% min thickness
        var(:,2)=iniTry(:,1)+iniTry(:,1).*0.5;% max thickness
    end
    var(end,1:2)=0;% halfspace
    var(:,4)=iniTry(:,2);var(:,5)=iniTry(:,2);% min and max Vp
    var(:,7)=iniTry(:,3);var(:,8)=iniTry(:,3);% min and max Vs
    var(:,10)=iniTry(:,4);var(:,11)=iniTry(:,4);% min and max density
    var(:,13)=0.1000;% min Poisson
    var(:,14)=0.4999;% max Poisson
    % Rounding and computing amplitudes
    var(:,[1 2 4 5 7 8 10 11])=round(var(:,[1 2 4 5 7 8 10 11]),2);
    var(:,13:14)=round(var(:,13:14),4);
    var(:,3:3:12)=round(var(:,2:3:11)-var(:,1:3:10),2);
    var(:,15)=round(var(:,14)-var(:,13),4);     
    handles2.TableParameters.Data=num2cell(var(:,[1 2 4 5 7 8 10 11 13 14]));
    limdd=sum(var(:,2)); %max. depth down to the halfspace
    LoadV('var',var,'limdd',limdd,'NC',NC_IM);
end
evalin('base', 'clear iniTry');
close;

%% GUI Function - Cancel initial model
function pushbutton2_Callback(~, ~,~)
evalin('base', 'clear iniTry');
close;

%% GUI Function edit elastic parameters by user (Initial Model Table parameters)
%% Rounds values modified in the table and address Vp Vs Poisson's ratio dependences
function uitable2_CellEditCallback(hObject,  ~, handles)
dataF = cell2mat(hObject.Data);
iniTry=GetV('iniTry');
change=(dataF(:,1)~=iniTry(:,1));
dataF(change,1)=round(dataF(change,1),2);% Round new Thickness in dataF
change=(dataF(:,2)~=iniTry(:,2));
dataF(change,2)=round(dataF(change,2),2);% Round new Vp in dataF
changeS=(dataF(:,3)~=iniTry(:,3));
dataF(changeS,3)=round(dataF(changeS,3),2);% Round new Vs in dataF
% Update Poisson for layes where either Vp or Vs changed:
change=change|changeS;
dataF(change,5)=round((0.5-(dataF(change,3)./dataF(change,2)).^2)./(1-(dataF(change,3)./dataF(change,2)).^2),4);
% Where neither Vp nor Vs changed, and Poisson changed -> Round Poisson and update Vp
change= ~change & (dataF(:,5)~=iniTry(:,5));
dataF(change,5)=round(dataF(change,5),4);
dataF(change,2)=round(real(dataF(change,3).*(sqrt((2*(1-dataF(change,5)))./(1-(2.*dataF(change,5)))))),2);
% Save dataF to iniTry and to the table:
dataF(end,1)=0;
handles.uitable2.Data=num2cell(dataF);
LoadV('iniTry',dataF);

%% Additional GUI Functions
function varargout = initialmodel_OutputFcn(~, ~, ~)
%varargout{1} = handles.output;
varargout={};
