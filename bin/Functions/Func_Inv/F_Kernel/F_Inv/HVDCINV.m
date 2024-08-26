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
%      #ok<*DEFNU>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for inversion methods and GUI control

function HVDCINV(handles)
%INPUT
%handles elements of GUI
%ER structure for information during inversion
%Cont contador de iteraciones
%%  Getting variables from the workspace
% A = Information of iterations
% A.open= Enable-disable parallel
% G.Pl= Type axes control
% C= Information of targets curves
% A.apsv= Controls a tiny imaginary part of frequency to stabilize body-wave integrals
% pmy= equalize misfits by the number of samples (1=yes)
% C.ww= weight of dispersion curve in joint inversion
% ini1= initial model (for some methods only)
% A.HHS= Hard Half Space option (1=yes or 0=no)
% A.ZLVs= Zone of Low S-wave Velocity
% A.ZLVp= Zone of Low P-wave Velocity
% ISDEPTH= deals with layer thicknesses (false) or depths (true)
% A.NC=%Number of layers
% A.var=%Range of parameters
% A.T0=%Method for temperature (1-> from A.EA, A.AP; or 2->from initial temp in A.AP)
% A.EA=;%Error aceptation (SA MSA)
% A.AP=%aceptaction Probability  (SA MSA)
% A.RT=%Reduce temperature
% Tinv= Type of inversion (1=Monte Carlo 2=SA 3=MSA 4=Downhill Simplex 5=Interior Point)
% ZOMBIE= Flag for automatic operation of the program 
%Get Variables from workspace
[A,A.open,G.Pl,C,A.apsv,pmy,C.ww,ini1,A.HHS,A.ZLVs,A.ZLVp,A.ISDEPTH,A.NC,A.var,A.T0,A.EA,A.AP,...
    A.RT,Tinv,ZOMBIE]=GetV('A','open','Pl','C','apsv','pmy','ww','ini1','HHS','ZLVs',...
    'ZLVp','ISDEPTH','NC','var', 'T0','EA','AP', 'RT','Tinv','ZOMBIE');

%% Calculation of weights for joint inversion when equalization is required
if pmy &&  ~isempty(C.HVFobs) && ~isempty(C.DCFobs)
    % The weight for the dispersion curve will be:
    C.ww=C.nmHV/(C.nmDC+C.nmHV);% = const/C.nmDC, with const=C.nmDC*C.nmHV/(C.nmDC+C.nmHV)
    % The weight for H/V will be (1-C.ww) = const/C.nmHV
end

% Variables for history control
%%%%%%%%%%%%%
ER.MODALL=[]; % Se usa para las gráficas de sensitiv. Todos los modelos evaluados.
              % To be saved in the base ws as MODALL
ER.L2ALL=[];
ER.DCALL=[];
ER.HVALL=[];
%%%%%%%%%%%%%
ER.err=[];
ER.err1=[];
ER.temer=[];
ER.kk=[];
Cont.iter=0;% Total iterations: sum methods, only one for each step in MSA    L2 es No NaN y No infinito
Cont.all=0;% Total iterations: sum methods, NW for each step in MSA     
Cont.valid=0;% Sum methods? L2 es no cero y no infinito
Cont.IterLastMethod=0;% Counting iterations for the curret method

%% flag for variables HV curve
if  ~isempty(C.HVFobs)
    %Informacion de la diferente contribucion de ondas para el calculo del
    %HV
    [C.NR,C.NL,C.dk,C.mdk]=GetV('ramasR','ramasL','kb_by_dk','mkb_by_dk');
end
%% flag for variables DC curve
if  ~isempty(C.DCFobs)
    % Obteniendo informacion del tipo de curva de dispersión
    [C.POL,C.VEL,C.MODE]=GetV('POL','VEL','MODE');
end
A.ST=0;% Cumulative number of sub-iterations (Temperatures). Cumulative for several inversions
if A.ISDEPTH
    A.limd=1.1*A.var(end-1,2);% Limit for plotting halfspace    
else
    A.limd= sum(A.var(:,2))+0.1*sum(A.var(:,2));% Limit for plotting halfspace
end
LoadV('limdd',A.limd);
if isempty(ini1)% If no initial model is given
    INPO=GetV('INPO');% Initial population for random search
    if INPO<1,INPO=1;LoadV('INPO',INPO);end
else
    INPO=1;
end
%% Control processes in parallel.
if A.open==0
    NUMW=1;
    INPOround=INPO;
    NBarSteps=20;% 20 means bar steps of 5%
    NModelsPerBarStep=round(INPO/NBarSteps);
else
    NUMW=GetV('NUMW');% Numero de nucleos habilitados para procesos en paralelo.
    NBarSteps=20;
    NModelsPerBarStep=round(NUMW*ceil(INPO/(NBarSteps*NUMW)));
    INPOround=NModelsPerBarStep*NBarSteps;
end
HVT=zeros(C.nmHV,INPOround);%HV curve
DCT=zeros(C.nmDC,INPOround);%Disp curve
HV=zeros(C.nmHV,1);%HV curve
DC=zeros(C.nmDC,1);%Disp curve
DCTp=zeros(C.nmDC,NModelsPerBarStep);%DC curve in parallel
HVTp=zeros(C.nmHV,NModelsPerBarStep);%HV curve in parallel
mod=nan(A.NC,4, INPOround);
modp=zeros(A.NC,4,NModelsPerBarStep)*nan; %model in parallel
L2i=nan(1, INPOround)*nan;% norm L2
L2ip=nan(1,NModelsPerBarStep);% norm L2 in parallel
A.flag=zeros(NUMW,1);% Flag for SAM in parallel
A.in=zeros(A.NC,4);% initial Model for perturbation

% Initializing auxiliar matrices for the 'Low Velocity Zones' and 'Hard Halfspace' conditions
% A.HHS==true means that the Vs of the halfspace must be the largest one. 
% This flag can be changed below to false for not implemented cases
if find(A.var(1:A.NC-1,7)>A.var(A.NC,8)),errordlg('Impossible to set hard-halfspace condition. Turning that condition to false');uiwait;A.HHS=false;end
if A.ZLVs&&A.ZLVp&&A.HHS % only case with HHS condition implemented
    [A.pol,A.limseg]=deal([]);
    [~,~,~,A.polHHS,A.facHHS,A.ProbLimSupSegHHS,A.LimVelSegHHS]=ini1_plus(A.NC,A.var,A.ZLVs,A.ZLVp,[],[],A.HHS);
else
    A.HHS=false;[A.polHHS,A.facHHS,A.ProbLimSupSegHHS,A.LimVelSegHHS]=deal([]);
    [~,A.pol,A.limseg]=ini1_plus(A.NC,A.var,A.ZLVs,A.ZLVp);
end

%% Generation of random models (random search) or evaluation of the initial model misfit
pass=1; % No error yet
% Waiting Message
if isempty(ini1);% (if an initial model is not given)
    BarId=waitbar(0,{'Initial random exploration of the parameters space','Please wait'},'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(BarId,'canceling',0);
    set(BarId, 'CloseRequestFcn','')
end
handles.STOP.Enable='off';
if A.open==0 || ~isempty(ini1), % Serial computation
%fprintf(1,'Serie\n');
    if isempty(ini1)
        LastBarStep=uint32(0);
        for i=1:INPO
            if getappdata(BarId,'canceling') % Stopped by user
                [handles.STOP.UserData, handles.START.Visible]=deal(1,'on');
                delete(BarId);
                pass=0;
                break;
            end
            [mod(:,:,i),HV,DC,~,L2i(:,i)]=CTF(A,C,i,1,G);
            [DCT(:,i),HVT(:,i)]=deal(DC,HV);
            if floor(i/NModelsPerBarStep)>LastBarStep
                waitbar(i/INPO,BarId);
                LastBarStep=LastBarStep+1;
            end
        end
    else    
        i=1;
        [mod(:,:,i),HV(:,i),DC(:,i),~,L2i(:,i)]=CTF(A,C,i,3,G);
        [DCT(:,i),HVT(:,i)]=deal(DC,HV);
        if isnan(L2i);
            errordlg(sprintf('Forward calculation failed for the initial model\n\nPlease, enter a new initial model'),'Error');
            [handles.STOP.UserData, handles.START.Visible]=deal(1,'on');
            pass=0;
        end
    end
%toc
%INPOround
else % Parallel computation
%fprintf(1,'paralelo\n');
%tic;
    kij=0;
    for ii=1:NBarSteps
        % Cancelel by user
        if getappdata(BarId,'canceling')
            [handles.STOP.UserData, handles.START.Visible]=deal(1,'on');
            delete(BarId);
            pass=0;
            break;
        end
        % processes in parallel
        parfor i=1:NModelsPerBarStep
            [modp(:,:,i),HV,DC,~,L2ip(:,i)]=CTF(A,C,i,1,G);
            DCTp(:,i)=DC;
            HVTp(:,i)=HV;
        end
        mod(:,:,((kij+1):(ii*NModelsPerBarStep)))=modp;
        L2i(((kij+1):(ii*NModelsPerBarStep)))=L2ip;
        % Cálculo para las curvas de dispersión modelos y misfit
        if  ~isempty(C.DCFobs)
            DCT(:,((kij+1):(ii*NModelsPerBarStep)))=DCTp;
        end
        %Calculo para las curvas HV modelos y misfit
        if  ~isempty(C.HVFobs)
            HVT(:,((kij+1):(ii*NModelsPerBarStep)))=HVTp;
        end
        kij=ii*NModelsPerBarStep;
        waitbar(ii/NBarSteps,BarId);
    end
    %toc
    %INPOround
end
if isempty(ini1);
    delete(BarId);
end
handles.STOP.Enable= 'on';% STOP button enabled
% Eliminando archivos de lectura
delete('etc\*.txt');
% Sorting
[~ ,b]=sort(L2i);% Sorting the misfit
aopt=mod(:,:,b(1));% Taking the best model
[DCopt,HVopt]=deal(DCT(:,b(1)),HVT(:,b(1)));% Best DC and HV
A.L2opt=L2i(b(1));% Taking the lower misfit. This copy will not be updated.
%Clean Memory
clear HVTp DCTp modp L2ip HVT DCT mod L2i

%% Preparing for subsequent algorithms
% L2i1= Auxiliar norm L2
% A.L2optimo= Best model misfit. To be updated by algorithms
[L2i1,A.L2optimo]=deal(A.L2opt);
handles.textmisfit.String=['#Evaluated models',' ','1',' ','Min. Misfit',' ',num2str(A.L2optimo)];
A.aini=aopt;% initial model
% A.NN1=A.NN/100;% Perturbation of range in percentage
% A.NN2=A.NN/100;% Perturbation of range in percentage
[A.NN1,A.NN2]=deal(A.NN/100);
if ~isempty(ini1)
    %[modini1,HVTini1,DCTini1,L2ini1]=deal(mod(:,:,1),HVT(:,1),DCT(:,1),L2i(:,1));
    [modini1,HVTini1,DCTini1,L2ini1]=deal(aopt,HVopt,DCopt,A.L2opt);
end

%% GUI graphing module
% GUI graphical error
[F1,F2,Perr]=GetV('F1','F2','Perr');
if pass==1
    if isempty(F1);
        handles.ploterror.NextPlot='add';
        [Perr,F1,F2]=plotyy(handles.ploterror,0,[0; 0],0,0);
        uistack(Perr(1),'top');
        %% Graphic settings
        % Escalas automaticas
        [Perr(2).YScale,Perr(1).YScale]=deal('Log');
        [F1(2).Parent.YTickMode,Perr(1).YLimMode,Perr(1).XLimMode,Perr(1).XTickMode,Perr(2).YLimMode,...
            Perr(2).XLimMode,Perr(2).XTickMode,handles.ploterror.XLimMode,handles.ploterror.YLimMode]=deal('auto');
        % Legendas de ejes
        Perr(1).YLabel.String='Misfit';
        if Tinv==2 || Tinv==3
            Perr(2).YLabel.String='Temperature';
        elseif Tinv==5 || Tinv==4 || Tinv==1
            [Perr(2).YLabel.String,Perr(2).YTickMode,Perr(2).YTick]=deal(' ','manual',[]);
        end
        %Legend for misfit evolution error graphic
        fgd1=legend(handles.ploterror,'Minimum misfit','Current misfit');
        %% Graphic settings
        % Type of lines, ...
        [Perr(1).Color,Perr(2).Color,Perr(2).YColor,F1(2).Parent.YColor,F1(1).LineWidth,...
            F1(1).Color,F1(1).DisplayName,F1(2).LineWidth,F1(2).Color,F1(2).DisplayName,...
            F2.LineWidth,F2.Color,fgd1.FontSize,fgd1.Box,handles.ploterror.XLabel.String]=...
            deal('none','w','b','r',2.5,'r','Minimum misfit',1.5,'m','Current misfit',2.5,'b',10,'off','Iteration');
        % Type interpreter Latex and normalized legend
        [fgd1.Interpreter,Perr(2).YLabel.Interpreter,handles.ploterror.XLabel.Interpreter,...
            Perr(1).YLabel.Interpreter]=deal('latex');
        [Perr(1).YLabel.Units,Perr(2).YLabel.Units,handles.ploterror.XLabel.Units,...
            handles.ploterror.FontUnits]=deal('normalized');
        uistack(Perr(2),'down');
        linkaxes([Perr(1),Perr(2)],'x');
        % Load Id graphics
        LoadV('F1',F1,'F2',F2,'Perr',Perr);
    else
        [F1(1).XData,F1(2).XData,F2.XData]=deal(ER.kk);% Number of iterations
        F1(1).YData=ER.err;% Minimum misfit
        F1(2).YData=ER.err1;% Current misfit
        F2.YData=ER.temer;% Temperature
        linkaxes([Perr(1),Perr(2)],'x');
        [handles.ploterror.YLimMode,handles.ploterror.XLimMode]=deal('auto');
        if Tinv==2 || Tinv==3
            Perr(2).YLabel.String='Temperature';
        elseif Tinv==4 || Tinv==1 || Tinv==5
            [Perr(2).YLabel.String,Perr(2).YTickMode,Perr(2).YTick]=deal(' ','manual',[]);
        end
    end
    
    %% GUI graphical HV,Disp curve & models
    delete(findobj(handles.PROFILE,'type','axes'))
    ppv=subplot(1,1,1,'Parent',handles.PROFILE);
    [BTA,ALFA,~,E]=MODELO(aopt,GetV('limdd'));
    G.ptlmod=plot(ppv,BTA,E,'-r',ALFA,E,':r','LineWidth',3);
    
    %% Graphic settings
    [ppv.YLim,ppv.XLabel.String,ppv.YLabel.String,ppv.YDir,ppv.FontUnits]=...
        deal([0 A.limd ],'$Velocity$ [${m} \over {s}$]','$Depth$ [$m$]','reverse','normalized');
    [ppv.YLabel.Interpreter,ppv.XLabel.Interpreter]=deal('latex');
    l3=legend(ppv,'V_s','V_p');
    l3.FontSize=8;
    %Drawing the H/V curve (if it exists)
    if  ~isempty(C.HVFobs)
        if length(G.Pl.HV.Children)>1
            delete(G.Pl.HV.Children(1));
        end
        G.phv=plot(G.Pl.HV,C.HVFobs,HVopt,'-r','LineWidth',4);
        l1=legend (G.Pl.HV,C.HVFile,'Best fitting HV');
        [l1.Interpreter,l1.FontSize,G.phv.Parent.YLimMode]=deal('none',8,'auto');
    end
    %Drawing dispersion curve (if it exists)
    if  ~isempty(C.DCFobs)
        if length(G.Pl.DC.Children)>1
            delete(G.Pl.DC.Children(1));
        end
        G.pdc=plot(G.Pl.DC,C.DCFobs,DCopt,'-r','LineWidth',4);
        l1=legend (G.Pl.DC,C.DCFile,'Best fitting DC');
        [l1.Interpreter,l1.FontSize,G.pdc.Parent.YLimMode]=deal('none',8,'auto');
    end
    drawnow
end

%%  Initialization of variables
[mod,out]=deal(zeros(A.NC,4,NUMW));% Modelos, Modelo perturbado de salida
[L2i2,deltaE]=deal(zeros(1,NUMW));% Norma, Delta de energia SA
DATERROR.lim=A.limd;% Limite de perfil Halfspace
iForStats=0; %Counter for iterations of the last Markov chain, used for statistics, or Monte Carlo
             %It will be reset if a new method is chosen.
Tfin=0;% Tiempo final
Tfin2=0;
resp=1;% Respuesta Usuario a si se continúa la inversión
HVT=zeros(C.nmHV,NUMW);% HV salida function CFT
DCT=zeros(C.nmDC,NUMW);% DC salida function  CFT
G.F1=F1;% Grafica Misfit
G.F2=F2;% Grafica Curren misfit

% Assigning initial value of temperature depending on user's settings (SA)
if A.T0==1
    [T1,Tini,TSAM]=deal(-A.EA*A.L2optimo/log(A.AP)); % Computing the temperature dependence from mismit
else
    [T1,Tini,TSAM]=deal(A.AP);
end
handles.textmisfit.Visible= 'on';

if Tinv~=3
    NUMW=1;
end

%% MAINBLOCK: APPLICATION (POSSIBLY REITERATED) OF INVERSION METHODS. 
while resp==1 % Application of a new inversion method
    if Tinv==1 || Tinv ==2 ||Tinv==3

        %% Kernel of SA, MSA and Monte Carlo methods
        if handles.STOP.UserData==1%User stop
            break;
        end        
        % Increasing storage in ER variables, used for the Misfit History Graphic        
        switch Tinv
            case 1, % Monte Carlo
                ItNew=A.Niterf;
            case 2,% SA    
                ItNew=(A.Nsubiter)*A.Niter+floor(A.Niterf);
                A.ST=A.ST+A.Nsubiter-1;
            otherwise % MSA
                ItNew=(A.Nsubiter)*A.Niter+floor(A.Niterf);
                A.ST=A.ST+A.Nsubiter;
        end
        ER.kk(end+1:end+ItNew,1)=NaN;%[ER.kk; nan(ItNew,1)];
        ER.err(end+1:end+ItNew,1)=NaN;%=[ER.err; nan(ItNew,1)];
        ER.err1(end+1:end+ItNew,1)=NaN;%=[ER.err1; nan(ItNew,1)];
        ER.temer(end+1:end+ItNew,1)=NaN;%=[ER.temer; nan(ItNew,1)];
        ER.DCALL(1:C.nmDC,end+1:end+ItNew)=0;%[ER.DCALL zeros(C.nmDC,ItNew)];
        ER.HVALL(1:C.nmHV,end+1:end+ItNew)=0;%[ER.HVALL zeros(C.nmHV,ItNew)];
        ER.MODALL(1:A.NC,1:4,end+1:end+ItNew)=0;%=cat(3,ER.MODALL, zeros(A.NC,4,ItNew));
        ER.L2ALL(end+1:end+ItNew,1)=NaN;%=[ER.L2ALL ;zeros(ItNew,1)];
        % Initializing storage for DATERROR (new by Antonio) for the Last
        % Markov Chain or iters with Monte Carlo
        DATERROR.mod=nan(size(A.aini,1),size(A.aini,2),A.Niterf);
        DATERROR.Mis=nan(1,A.Niterf);
        % Do the iterations. 
        for j=1:A.Nsubiter+1
            if handles.STOP.UserData==1 %User stop
                break;
            end
            %% If we are in the final interation (Last Markov Chain), updating A.Niter and other parameters
            if j==A.Nsubiter+1
                T1=TSAM;
                A.aini=aopt;
                L2i1=A.L2optimo;
            end
            %% Loop for iterations
            for i = 1:max((j<=A.Nsubiter)*A.Niter,(j>A.Nsubiter)*A.Niterf);
                if handles.STOP.UserData==1%User stop
                    break;
                end
                %% Semi-random generation of HV or Disp curves ( direct problem )
                tic;% Inicio del tiempo de cálculo para estadistica
                if NUMW==1, % Non-parallel ¿NUMW=1?
                        [mod(:,:,1),HV,DC,out(:,:,1),L2i2(:,1)]=CTF(A,C,1,2,G);% forward calculations
                        DCT(:,1)=DC;
                        HVT(:,1)=HV;
                        deltaE=L2i2(:,1)-L2i1;
                else % parallel for MSA if available
                    parfor ij=1:NUMW;
                        [mod(:,:,ij),HV,DC,out(:,:,ij),L2i2(:,ij)]=CTF(A,C,ij,2,G);
                        DCT(:,ij)=DC;
                        HVT(:,ij)=HV;
                        deltaE(ij)=L2i2(:,ij)-L2i1;
                    end
                end
                Tfin=Tfin+toc; %tiempo de cálculo para estadistica
                %% Assigning variables and control parameters
                A.flag(1,1)=0;
                A.NN1=A.NN2;
                A.NN1=A.NN/100;
                Cont.iter=Cont.iter+1;
                Cont.IterLastMethod=Cont.IterLastMethod+1;                
                [~,b2]=sort(deltaE);
                for ij=1:NUMW
                    % Assigning models
                    Cont.all=Cont.all+1;
                    ER.MODALL(:,:,Cont.all)=mod(:,:,NUMW);
                    ER.L2ALL(Cont.all)=L2i2(:,NUMW);
                    ER.DCALL(:,Cont.all)=DCT(:,b2(1));
                    ER.HVALL(:,Cont.all)=HVT(:,b2(1));
                end
                %% Acceptance condition in SA, MSA and Metropoli
                modapro=ceil(rand()*length(b2));
                if deltaE(b2(1)) < 0,% Always accepted
                    A.aini=mod(:,:,b2(1));
                    L2i1=L2i2(b2(1));
                elseif Tinv~=1 && rand < exp(-1*deltaE(modapro)/T1) %% Probabilistic acceptance in SA & MSA
                    A.aini=mod(:,:,b2(modapro));
                    L2i1=L2i2(:,b2(modapro));
                    %% Acceptance condition in Metropoli
                elseif Tinv==1 && rand < exp(-0.5*deltaE(modapro))
                    A.aini=mod(:,:,b2(modapro));
                    L2i1=L2i2(:,b2(modapro));
%                 elseif Tinv==1 %% Acceptance condition in Metropoli, step by step
%                     deltaE
%                     modapro
%                     kkkk=exp(-0.5*deltaE(modapro));
%                     prueba=(rand<kkkk);
%                     fprintf(1,'aceptado si menor de %f, %d\n',kkkk,prueba);
%                     if prueba
%                         A.aini=mod(:,:,b2(modapro));
%                         L2i1=L2i2(:,b2(modapro));
%                     end
                end
                %% Optimal model and reallocation of control parameters
                if L2i2(:,b2(1)) < A.L2optimo,% Hemos mejorado el mejor modelo
                    if Tinv==3
                        A.flag(1,1)=1;
                        A.NN1=A.NN2/2;
                    end
                    % Asignando nuevos valores para la siguiente iteracion
                    A.in=out(:,:,b2(1));
                    aopt=mod(:,:,b2(1));
                    A.aini=mod(:,:,b2(1));
                    A.L2optimo=L2i2(:,b2(1));
                    L2i1=A.L2optimo;
                    DCopt=DCT(:,b2(1));
                    HVopt=HVT(:,b2(1));
                    if Tinv==3 && j~=A.Nsubiter+1, % Temperature reduction for MSA
                        T1=T1*A.RT;
                    end
                    %% Plotting THE BEST model and curves using Monte Carlo, SA & MSA
                    [BTA,ALFA,~,E]=MODELO(aopt,GetV('limdd'));
                    [G.ptlmod(1).YData,G.ptlmod(2).YData]=deal(E);
                    G.ptlmod(1).XData=BTA;
                    G.ptlmod(2).XData=ALFA;
                    if  ~isempty(C.HVFobs)
                        G.phv.YData=HVopt;
                    end
                    if  ~isempty(C.DCFobs)
                        G.pdc.YData=DCopt;
                    end
                    %% Plotting Misfit History 
                    [G.F1(1).XData,G.F1(2).XData,G.F2.XData]=deal(ER.kk);
                    G.F1(1).YData=ER.err;
                    G.F1(2).YData=ER.err1;
                    G.F2.YData=ER.temer;
                    G.F1(1).Parent.XLimMode='Auto';
                    drawnow
                end
                %% Ploting error, misfit & temperature
                if j==A.Nsubiter+1 || Tinv==1% We are in Last Markov chain, used for statistics
                    iForStats=iForStats+1;% contador para la grafica del misfit
                    DATERROR.mod(:,:,iForStats)=A.aini;
                    DATERROR.Mis(iForStats)=L2i1;
                    DATERROR.temp=TSAM;
                    if floor(i/(A.Niterf/10))~=floor((i-1)/(A.Niterf/10))                   
                        %% Plotting Misfit History 
                        [G.F1(1).XData,G.F1(2).XData,G.F2.XData]=deal(ER.kk);
                        G.F1(1).YData=ER.err;
                        G.F1(2).YData=ER.err1;
                        G.F2.YData=ER.temer;
                        G.F1(1).Parent.XLimMode='Auto';
                        drawnow
                    end
                end
                %% Graphing evolution of the error
                ER.kk(Cont.iter)=Cont.iter;
                ER.err1(Cont.iter)=L2i1;
                ER.err(Cont.iter)=A.L2optimo;
                if Tinv~=1
                    ER.temer(Cont.iter)=T1;
                else %NaN for Monte Carlo
                    ER.temer(Cont.iter)=nan;
                end
                %Adjusting Control Parameters of HV curve
                if  ~isempty(C.HVFobs)
                    C.dk=C.dk+1;% Aumentando gradualmente el parametro de dk para el HV
                    if C.dk >C.mdk % dk min supera al dk maximo, se quedan igual
                        C.dk=C.mdk;
                    end
                end
                % Pintando el numero de iteraciones.
                set(handles.textmisfit,'String',['#Evaluated models','  ',num2str(Cont.all),'  ','Min. Misfit','  ',num2str(A.L2optimo)]);
                drawnow
                Tfin2=Tfin2+toc; %tiempo final de cálculo para estadistica
            end
            %% settings of temperature
            if T1<TSAM
                TSAM=T1;
            end
            if Tinv==3 && j~=A.Nsubiter+1;
                if A.Nsubiter~=0;
                    T1=T1*1.5;
                end
            elseif Tinv==2 && j~=A.Nsubiter+1;
                if A.Nsubiter~=0;
                    T1=T1*A.RT;
                end
            end
            %% Plotting Misfit History before the next subiter
            [G.F1(1).XData,G.F1(2).XData,G.F2.XData]=deal(ER.kk);
            G.F1(1).YData=ER.err;
            G.F1(2).YData=ER.err1;
            G.F2.YData=ER.temer;
            G.F1(1).Parent.XLimMode='Auto';
            drawnow
        end
  
    else
        
%% Local Methods: Simplex Downhill & Interior point
        if handles.STOP.UserData==1
            break;
        end
        [aopt(:,:),optimValues,BestL2f,ER,HVopt,DCopt,Cont,TfinL,Tfin2L]=Optimize(A,C,1,G,ER,handles,Tinv,Cont);
        Tfin2=Tfin2+Tfin2L; %Tiempo final de cálculo acumulativo (exluye las tomas de decisiones de los algoritmos locales)
        Tfin=Tfin+TfinL;% Tiempo de cálculo para el algoritmo de Optimize (excluye las tomas de decisiones y algunas acciones gráficas, poco mayor que la suma de cálculos directos)
        if ~ZOMBIE.ON,questdlg(optimValues.message,' ','Ok','Ok');end
        A.L2optimo=BestL2f; % mejor modelo encontrado
    end
    
    %% Graphing final models and curves
    [BTA,ALFA,~,E]=MODELO(aopt,GetV('limdd'));
    %Graficando modelos con el mejor ajuste (Curvas HV, DC y perfiles de velocidad
    [G.ptlmod(2).XData,G.ptlmod(2).YData,G.ptlmod(1).XData,G.ptlmod(1).YData]=...
        deal(ALFA,E,BTA,E);
    if  ~isempty(C.HVFobs)
        G.phv.YData=HVopt;
    end
    if  ~isempty(C.DCFobs)
        G.pdc.YData=DCopt;
    end
    handles.textmisfit.String=['# Evaluated models',' ',num2str(Cont.all),' ','Min. Misfit',' ',num2str(A.L2optimo)];
    [G.F1(1).XData,G.F1(2).XData,G.F2.XData]=deal(ER.kk);
    [G.F1(1).YData,G.F1(2).YData,G.F2.YData]=deal(ER.err,ER.err1,ER.temer);
    G.F1(1).Parent.XLimMode='Auto';
    drawnow
    
    %% Concatenation of data
    %Concatenacion de los datos listos para otro tipo de inversion.
    indOK=find(ER.L2ALL~=0 & ~isinf(ER.L2ALL));
    Cont.iter=length(find(~isnan(ER.err) & ~isinf(ER.err)));
    A.aini=aopt;
    Cont.valid=length(indOK);
    ER.HVALL=ER.HVALL(:,indOK);
    ER.DCALL=ER.DCALL(:,indOK);
    ER.MODALL=ER.MODALL(:,:,indOK);
    ER.L2ALL=ER.L2ALL(indOK);

    %% Continuation of inversion with a different method
    if handles.STOP.UserData==0
        if ZOMBIE.ON
            choice='No';
        else
            choice = questdlg('Continue Inversion?',' ','Yes','No','Yes');
        end
        switch choice
            case 'Yes'
                resp=1;
                handles.STOP.UserData=0;
                %NN=num2str(A.NN);
                AA=selectmethod;
                uiwait(AA)
                Tinv=evalin('base','Tinv');
                % settings depending of inversion type:
                if Tinv<0, % The user closed method selection window. 
                           % End inversion and keep the last method used as selected method
                    resp=0;
                    LoadV('Tinv',-Tinv);% The last method used is restored                 
                    %A.Niter=100;
                elseif Tinv==1;% Monte Carlo Sampling (Metropoli)
                    NUM=str2double(inputdlg({'Number of iterations','Perturbation range (%)'},'',1,{'100','5'}));
                    if ~isempty(NUM)
                        NUM(3)=NUM(1);
                        NUM(1)=0;
                        NUM(4)=NUM(2);
                        NUM(2)=0;                        
                        resp=1;
                        Perr(2).YLabel.String=' ';
                        Perr(2).YTickMode='manual';
                        Perr(2).YTick=[];                     
                    else
                        resp=0;
                    end
                elseif Tinv==2; % SA
                    prompt={['M',char(225),'rkov`s chain length'],'Number of temperatures',['Last M',char(225),'rkov`s chain length'],'Perturbation range (%)'};
                    NUM=str2double(inputdlg(prompt,'',1,{'100','1','100','5'}));
                    if ~isempty(NUM)
                        resp=1;
                        NUM(2)=max(NUM(2),1); 
                        Perr(2).YLabel.String='Temperature';
                        Perr(2).YTickMode='auto';
                        if A.T0==1
                            [Tini,T1]=deal(-A.EA*A.L2optimo/log(A.AP));
                        else
                            [Tini,T1]=deal(A.AP);
                        end
                    else
                        resp=0;
                    end                    
                elseif Tinv==3; % Modified Simulated Annealing (MSA)
                    prompt={'Number of iterations','Number of Reheatings','Number of Last iterations','Perturbation range (%)'};
                    NUM=str2double(inputdlg(prompt,' ',1,{'100','1','100','5'}));
                    if ~isempty(NUM)
                        [resp,Perr(2).YLabel.String, Perr(2).YTickMode]=deal(1,'Temperature','auto');
                        if A.T0==1
                            [Tini,T1]=deal(-A.EA*A.L2optimo/log(A.AP));
                        else
                            [Tini,T1]=deal(A.AP);
                        end
                        if A.open==0
                            choice = questdlg(sprintf('The inversion method is designed to run in parallel.\n\nTo use this method in parallel?'));
                            if strcmp('Yes',choice);
                                pop=GetV('NUMW');
                                NUMW=str2double(cell2mat(inputdlg({sprintf('Enter number of Workers (Labs): (>=2 & <512)')},'',1,{num2str(pop)})));
                                if isnan(NUMW) || NUMW<2
                                else
                                    %% Parallel settings & GUI
                                    mss=msgbox({sprintf('The Parallelization Settings will take a few seconds.\n\nPlease Wait')},'HV-inv');
                                    delete(findobj(mss,'string','OK'));
                                    delete(findobj(mss,'style','frame'));
                                    set(mss, 'CloseRequestFcn','')
                                    %edit GUI
                                    [handles.Poff.Checked,handles.edit1.Enable,handles.TDC.Enable,...
                                        handles.THV.Enable,handles.FHV.Enable,handles.THVDC.Enable,...
                                        handles.PS.Enable,handles.Models.Enable,handles.Inversion.Enable,...
                                        handles.LP.Enable,handles.SPara.Enable,handles.START.Enable,...
                                        handles.STOP.Enable,handles.Menu.Enable]=deal('off');
                                    % Enable kernels (workers)
                                    parallel(NUMW);
                                    % Edit GUI
                                    [handles.textNumwork.Visible,handles.Pon.Checked,handles.TDC.Enable,...
                                        handles.THV.Enable,handles.FHV.Enable,handles.THVDC.Enable,...
                                        handles.PS.Enable,handles.Models.Enable,handles.Inversion.Enable,...
                                        handles.LP.Enable,handles.SPara.Enable,handles.START.Enable,...
                                        handles.STOP.Enable,handles.Menu.Enable,handles.edit1.Enable]=deal('on');
                                    [handles.textNumwork.String,A.open,A.flag,NUMW]=...
                                        deal(['Connection with',' ',num2str(NUMW),' ', 'Labs.'],1,zeros(NUMW,1),NUMW);
                                    LoadV('NUMW',NUMW,'open',1);
                                    delete(mss);
                                end
                            end
                        else
                            NUMW=GetV('NUMW');
                            A.flag=zeros(NUMW,1);
                        end
                    else
                        resp=0;
                    end
                elseif Tinv==4; % Downhill Simplex
                    prompt={'Maximum number of iterations allowed','Maximum number of evaluations allowed','Termination tolerance on the function value (TolFun)','Termination tolerance on parameters (TolX)'};
                    NUM=str2double(inputdlg(prompt,' ',1,{'500','500','1.e-34','1.e-34'}));
                    if ~isempty(NUM)
                        resp=1;
                        Perr(2).YLabel.String=' ';
                        Perr(2).YTickMode='manual';
                        Perr(2).YTick=[];
                    else
                        resp=0;
                    end
                elseif Tinv==5; % Interior point
                    prompt={'Maximum number of iterations allowed','Maximum number of evaluations allowed','Termination tolerance on the function value (TolFun)','Termination tolerance on parameters (TolX)'};
                    NUM=str2double(inputdlg(prompt,',',1,{'500','500','1.e-34','1.e-34'}));
                    if ~isempty(NUM)
                        [resp,Perr(2).YLabel.String,Perr(2).YTickMode,Perr(2).YTick]=deal(1,' ','manual',[]);
                    else
                        resp=0;
                    end
                end
            otherwise % User pressed 'No more inversions'
                 resp=0;           
        end
        
        %% reassignment of settings for inversion type
        if resp==1
            A.aini=aopt;
            A.Niter=NUM(1);
            A.Nsubiter=NUM(2);
            A.Niterf=NUM(3);
            A.NN=NUM(4);           
            %Remove data for statistics (standard errors)
            if isfield(DATERROR, 'mod')
                DATERROR = rmfield(DATERROR,'mod');
                DATERROR = rmfield(DATERROR,'Mis');
                DATERROR = rmfield(DATERROR,'temp');
            end
            iForStats=0;% Inicia contador del misfit para las ultimas iteraciones (estadísticas) ¿¿¿necesario si vamos de MC a MC???
            Cont.IterLastMethod=0;% Resetea el contador de iteraciones del último metodo, para usarse en el report final
        else
            %Enable buttons of GUI
            [handles.STOP.Visible,handles.smodel.Visible,handles.START.Visible]= deal('off','on','on');
            % Message: minimum misfit
            if ~ZOMBIE.ON
                msgbox(['Finished program. ',' The minimun misfit is ' num2str(A.L2optimo)],' ');
            end
        end
    end
end

%% End of inversion
% Assigning variables in the workspace
if ~isempty(ini1), %Agregando datos del modelo inicial 
    ER.MODALL(:,:,end+1)=modini1;
    ER.HVALL(:,end+1)=HVTini1;
    ER.DCALL(:,end+1)=DCTini1;
    ER.L2ALL(end+1)=L2ini1;
end

%% Valid models
acceptedmodels=(find(~isnan(ER.L2ALL)));
A.aini=aopt; % Modelo optimo
% Assigning variables in the workspace
LoadV('aopt',aopt,'A',A,'HVALL',ER.HVALL(:,acceptedmodels),'DCALL',ER.DCALL(:,acceptedmodels),...
    'MODALL',ER.MODALL(:,:,acceptedmodels),'L2ALL',ER.L2ALL(acceptedmodels),'DATAERROR',DATERROR);   
%% Reports of inversion
if Tinv==1
    %Final report for Metropoli Walk
    if pass==1
        ST1='Monte Carlo';
        handles.Report.Enable= 'on';
        handles.Report.String=[{'Inversion report'};{' '};{['Inversion Method : ',ST1 ]};{['Number of iterations : ',num2str(Cont.IterLastMethod)]};...
            {['Perturbation range : ',num2str(A.NN),'%']};{' '};...
            {['Total iterations : ',num2str(Cont.iter)]};{['Total models : ',num2str(Cont.all)]};{['Accepted models  : ',num2str(length(acceptedmodels))]} ;{['Models for statistics : ',num2str(iForStats)]};{'  '};...
            {['Maximum misfit : ',num2str(max(ER.L2ALL))]};{['Minimum misfit : ',num2str(min(ER.L2ALL))]};{'  '};...
            {['Total time (s): ',sprintf('%.2f',Tfin2)]};{['Avg time per model (s): ',sprintf('%.2f',Tfin/Cont.all)]}];
    end
elseif Tinv==2
    ST1='SA';
    KGF=['Number of temperatures : ',num2str(A.Nsubiter)];
    KGF2=['Temp. drop / reheatings : ',num2str(A.ST)];    
elseif Tinv==3
    ST1='Modified SA';
    KGF=['Number of reheatings : ',num2str(A.Nsubiter)];
    KGF2=['Temp. drop / reheatings : ',num2str(A.ST)];
end
if Tinv==2||Tinv==3
    if pass==1
        %Final report for SA and MSA
        handles.Report.Enable= 'on';
        handles.Report.String=[{'Inversion report'};{' '};{['Inversion Method : ',ST1 ]};{['Number of Iterations : ' ,num2str(Cont.IterLastMethod)]};{KGF};{['Models for statistics : ',num2str(iForStats)]};...
            {['Initial Temperature : ',num2str(Tini)]};{['Final Temperature : ',num2str(T1)]};{['Temperature ratio : ', num2str(A.RT)]} ;{['Perturbation range : ',num2str(A.NN),'%']};{' '};...
            {['Total iterations : ',num2str(Cont.iter)]};{KGF2};{['Evaluated models : ',num2str(Cont.all)]};{['Valid models  : ',num2str(length(acceptedmodels))]} ;{'  '};...
            {['Maximum misfit : ',num2str(max(ER.L2ALL))]};{['Minimum misfit : ',num2str(min(ER.L2ALL))]};...
            {['Total time (s): ',sprintf('%.2f',Tfin2)]};{['Avg time per model (s): ',sprintf('%.2f',Tfin/Cont.all)]}];
    end
end
%Final report for local methods
if Tinv==4
    ST1='Downhill simplex ';
elseif Tinv==5
    ST1='Interior-point';
end
if pass==1 && (Tinv==4 || Tinv==5);
    handles.Report.Enable= 'on';
    handles.Report.String=[{'Inversion report'};{' '};{['Inversion Method : ',ST1 ]};{['Number of iterations : ',num2str(Cont.IterLastMethod)]};{'  '};...
        {['Total iterations : ',num2str(Cont.iter)]};{['Valid models : ',num2str(length(acceptedmodels))]};{['Evaluated models  : ',num2str(Cont.all)]} ;{'  '};...
        {['Maximum misfit : ',num2str(max(ER.L2ALL))]};{['Minimum misfit : ',num2str(min(ER.L2ALL))]};{'  '};...
        {['Total time : ',sprintf('%.2f',Tfin2)]};{['Avg time per model (s): ',sprintf('%.2f',Tfin/Cont.all)]}];
end

