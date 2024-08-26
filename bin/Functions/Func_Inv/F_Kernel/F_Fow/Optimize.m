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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function applies local inversion methods: Simplex down-hill (Method #4)
% and Interior point (Method #5).
% It prepares input files, runs HVf.exe and updates the GUI

function[Bestmodel,optimValues,BestL2f,ER,BestHVT,BestDCT,Cont,TfinL,Tfin2L]=Optimize(A,C,jj,G, ER,handles,method,Cont)
%INPUT
% A Iteration Information
% C Curves information (Frequencies, velocities, amplitud)
% jj % Number of control files
% G  Structure for plotting H/V and DC
% ER  Structure for plotting misfit and iterations
% handles  Elements of  control GUI
% method: Method for optimization to be used: Simplex down-hill (Method #4) or Interior point (Method #5)
% Cont % Contador de iteraciones y modelos.
%OUTPUT
%Bestmodel: The best model
%optimValues %Id for stop inversion
% BestL2f Best misfit
% ER% Struc for model and misfit history
% BestHVT % Best curve HV
% BestDCT % Best curve DC
% Cont Contador de modelos e iteraciones
% TfinL,% time of call to forward calculation
% Tfin2L % total time, includes TfinL and updating of some figures

[TfinL,Tfin2L]=deal(0); %tiempo para estadistica;
Cont.IterLastMethod=Cont.iter;% backup of the total number of iterarions in the begining
if ispc
    barra='\';
    folder='exe_Win';
elseif ismac
    barra='/';
    folder='exe_Mac';
elseif isunix
    barra='/';
    folder='exe_Linux';
end
if method==4,% Simplex down-hill (fminsearch). In this case, the unknowns are those particular model parameters
    % for which maximum value is strictly larger than the minimum
    indexlist=[A.var(:,3);A.var(:,6);A.var(:,9);A.var(:,12)]>0;%Each row represent a layer. Four column representing thickness/depth, Vp, Vs, density.
    %Set actual variables (true) are those with positive increments (maximum - minimum)
elseif method==5 % Interior point (fmincon): In this case we remove from the list of unknowns those parameters which are constant along the entire structure
    % If a parameter is constant for a subset of layers then it will be treated as a constraint for fmincon
    indexlist=logical([]);
    % Accepts a property as a variable (set its position as true) if it is a variable for some layers.
    if find(A.var(:,3)>0,1), indexlist=[indexlist;true(A.NC-1,1);false];else indexlist=[indexlist;false(A.NC,1)];end% Thickness
    if find(A.var(:,6)>0,1), indexlist=[indexlist;true(A.NC,1)];else indexlist=[indexlist;false(A.NC,1)];end % Vp
    if find(A.var(:,9)>0,1), indexlist=[indexlist;true(A.NC,1)];else indexlist=[indexlist;false(A.NC,1)];end % Vs
    if find(A.var(:,12)>0,1),indexlist=[indexlist;true(A.NC,1)];else indexlist=[indexlist;false(A.NC,1)];end % rho
end

% Generates 'a' (current full model, internal variable) and x0 ("a" restricted to actual variables).
% Note that both of them may list thicknesses or depths
if A.ISDEPTH
    % copy of A.aini to store the full optimum model. Elements of positions "indexlist" will be updated later.
    % Conversion from thickness to depth is required within this function
    a=[[round(cumsum(A.aini(1:end-1,1)),2);0],A.aini(:,2:4)];
else
    % copy of A.aini to store the full optimum model. Elements of positions "indexlist" will be updated later.
    a=A.aini;
end
x0=a(indexlist);% Take the actual variables of the initial model

%% Below, conditions over thickness, Vp, Vs and density are coded. They involve all these model parameter for the moment (not only the actual unknowns).

% Take lower bound and upper bound (Poisson's ratio excluded)
lb=[A.var(:,1);A.var(:,4);A.var(:,7);A.var(:,10)];
ub=[A.var(:,2);A.var(:,5);A.var(:,8);A.var(:,11)];

% Take conditions of downward-increasing variables
CREC=zeros(0,4*A.NC); % Empty array. (It will support evaluation of CREC(:,indexlist) even if CREC remains empty)
if ~A.ZLVp ||~A.ZLVs||(A.ISDEPTH && A.NC>2), % In this case Vp or Vs or layer depths are assumed to increase as depth increases
    MASK=eye(A.NC,A.NC-1);
    MASK(2:A.NC+1:end)=-1;
    MASK=MASK';
    if ~A.ZLVp,CREC=[CREC;[zeros(1*A.NC,A.NC-1);MASK';zeros(2*A.NC,A.NC-1)]'];end% Includes "Vp increases downward" in the matrix condition CREC*x <= 0
    if ~A.ZLVs,CREC=[CREC;[zeros(2*A.NC,A.NC-1);MASK';zeros(1*A.NC,A.NC-1)]'];end% Includes "Vs increases downward" in the matrix condition CREC*x <= 0
    if(A.ISDEPTH && A.NC>2)
        % If we use depth instead of thickness we also code that depth increases as the layer index increases:
        % We can reduce the number of conditions excluding trivial cases
        % For depths, the A.NC-1-th row of MASK (the last one) makes no
        % sense. We have  A.NC-2 conditions maximum
        NonTrivial=A.var(1:end-2,2)>A.var(2:end-1,1);
        CREC=[CREC;[MASK(NonTrivial,:)';zeros(3*A.NC,sum(NonTrivial))]'];
    end
end

HARDSEMI=zeros(0,4*A.NC);
if A.ZLVs && A.HHS, %If low S-wave vel. zones are allowed, the halfspace can still be defined as the layer for which Vs is maximum (HARDSEMI*x <= 0)
    MASK=eye(A.NC-1,A.NC);MASK(:,A.NC)=-1;
    %HARDSEMI=[HARDSEMI;[zeros(1*A.NC,A.NC-1);MASK';zeros(2*A.NC,A.NC-1)]'];% Halfspace has the highest Vp
    HARDSEMI=[HARDSEMI;[zeros(2*A.NC,A.NC-1);MASK';zeros(1*A.NC,A.NC-1)]'];% Halfspace has the highest Vs
    %HARDSEMI=[HARDSEMI;[zeros(3*A.NC,A.NC-1);MASK']'];% Halfspace has the highest density
end

% We write the condition of maximum and minimum Poisson's ratio as POIMAX*x<=0
[POIMIN,POIMAX]=deal(zeros(A.NC,4*A.NC));
[POIMIN(A.NC*A.NC+1 : A.NC+1:   A.NC*A.NC+1 + (A.NC-1)*(A.NC+1)),POIMIN(2*A.NC*A.NC+1 : A.NC+1: 2*A.NC*A.NC+1 + (A.NC-1)*(A.NC+1)),...
    POIMAX(  A.NC*A.NC+1 : A.NC+1:   A.NC*A.NC+1 + (A.NC-1)*(A.NC+1)),POIMAX(2*A.NC*A.NC+1 : A.NC+1: 2*A.NC*A.NC+1 + (A.NC-1)*(A.NC+1))]=...
    deal(-1,sqrt((2*A.var(:,13)-2)./(2*A.var(:,13)-1)),1,-sqrt((2*A.var(:,14)-2)./(2*A.var(:,14)-1)));
%%
% Characteristics of the current best model
BestHVT=[];
BestDCT=[];
if A.ISDEPTH
    % copy of A.aini to store the full optimum model. Elements of positions "indexlist" will be updated later.
    % Conversion from thickness to depth is required within this function
    Bestmodel=[[round(cumsum(A.aini(1:end-1,1)),2);0],A.aini(:,2:4)];
else
    Bestmodel=A.aini;
end
BestL2f=inf;

% Enlarge the lists of tested model characteristics in order to store A.Niter more models
ER.DCALL=[ER.DCALL zeros(C.nmDC,A.Niter)];
ER.HVALL=[ER.HVALL zeros(C.nmHV,A.Niter)];
ER.MODALL=cat(3,ER.MODALL, zeros(A.NC,4,A.Niter*2));
ER.L2ALL=[ER.L2ALL ;zeros(A.Niter*2,1)];
ER.err=[ER.err; nan(A.Niter,1)];
ER.err1=[ER.err1; nan(A.Niter,1)];
ER.temer=[ER.temer; nan(A.Niter,1)];
ER.kk=[ER.kk ;nan(A.Niter,1)];

% Get limdd. The maximum depth for representation will be 1.1*limdd
limdd=GetV('limdd');

if method==4,% Call Simplex-downhill method
    OPTIONS = optimset('MaxIter',A.Niter,'MaxFunEvals',A.Nsubiter,'OutputFcn', @outfun,'TolFun',A.Niterf,'TolX',A.NN);
    LIN_INEQ=[CREC;HARDSEMI;POIMIN;POIMAX];
    [x,~,~,optimValues] = fminsearch(@misfitfun,x0,OPTIONS);
elseif method==5, % Call Interior-Point method
    % In this case we restrict the conditions to the actual unknowns (true in indexlist matrix).
    % Only complete columns in ub, lb or LIN_INEQ can be removed
    ub=ub(indexlist);lb=lb(indexlist);
    LIN_INEQ=[CREC(:,indexlist);HARDSEMI(:,indexlist);POIMIN(:,indexlist);POIMAX(:,indexlist)];
    % Set options for fmincon
    options = optimoptions(@fmincon,'MaxFunEvals',A.Nsubiter,'MaxIter',A.Niter,'TolFun',A.Niterf,'TolX',A.NN,'TolCon',1e-50,'TypicalX',x0,...
        'UseParallel',false,'FinDiffRelStep',1e-5,'FinDiffType','central','OutputFcn',@outfun); % add 'Diagnostics','on' for debug
    [x,~,~,optimValues] = fmincon(@misfitfun,x0,LIN_INEQ,zeros(size(LIN_INEQ,1),1),[],[],lb,ub,[],options);
end
a(indexlist)=x; % Raplace optimum values of the unknowns in the whole model a


    function L2f = misfitfun(x)
        % Computes the misfit of model "a" after replacing the unknowns at positions "indexlist" with x. This is a nested function for Optimize.
        tic;
        FC=1;% 0 = successful calculation; 1 = errors in forward calculation
        DCT=zeros(C.nmDC,1);
        HVT=zeros(C.nmHV,1);
        a(indexlist)=x;
        L2f=inf;
        if method==5  || ( sum(max(LIN_INEQ*a(:),0))==0 && sum(a(:)<lb|a(:)>ub)==0 ),
            % sum(max(LIN_INEQ*a(:),0)) is the number of violated linear inequalities
            % sum(a(:)<lb|a(:)>ub) is the number of not fulfilled bounds
            % Here we are impose a drastic rejection of invalid models for the
            % simplex-downhill method, taking the default misfit L2f=inf for them
            
            %% Writing imput files for hvf.exe.
            % Model file
            ai=num2str(jj);
            modl=['etc',barra,ai,'mod.txt'];
            dlmwrite(modl,A.NC);
            if A.ISDEPTH
                dlmwrite(modl,[[a(1,1);diff(a(1:end-1,1));0],a(:,2:end)],'-append','delimiter','\t','precision','%8.2f');
            else
                dlmwrite(modl,a,'-append','delimiter','\t','precision','%8.2f');
            end
            
            if ~isempty(C.DCFobs), % Computation of dispersion curve is required
                
                ffrec=['etc',barra,ai,'freDC.txt'];% file of frequencies
                dlmwrite(ffrec,C.DCFobs,'delimiter','\t','precision',8);
                NM=num2str(C.MODE+1);
                
                %% Preparing command line for the OS
                if C.POL==1 && C.VEL==1
                    gr=['exe',barra,folder,barra,'HVf',' -nmr ',NM,' -ff ',ffrec, ' -ph -f ',modl];
                elseif C.POL==1 && C.VEL==2
                    gr=['exe',barra,folder,barra,'HVf',' -nmr ',NM,' -ff ',ffrec, ' -gr -f ',modl];
                elseif C.POL==2 && C.VEL==1
                    gr=['exe',barra,folder,barra,'HVf',' -nml ',NM,' -ff ',ffrec, ' -ph -f ',modl];
                else
                    gr=['exe',barra,folder,barra,'HVf',' -nml ',NM,' -ff ',ffrec, ' -gr -f ',modl];
                end
                
                %% Running HVf.exe
                [~,OUT]=system(gr);
                
                %% Load dispersion curve computed by HVf.exe
                DCC = sscanf(OUT,'%f',[1 inf] )';
                DC=DCC(3:end);
                if sum(DC)~=0
                    if length(DC)==length(C.DCVobs)
                        if isnan(sum(DC)) || isinf(sum(DC))
                            imal=isnan(DC) | isinf(DC);
                            if length(imal)>length(DC)*0.1
                                FC=1; %Dispersion curve could not be computed
                            else
                                DC(imal)=interp1(C.DCFobs(~imal),DC(~imal),C.DCFobs(imal),'linear','spline'); % DC at a small number of frequencies can be fixed
                            end
                        else
                            FC=0;
                        end
                        if find(DC<=0),FC=1;end
                        DCT=1./DC;% from slowness to velocity
                    else
                        FC=1;
                    end
                else
                    FC=1;
                end
            end
            
            if ~isempty(C.HVFobs) % Computation of H/V is required
                ffrec=['etc',barra,ai,'freHV.txt'];
                dlmwrite(ffrec,C.HVFobs,'delimiter','\t','precision',8);
                %% Running HVf.exe
                [~,HVC]=system(['exe',barra,folder,barra,'HVf',' -nmr ',num2str(C.NR),' -nml ',num2str(C.NL),' -nks ',num2str(C.dk),' -apsv ',num2str(A.apsv),' -ash 0.01 ',' -hv ',' -ff ',ffrec,' -f ',modl]);
                
                %% Load H/V curve computed by HVf.exe
                OUT = sscanf(HVC,'%f',[2 inf] )';
                HVTT=OUT(:,2);
                if sum(HVTT)~=0
                    if length(HVTT(:,1))==length(C.HVAobs)
                        if isnan(sum(HVTT))|| isinf(sum(HVTT)) ;
                            imal=isnan(HVTT) | isinf(HVTT);
                            if length(imal)>length(HVTT)*0.1
                                FC=1;
                            else
                                HVTT(imal)=interp1(C.HVFobs(~imal),HVT(~imal),C.HVFobs(imal),'linear','spline');
                            end
                        else
                            FC=0;
                        end
                        if find(HVTT<=0) ;FC=1;end
                        HVT=HVTT;
                    else
                        FC=1;
                    end
                else
                    FC=1;
                end
            end
            
            if FC==0,% Compute Misfit. Otherwise we keep the default infinite misfit (L2f=inf)
                if ~isempty(C.HVFobs) && ~isempty(C.DCFobs)
                    L2f=2*(1-C.ww)*sum(((HVT-C.HVAobs).^2)./C.HVDS.^2) + 2*C.ww*sum(((DCT-C.DCVobs).^2)./C.DCDS.^2) ;% Joint inversion
                elseif ~isempty(C.HVFobs)
                    if ~isempty (C.Cov)
                        %                         tic
                        %                         L2f=(C.HVAobs-HVT).'*(C.Cov\((C.HVAobs-HVT)));
                        %                         a1=toc;
                        %                         tic
                        L2f=(C.HVAobs-HVT).'*C.CovInv*(C.HVAobs-HVT);%log(C.det)+
                        %                         a2=toc;
                        %                         disp([a1 a2 L2f L2f2])
                    else
                        L2f=sum(((HVT-C.HVAobs).^2)./C.HVDS.^2);
                    end
                    
                elseif  ~isempty(C.DCFobs)
                    L2f=sum(((DCT-C.DCVobs).^2)./C.DCDS.^2);
                end
            end
            TfinL=TfinL+toc;
        end
        
        if L2f<BestL2f, % In this case, state the current model as the best model.
            if ~isempty(C.HVFobs)
                BestHVT=HVT;
                G.phv.YData=BestHVT;
            end
            if ~isempty(C.DCFobs)
                BestDCT=DCT;
                G.pdc.YData=BestDCT;
            end
            BestL2f=L2f;
            Bestmodel=a;
            % Updates the represented best model
            if A.ISDEPTH
                H(2:2:2*size(a,1))=[a(1:end-1,1);1.1*limdd]';
                H(3:2:end-1)=H(2:2:end-2)';
            else
                H(2:2:2*size(a,1))=[round(cumsum(a(1:end-1,1)),2);1.1*limdd]';
                H(3:2:end-1)=H(2:2:end-2)';
            end
            VEL(2:2:2*size(a,1))=a(:,2)';
            VEL(1:2:end-1)=a(:,2)';
            VEL2(2:2:2*size(a,1))=a(:,3)';
            VEL2(1:2:end-1)=a(:,3)';
            G.ptlmod(2).XData=VEL;
            G.ptlmod(1).XData=VEL2;
            [G.ptlmod(2).YData,G.ptlmod(1).YData]=deal(H);
            [G.F1(1).XData,G.F1(2).XData]=deal(ER.kk);
            G.F1(1).YData=ER.err;
            G.F1(2).YData=ER.err1;
            G.F1(1).Parent.XLimMode='Auto';
            drawnow
        end
        % Save model, save forward calculation, update message
        Cont.all=Cont.all+1;
        ER.DCALL(:,Cont.all)=DCT;
        ER.HVALL(:,Cont.all)=HVT;
        ER.L2ALL(Cont.all)=deal(L2f);
        if A.ISDEPTH
            ER.MODALL(:,:,Cont.all)=[[a(1,1);diff(a(1:end-1,1));0],a(:,2:end)];
        else
            ER.MODALL(:,:,Cont.all)=a;
        end
        set(handles.textmisfit,'String',['# Evaluated models','  ',num2str(Cont.all),'  ','Min. Misfit','  ',num2str(BestL2f)]);
        drawnow
        Tfin2L=Tfin2L+toc;
        return;
    end % misfitfun ends

    function stop = outfun(~, optimValues,~)
        % After each iteration this function is called from the inversion routine.
        % Requests the Matlab inversion algorithm to stop if the users pressed the stop button
        if get(handles.STOP,'Userdata')==1
            stop=true;
        else
            stop=false;
        end
        % Updates Misfit History panel
        if optimValues.iteration~=0 && optimValues.iteration~=1
            Cont.iter=Cont.iter+1;
            ER.kk(Cont.iter)=Cont.iter;
            [ER.err(Cont.iter),ER.err1(Cont.iter)]=deal(optimValues.fval);
        end
        
        return;
    end % outfun ends

if A.ISDEPTH
    % From depth to thickness
    Bestmodel(2:end-1,1)=diff(Bestmodel(1:end-1,1));
end
Cont.IterLastMethod=Cont.iter-Cont.IterLastMethod;% Calculation of the actual number or iterations with the last local method
return;

end% Optimize function ends



