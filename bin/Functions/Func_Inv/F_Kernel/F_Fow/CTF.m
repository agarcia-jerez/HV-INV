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

function[a,HVT,DCT,out,L2f]=CTF(A,C,jj,typ,~)
% This function generates a random model (to be stored in ouput "a"), calls
% HVf.exe for computation of its forward H/V and/or dispersion curve
% (HVT,DCT) and the corresponding misfit L2f.
%
% out is an auxiliary variable related with pretubation directios used by
% the Modified S A method only (experimental algorithm by J. Piña)
%
% Inputs:
% A contains problem settings
% C contains experimental curves
% jj an integer to identify the call number to this function
% typ controls the method used for generation of the random model
% typ = 1 for random models fulfilling global ranges/constraints (Random
%         Search)
% typ = 2 Take a perturbation range around the model A.aini as additional
%         constraint (SA, MSA & MC sampling methods)
% typ = 3 Use particular model from ini1 (base WS) instead of random model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Select either \ or / depending on the computer OS (Windows-Linux-Mac)
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

%% Initialize some variables
%L2f=1; %Misfit
[FC,L2f]=deal(1);% it will be set to 0 after a successfully evaluated misfit. Otherwise, alternative models can be tried
[DCT,HVT]=deal(zeros(1,C.nmDC),zeros(1,C.nmHV)); %Dispersion curve; %H/V curve

if typ==2
    flag=A.flag(jj); % Only for MSA (experimental Piña's version of Simulated Annealing).
    % If flag = 1 the direction of the model variations are fixed in the current stage depending on the
    % behavior in the previous steps.For SA and MC, flag is always zero.
end

while FC==1, % While we don't have a successfull misfit evaluation
    FC=0;
    %% Generation of random models
    if typ==1
        %% Random models in the full ranges (encoded in A.var) with additional constraints (depending on A.ZLVs,A.ZLVp,A.HHS)
        % A.pol,A.limseg,A.polHHS,A.facHHS,A.ProbLimSupSegHHS,A.LimVelSegHHS are inner variables of the random model
        % generation algorithm, exported to improve performance.
        a=ini1_plus(A.NC,A.var,A.ZLVs,A.ZLVp,A.pol,A.limseg,A.HHS,A.polHHS,A.facHHS,A.ProbLimSupSegHHS,A.LimVelSegHHS,A.ISDEPTH);
        out=0;% not used for this type
    elseif typ==2
        %% Random models around A.aini with preturbation range A.NN1. The full ranges (encoded in A.var) and the additional
        %  constraints (A.ZLVs,A.ZLVp,A.HHS) are also fulfilled.
        [a,out]=ini2A(A.aini,A.NC,A.var,A.NN1,A.ZLVs,A.ZLVp,flag,A.in,A.HHS,A.ISDEPTH);
        if ~isempty(find(a(1:end-1,:)==A.aini(1:end-1,:),1)) || ~isempty(find(a(end,2:end)==A.aini(end,2:end),1))
            % Some property colapsed to a boundary in MSA (Piña). Release the constraints about directions of parameter variations and try again.
            [a,out]=ini2A(A.aini,A.NC,A.var,A.NN1,A.ZLVs,A.ZLVp,0,A.in,A.HHS,A.ISDEPTH);
        end
    elseif typ==3
        a=GetV('ini1');% take initial model
        out=0;% not used for this type
    end
    
    %% Writting model file
    ai=num2str(jj);
    mod=['etc',barra,ai,'mod.txt'];
    dlmwrite(mod,A.NC);            
%    if A.ISDEPTH
%        dlmwrite(mod,[[diff([0;a(1:end-1,1)],1,1);0],a(:,2:end)],'-append','delimiter','\t','precision',8);
%    else
        dlmwrite(mod,a,'-append','delimiter','\t','precision','%8.2f');
%    end
    arch=[' -f ',mod];
    
    %% Forward calculation of dispersion curve if required
    if  ~isempty(C.DCFobs)
        % Write frequency file
        ffrec=['etc',barra,ai,'freDC.txt'];
        dlmwrite(ffrec,C.DCFobs,'delimiter','\t','precision',8);
        %Number of mode (Only fundamental mode)
        NM=num2str(C.MODE+1);
        % Generation of inline code for HVf.exe forward calculation program. Dispersion curves
        if C.POL==1 && C.VEL==1
            gr=['exe',barra,folder,barra,'HVf',' -nmr ',NM,' -ff ',ffrec, ' -ph ', arch ];% RALEIGH PHASE VELOCITY
        elseif C.POL==1 && C.VEL==2
            gr=['exe',barra,folder,barra,'HVf',' -nmr ',NM,' -ff ',ffrec, ' -gr ', arch ]; % RAYLEIGH GROUP VELOCITY
        elseif C.POL==2 && C.VEL==1
            gr=['exe',barra,folder,barra,'HVf',' -nml ',NM,' -ff ',ffrec, ' -ph ', arch ]; % LOVE PHASE VELOCITY
        else
            gr=['exe',barra,folder,barra,'HVf',' -nml ',NM,' -ff ',ffrec, ' -gr ', arch ]; % LOVE GROUP VELOCITY
        end
        % Run HVf.exe for synthetic H/V
        [~,OUT]=system(gr);
        % Read forward computation from HVf.exe output
        DCC = sscanf(OUT,'%f',[1 inf] )';
        DC=DCC(3:end);
        if length(DC)==length(C.DCVobs)
            % Condition of stability
            if isnan(sum(DC)) || isinf(sum(DC)) ;
                imal=isnan(DC) | isinf(DC);
                if length(imal)>length(DC)*0.1
                    FC=1; % invalid model
                    flag=0;% For MSA method, release constraints about variation directions to favour valid model finding
                    if typ==3
                        FC=0;
                        [DC,L2f]=deal(nan);
                    end
                else
                    DC(imal)=interp1(C.DCFobs(~imal),DC(~imal),C.DCFobs(imal),'linear','spline');
                end
            else
                FC=0;
            end
            % Another quality control
            if find(DC==0)
                FC=1;% invalid model
                flag=0;% For MSA method, release constraints about variation directions to favor valid model finding
                if typ==3
                    FC=0;
                    [DC,L2f]=deal(nan);
                end
            end
            % Dispersion curve output
            DCT=1./DC;
        else
            FC=1;
            flag=0;% For MSA method, release constraints about variation directions to favor valid model finding
            if typ==3
                FC=0;
                [DCT,L2f]=deal(nan);
            end
        end
    end
    
    %% Forward calculation of H/V curve if required
    if  ~isempty(C.HVFobs) && ~isnan(L2f);
        % Write frequency file
        ffrec=['etc',barra,ai,'freHV.txt'];
        dlmwrite(ffrec,C.HVFobs,'delimiter','\t','precision',8);
        % Number of modes
        %nmR=; %Rayleigh
        %nmL=num2str(C.NL); %Love
        % Run HVf.exe for computation of H/V
        [~,HVC]=system(['exe',barra,folder,barra,'HVf',' -nmr ',num2str(C.NR),' -nml ',num2str(C.NL),' -nks ',num2str(C.dk),' -apsv ',num2str(A.apsv),' -ash 0.01 ',' -hv ',' -ff ',ffrec,arch]);
        % Extract H/V from program output
        OUT = sscanf(HVC,'%f',[2 inf] )';
        HVTT=OUT(:,2);
        if length(HVTT(:,1))==length(C.HVAobs)
            % Check quality and interpolate for missing frequencies
            if isnan(sum(HVTT))|| isinf(sum(HVTT)) ;
                imal=isnan(HVTT) | isinf(HVTT);
                if length(imal)>length(HVTT)*0.1
                    FC=1;% Invalid model
                    flag=0;% For MSA method, release constraints about variation directions to favor valid model finding
                    if typ==3
                        FC=0;
                         [HVTT,L2f]=deal(nan);
                    end
                else
                    HVTT(imal)=interp1(C.HVFobs(~imal),HVT(~imal),C.HVFobs(imal),'linear','spline');
                end
            else
                FC=0;
            end
            % Another quality control
            if find(HVTT==0);
                FC=1;% invalid model
                flag=0;% For MSA method, release constraints about variation directions to favor valid model finding
                if typ==3
                    FC=0;
                     [HVTT,L2f]=deal(nan);
                end
            end
            % Output
            HVT=HVTT;
        else
            FC=1;
            flag=0;% For MSA method, release constraints about variation directions to favor valid model finding
            if typ==3
                FC=0;
                 [HVT,L2f]=deal(nan);
            end
        end
    end
    
    %% Compute misfit
    if FC==0 && ~isnan(L2f)
        % Joint inversion case
        if  ~isempty(C.HVFobs) &&  ~isempty(C.DCFobs)
            L2f=2*(1-C.ww)*sum(((HVT-C.HVAobs).^2)./C.HVDS.^2) + 2*C.ww*sum(((DCT-C.DCVobs).^2)./C.DCDS.^2) ;
            % H/V case
        elseif  ~isempty(C.HVFobs) &&  isempty(C.DCFobs)            
            if ~isempty (C.Cov)
%                tic
%               L2f=(C.HVAobs-HVT).'*(C.Cov\((C.HVAobs-HVT)));
%                a1=toc;
%                tic
                L2f=(C.HVAobs-HVT).'*C.CovInv*(C.HVAobs-HVT);
%                a2=toc;
%                disp([a1 a2 L2f L2f2])
            else
                L2f=sum(((HVT-C.HVAobs).^2)./C.HVDS.^2);    
            end
            
            % Dispersion curve case
        elseif  isempty(C.HVFobs) && ~isempty(C.DCFobs)
            L2f=sum(((DCT-C.DCVobs).^2)./C.DCDS.^2);
        end
    end
end
