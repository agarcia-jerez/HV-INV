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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates a random model in the specified ranges of parameter fulfilling some constraints (specified by means of variables ZLVs,ZLVp,HHS,ISDEPTH)
% Details can be found in Garcia-Jerez et al. (2016, Computers & Geosciences)

function [a, pol,limseg, polHHS,facHHS,ProbLimSupSegHHS,LimVelSegHHS]=ini1_plus(NC,var, ZLVs,ZLVp,varargin) 

% Usage:
% [a, pol,limseg, polHHS,facHHS,ProbLimSupSegHHS,LimVelSegHHS]=ini1_plus(NC,var, ZLVs,ZLVp)
% [a, pol,limseg, polHHS,facHHS,ProbLimSupSegHHS,LimVelSegHHS]=ini1_plus(NC,var, ZLVs,ZLVp,pol,limseg)    
% [a, pol,limseg, polHHS,facHHS,ProbLimSupSegHHS,LimVelSegHHS]=ini1_plus(NC,var, ZLVs,ZLVp,pol,limseg, HHS)
% [a, pol,limseg, polHHS,facHHS,ProbLimSupSegHHS,LimVelSegHHS]=ini1_plus(NC,var, ZLVs,ZLVp,pol,limseg, HHS,polHHS,facHHS,ProbLimSupSegHHS,LimVelSegHHS)
% [a, pol,limseg, polHHS,facHHS,ProbLimSupSegHHS,LimVelSegHHS]=ini1_plus(NC,var, ZLVs,ZLVp,pol,limseg, HHS,polHHS,facHHS,ProbLimSupSegHHS,LimVelSegHHS,ISDEPTH)
% [a]=ini1_plus(NC,var, ZLVs,ZLVp,___)    
% [a, pol,limseg]=ini1_plus(NC,var, ZLVs,ZLVp,___)
% [a, pol,limseg, polHHS,facHHS,ProbLimSupSegHHS,LimVelSegHHS]=ini1_plus(NC,var, ZLVs,ZLVp,pol,limseg, HHS,___)
%
% Meaning of the variables:
% NC  = Number of layers.
% var = Ranges for the parameters. For each layer (row of the matrix) it contains Min Thickness,Max Thickness, amplitude of interval in thickness,
%       min Vp, max Vp, Ampl in Vp, min Vs, max Vs, Ampl in Vs, Min density, Max density, Ampl in desity, Min Poisson ratio,Max Poisson ratio, 
%       Ampl in Poisson ratio
%
% ZLVs = Are Low Velocity Zones in Vs allowed? 1=Yes, 0=No
% ZLVp = Are Low Velocity Zones in Vp allowed? 1=Yes, 0=No
% Currently, only one of these two conditions can be set to "No" at most.
%
% HHS  = "Hard HalfSpace" Has the halfspace the higher Vs? Currently, only available if ZLVs=1 && ZLVp=1. The default is FALSE (=0)
%
% ISDEPTH=1 Obtain random depths (1) or thicknesses (0). Output model uses thickness
%
% pol & limseg. Variables used to describe probability densities related with ZLVs==0 or ZLVp==0 conditions.
% They can be saved from the output and provided in subsequent calls to the function to speed up the calculations
% (while the ranges and conditions remain unchanged)
%
% polHHS, facHHS, ProbLimSupSegHHS & LimVelSegHHS: Variables used to describe probability densities for the HHS=1 condition.
% They can be saved from the output and provided in subsequent calls to the function to speed up the calculations
% (while the ranges and conditions remain unchanged)
%

% Deal with optional inputs and outputs
na_out=nargout;
na_in=nargin;
if na_in>4
    pol=varargin{1};
    limseg=varargin{2};
end
if na_in>6
    HHS=varargin{3};
end
if na_in>7&& na_in>10
    polHHS=varargin{4};
    facHHS=varargin{5};
    ProbLimSupSegHHS=varargin{6};
    LimVelSegHHS=varargin{7};
else
    polHHS=[];
    facHHS=[];
    ProbLimSupSegHHS=[];
    LimVelSegHHS=[];        
end
if na_in>11
    ISDEPTH=varargin{8};
else
    ISDEPTH=false;% Default is thickness
end

epsil=1e-3;
if ~isempty(find(diff(var(:,7))<0,1))&&~ZLVs,errordlg('INI1: Minima of Vs are not increasing downwards');end
if ~isempty(find(diff(var(:,8))<0,1))&&~ZLVs,errordlg('INI1: Maxima of Vs are not increasing downwards');end    
if ~isempty(find(diff(var(:,4))<0,1))&&~ZLVp,errordlg('INI1: Minima of Vp are not increasing downwards');end
if ~isempty(find(diff(var(:,5))<0,1))&&~ZLVp,errordlg('INI1: Maxima of Vp are not increasing downwards');end
if ISDEPTH
    if ~isempty(find(diff(var(1:end-1,1))<0,1)),errordlg('INI1: Minima of depth are not increasing downwards');end
    if ~isempty(find(diff(var(1:end-1,2))<0,1)),var(1:end-1,2),errordlg('INI1: Maxima of depth are not increasing downwards');end
    %                 (max i)   >  min (i+1)                             (min i)   <  max (i-1)
    %LinkedInter=[(var(1:end-2,2)>var(2:end-1,1));false] | [false;(var(2:end-1,1)<var(1:end-2,2))>0];            
    % Equivalentemente:
    LinkedInter=[var(1:end-2,2)>var(2:end-1,1);false];
    LinkedInter=LinkedInter|[false;LinkedInter(1:end-1)];          
end

if ZLVp 
    if find(max(var(:,7),var(:,4).*sqrt((1-(2*var(:,14)))./(2*(1-var(:,14)))))-min(var(:,8),var(:,5).*sqrt((1-(2*var(:,13)))./(2*(1-var(:,13)))))>epsil)
        errordlg('INI1: Problems with Poisson ratio')
        %max(var(:,7),var(:,4).*sqrt((1-(2*var(:,14)))./(2*(1-var(:,14)))))>min(var(:,8),var(:,5).*sqrt((1-(2*var(:,13)))./(2*(1-var(:,13)))))
    end
else
    if find(max(var(:,4),var(:,7).*sqrt((2*(1-var(:,13)))./(1-(2*var(:,13)))))-min(var(:,5),var(:,8).*sqrt((2*(1-var(:,14)))./(1-(2*var(:,14)))))>epsil)
        errordlg('INI1: Problems with Poisson ratio')
        %max(var(:,4),var(:,7).*sqrt((2*(1-var(:,13)))./(1-(2*var(:,13)))))-min(var(:,5),var(:,8).*sqrt((2*(1-var(:,14)))./(1-(2*var(:,14)))))
    end
end

randm=rand(NC,4);% collection of random numbers

%% Initializing random parameters
[d11,alfa11,bta11,ro11]=deal(zeros(1,NC));
% Take random values of the "independent" velocity, fulfilling the constraints
% If ZLVs=0, the independent velocity is Vs
% If ZLVp=0, the independent velocity is Vp
% If ZLVs=1 and ZLVp=1, the independent velocity is Vs
if ZLVs && ZLVp, % Low velocity zones are allowed in both Vs and Vp.
    if na_out>1 && na_in<6, pol=[];limseg=[];end
    if na_in<7  || ~HHS, % really independent variables. HARD HALFSPACE DEFAULT : FALSE
   %if na_in>=7 && ~HHS, % really independent variables. HARD HALFSPACE DEFAULT : TRUE
        for kl=1:NC 
            if var(kl,9)
                bta11(kl)=var(kl,7)+(var(kl,9)*randm(kl,3));
            else
                bta11(kl)=var(kl,7);
            end
        end
    else % The halfspace has to have higher Vs 
        if (na_in<11 || isempty(polHHS))
            if na_out>3
                [bta11(NC),polHHS,facHHS,ProbLimSupSegHHS,LimVelSegHHS]=P_inv_HS(randm(NC,3),var(:,7:8));
            else
                bta11(NC)=P_inv_HS(randm(NC,3),var(:,7:8));
            end
        else
            bta11(NC)=P_inv_HS(randm(NC,3),var(:,7:8),polHHS,facHHS,ProbLimSupSegHHS,LimVelSegHHS);
        end
        dif=min(var(:,8),bta11(NC))-var(:,7);
        for kl=1:NC-1 
            if dif(kl)>0
                bta11(kl)=var(kl,7)+(dif(kl)*randm(kl,3));
            elseif dif(kl)==0
                bta11(kl)=var(kl,7);
            else
                errordlg(['INI1: Error while setting Hard Halfspace condition in layer ',num2str(kl)]);bta11(kl)=var(kl,7);pause;
            end
        end
    end
elseif ~ZLVs % Low velocity zones in Vs are not allowed.
    if na_in<5 || isempty(pol) % pol,limseg are not available
        if na_out>1
            [bta11,pol,limseg]=P_inv(randm(:,3),var(:,7:8));
        else
            bta11=P_inv(randm(:,3),var(:,7:8));
        end
    else
        bta11=P_inv(randm(:,3),var(:,7:8),pol,limseg);        
    end
elseif ~ZLVp %  Low velocity zones in Vp are not allowed.
    if na_in<5 || isempty(pol)% pol,limseg are not available
        if na_out>1
            [alfa11,pol,limseg]=P_inv(randm(:,2),var(:,4:5));
        else
            alfa11=P_inv(randm(:,2),var(:,4:5));
        end
    else
        alfa11=P_inv(randm(:,2),var(:,4:5),pol,limseg);        
    end
end

% Set THICKNESS of DEPTH
if ISDEPTH
    % Interface depths    
    d11(~LinkedInter)=var(~LinkedInter,1)+(var(~LinkedInter,3).*randm(~LinkedInter,1));
    if any(LinkedInter)
        d11(LinkedInter)=P_inv(randm(LinkedInter,1),var(LinkedInter,1:2));
    end
else
    % Random thickness and density with uniform probability
    for kl=1:NC-1
        if var(kl,3)
            d11(kl)=var(kl,1)+(var(kl,3)*randm(kl,1));
        else
            d11(kl)=var(kl,1);
        end
    end
end
%d11(NC)=0;

% Take random values of the "dependent" velocity. A uniform distribution is
% used within the limits determined by the Poisson's ratio and the
% dependent velocity previously obtained.
% Here we also random density with uniform probability within the whole 
% ranges provided by the user 
for kl=1:NC 
   
    % SET VS
    if ~ZLVp
        if var(kl,14)<0.4999
            aux1=max(var(kl,7),alfa11(kl).*sqrt((1-(2*var(kl,14)))./(2*(1-var(kl,14)))));% Minimum vS
        else
            aux1=var(kl,7);
        end
        aux2=min(var(kl,8),alfa11(kl).*sqrt((1-(2*var(kl,13)))./(2*(1-var(kl,13)))))-aux1;% Amplitude of  vS range
        if aux2<-epsil,
            errordlg(['INI1: INCOMPATIBLE RANGES OF VP, VS AND POISSON RATIO IN LAYER ',num2str(kl)])
            return;
        elseif aux2>0
            bta11(kl)=aux1+(aux2*randm(kl,3));
        else
            bta11(kl)=aux1;
        end
    end    
    
    % SET VP
    if ZLVp
        aux1=max(var(kl,4),bta11(kl)*sqrt((2*(1-var(kl,13)))/(1-(2*var(kl,13)))));% Minimum vP
        aux2=min(var(kl,5),bta11(kl)*sqrt((2*(1-var(kl,14)))/(1-(2*var(kl,14)))))-aux1;% Amplitude of  vP range
        if aux2<-epsil,
            errordlg(['INI1: INCOMPATIBLE RANGES OF VP, VS AND POISSON RATIO IN LAYER ',num2str(kl)])
            return;
        elseif aux2>0
            alfa11(kl)=aux1+(aux2*randm(kl,2));
        else
            alfa11(kl)=aux1;
        end
    end
    
    % DENSITY
    if var(kl,12),
        ro11(kl)=var(kl,10)+(var(kl,12)*(randm(kl,4)));
    else
        ro11(kl)=var(kl,10);
    end
end

%% Final Model Matrix (uses thicknesses always)
if ISDEPTH
    a=[[diff([0;d11(1:end-1)'],1,1);0] alfa11' bta11' ro11'];
else
    a=([ d11' alfa11' bta11' ro11' ]);
end

if find(a(:,3)./a(:,2)<sqrt((1-(2*var(:,14)))./(2*(1-var(:,14))))& a(:,2)./a(:,3)<sqrt((2*(1-var(:,14)))./(1-(2*var(:,14))))),fprintf(1,'INI1 :se supera el poi maximo');end
if find(a(:,3)./a(:,2)>sqrt((1-(2*var(:,13)))./(2*(1-var(:,13))))& a(:,2)./a(:,3)>sqrt((2*(1-var(:,13)))./(1-(2*var(:,13))))),fprintf(1,'INI1 :no se llega al poi minimo');end
if find(a(:,3)<var(:,7)),errordlg('INI1: Vs does not reach the minimum');end
if find(a(:,3)>var(:,8)),errordlg('INI1: Vs is larger than the maximum');end
if find(a(:,2)<var(:,4)),errordlg('INI1: Vp does not reach the minimum');end
if find(a(:,2)>var(:,5)),errordlg('INI1: Vp is larger than maximum');end  

return

function [x,pol,limseg]=P_inv(P,val,varargin)
% Generation of a random increasing-downwards velocity model 
% Garcia-Jerez et al. (2016, Computers & Geosciences, Appendix 1)
%
% Mandatory inputs:
% P = random values between 0 and 1 obtained from a uniform distribution.
%     LENGTH(P) is the number of layers.
% val = matrix of parameters ranges.
%       val(index,1) is the minimum velocity for the index-th layer (as an independent variable)
%       val(index,2) is the maximum velocity for the index-th layer (as an independent variable)
%
% Output:
% x = output velocity model. LENGTH(x) is the number of layers.
%     x(index) solves P(index) = Pacum(x(index)) (Press et al. 2007, section 7.3.2)
%     where Pacum is increasing function with values from 0 and 1 and 
%     is defined in the range allowed for the velocity of the index-th layer, 
%     taking into account the constraint of increasing velocity and the ramdom
%     velocities already generated for other layers (if any).
% pol = polinomials describing probability densities for each layer and range of velocity values 
% limseg = limits of probability intervals at which a specific polinomial applies
%
% Other:
% NC = number of layers
% vel = minimum velocity for the current layer (compatible with current constraints)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if nargin>2,
    pol=varargin{1};
    limseg=varargin{2};
end
NC=size(val,1);
x=zeros(1,NC);

% We first compute matrices pol and limseg, which depend only on val 
% (ranges of the parameters). They can be saved and used as input in sucesive 
% calls to this function. This makes sense for random exploration of fixed
% parameter intervals
if nargin<4
    pol=zeros(NC+1,NC,NC+1);
    limseg=zeros(NC,NC);% (layer_index,range_index) Probability of having a velocity higher than the minimun of the zone
    %polnew=zeros(1,NC+1);
    pol(NC+1,NC,NC+1)=1;% Dimensions mean: layer, range, monomial. The polinomial in the velocity of the layer index are stored at pol(index+1,:,:) 
    for index=NC:-1:1,  % We are integrating the probability density of the index-th layer velocity. 
        acum_seg=0;  
        for indx=NC:-1:index, % ranges with a specific polinomial for the current layer. We start from the highest velocity range
            % Take maximum and minimum velocity for the indx-th segment (xmaxseg,xminseg)
            if indx==NC, % Last interval of the current layer
                xmaxseg=val(index,2);% Max. velocity of the current (index) layer
                xminseg=val(indx,1);
            else
                xminseg=val(indx,1);
                xmaxseg=min(val(indx+1,1),val(index,2));            
            end
            if xminseg<=xmaxseg            
                if diff(val(index,:))<1,% especial if the interval is too narrow
                    if indx<NC && limseg(index,indx+1)>0
                        pol(index,indx,NC+1)=pol(index,indx+1,NC+1);                    
                        limseg(index,indx)=limseg(index,indx+1);
                    else
                        limseg(index,indx)=polyval(squeeze(pol(index+1,indx,index+1:NC+1)),xmaxseg);
                        pol(index,indx,NC+1)=limseg(index,indx);                 
                    end
                    acum_seg=limseg(index,indx);
                else % Regular behavior: integrate
                    polnew=0;        
                    polnew(NC+1)=polyval([squeeze(pol(index+1,indx,index+1:NC+1))'./(NC-index+1:-1:1),0],xmaxseg);      
                    polnew(index:NC)=-squeeze(pol(index+1,indx,index+1:NC+1))'./(NC-index+1:-1:1);        
                    polnew(index:NC+1)=polnew(index:NC+1)/diff(val(index,:),1,2);
                    polnew(NC+1)=polnew(NC+1)+acum_seg;
                    pol(index,indx,index:NC+1)=polnew(index:NC+1);% save the integrated polinomial
                    acum_seg=polyval(polnew(index:NC+1),xminseg);
                    limseg(index,indx)=acum_seg;                    
                end  
            else
                limseg(index,indx)=0;
                pol(index,indx,index:NC+1)=0;
            end
        end
        if index>1,
            pol(index,index-1,NC+1)=limseg(index,index);
        end
    end
end

% Now, we build up the random model from the upper to the bottom layer
vel=val(1,1);% Minimum allowed velocity for the current layer. Consider the upper one now.   
for lyr=1:NC    
    
    vel=max(vel,val(lyr,1)); %Minimum velocity for the layer
    
    % Normalizing the "reverse" cumlative probability for the layer lyr. It has to be 1 for the minimum velocity vel
    isegmin=find(vel>=val(:,1),1,'last');
    limseg(lyr,:)=limseg(lyr,:)/polyval(squeeze(pol(lyr,isegmin,lyr:NC+1)),vel);
    pol(lyr,:,lyr:NC+1)=pol(lyr,:,lyr:NC+1)/polyval(squeeze(pol(lyr,isegmin,lyr:NC+1)),vel);
    
    % Find out which range contains the provided uniform pobability P.
    isegmax=find((1-P(lyr))<=limseg(lyr,:),1,'last');
        
    % Refer the cumulative probability to vel instead of to the maximum
    % of the layer: P(x<vmax) -> P(vel<x)
    FinalPol=-squeeze(pol(lyr,isegmax,lyr:NC+1));
    FinalPol(end)=FinalPol(end)+1;
    
    % Find out which velocity limits correspond to the isegmax range.
    if isegmax==NC,
        xmaxseg=val(lyr,2);
        xminseg=val(isegmax,1);
    else
        xminseg=val(isegmax,1);
        xmaxseg=min(val(isegmax+1,1),val(lyr,2));
    end
    
%     figure;
%     axes;hold on;
%     for indx=isegmin:isegmax
%         if indx==NC,
%             xmaxsegplot=val(lyr,2);
%             xminsegplot=val(indx,1);
%         else
%             xminsegplot=val(indx,1);
%             xmaxsegplot=min(val(indx+1,1),val(lyr,2));
%         end        
%         xx=xminsegplot+[0:0.01:1]*(xmaxsegplot-xminsegplot);
%         plot(xx,polyval(squeeze(pol(lyr,indx,:)),xx)); 
%     end
    
    % Solve P = Pacum(x) = Prob[vel< layer_velocity <x]
    if xminseg==xmaxseg
        x(lyr)=xminseg;        
    elseif FinalPol==0,
        x(lyr)=xmaxseg;
    else       
        x(lyr)=fzero(@(y)polyval([FinalPol(1:end-1)',FinalPol(end)-P(lyr)],y),[xminseg,xmaxseg]);
    end
    vel=x(lyr);% prepared for the next (deeper) layer
    
end
return

function [x,pol,fac,ProbLimSupSeg,LimVelSeg]=P_inv_HS(P,val,varargin) 

% Generation of a random models with the halfspace has the higher velocity 
% Garcia-Jerez et al. (2016, Computers & Geosciences, Appendix 1)
% P = random value between 0 and 1 obtained from a uniform distribution.
% x = output velocity for the halfspace.
if nargin>2,
    pol=varargin{1};
    fac=varargin{2};
    ProbLimSupSeg=varargin{3};
    LimVelSeg=varargin{4};
end
NC=size(val,1);
if val(NC,1)==val(NC,2),
    x=val(NC,1);
    pol=[];fac=[];ProbLimSupSeg=[];LimVelSeg=[];
    return;
end

% Normalization of val to improve stability
minval=max(val(:));maxdif=max(diff(val,1,2));val=(val-minval)/maxdif;
%minval=0;maxdif=1;%skip normalization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate reusable matrices pol,fac,ProbLimSupSeg,LimVelSeg
if nargin<6
    LimVelSeg=unique(val(val(:)>=val(NC,1)&val(:)<=val(NC,2)),'sorted');
    pol=zeros(length(LimVelSeg)-1,NC+1);%Dimmensions correspond to segment(range) and monomial, respectively
    pol(:,NC+1)=1;
    fac=ones(1,length(LimVelSeg)-1);
    for index=NC-1:-1:1
         % pol(indx,index+1:NC+1) is a polynomial giving the product of independent probabilities of having,
         %                        from the last layer ,(NC-1), to the layer index a velocity smaller than x,
         %                        belonging to the segment indx
         % imax is the maximum segment interval id to be considered for the current layer (index)
         % imin is the minimum segment interval id to be considered for the current layer (index)
        if val(index,2)<LimVelSeg(1)
            continue;% Probability is 1 for any x in the interval; skip trivial multiplication by 1
        elseif val(index,2)>LimVelSeg(end)
            imax=length(LimVelSeg)-1;            
        else
            imax=find(LimVelSeg==val(index,2),1)-1;% Probability is one in segments > imax 
        end        
        if val(index,1)<LimVelSeg(1)
            imin=1;
        elseif val(index,1)>=LimVelSeg(end)            
            fprintf(1,'panic\n');pause
        else
            imin=find(LimVelSeg==val(index,1),1);% Probability is zero in segments < imin 
        end
        pol(1:imin-1,:)=0;       
        for indx=imin:imax,
            % Compute the product of independent probabilities of having, from the last layer (NC-1) to the
            % current layer index, a velocity smaller than x, belonging to the segment indx. It is polyval(pol(indx,index+1:NC+1),x)
            % (x - val(index,1))/diff(val(index,:) is the probability for velocity of layer index < x,
            % where x is in any segment from imin to imax
            pol(indx,index+1:NC+1)=conv(pol(indx,index+2:NC+1),[1 -val(index,1)]);
            fac(indx)=fac(indx)/diff(val(index,:));% for stability, this cumulative factor missing in  pol(indx,:), for segment indx
        end
    end
    
%     figure;
%     subplot(2,1,1);hold on;
%     for indx=1:length(LimVelSeg)-1
%         xx=LimVelSeg(indx)+[0:0.01:1]*(LimVelSeg(indx+1)-LimVelSeg(indx));
%         plot([0 0]+maxdif*LimVelSeg(indx)+minval,[0 1])
%         plot(maxdif*xx+minval,fac(indx)*polyval(squeeze(pol(indx,:)),xx));
%     end

    
    % Integrate the probability density of the halfspace
    ProbLimSupSeg=zeros(1,length(LimVelSeg)-1);
    for indx=1:length(ProbLimSupSeg)
        pol(indx,:)=[pol(indx,2:NC+1)./(NC:-1:1) 0];% Integrate. pol(indx,NC+1)=0.      
        pol(indx,NC+1)=-polyval(pol(indx,:),LimVelSeg(indx));%   We try to improve stability shifting to zero the
                                                                  %   polynomial at the lower boundary of any interval
                                                                  %   This behavior is kept
        ProbLimSupSeg(indx)= polyval(pol(indx,:),LimVelSeg(indx+1));% This is the shift in the integrated probability
                                                                    % between LimVelSeg(indx) and  LimVelSeg(indx+1)
                                                                    % fac(indx) is missing
    end
    ProbLimSupSeg=cumsum(fac.*ProbLimSupSeg);% Sets ProbLimSupSeg to make the cumulative probability continuous.
    % The polynomials (via fac) and ProbLimSupSeg are going to be rescaled to set total probability to one,
    % according to the normalization in Eq A5 (Garcia-Jerez et al 2016 CAGEO)
    fac=fac/ProbLimSupSeg(end);    
    ProbLimSupSeg=ProbLimSupSeg/ProbLimSupSeg(end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% subplot(2,1,2);hold on;
% for indx=1:length(LimVelSeg)-1
%     xx=LimVelSeg(indx)+[0:0.01:1]*(LimVelSeg(indx+1)-LimVelSeg(indx));
%     if indx==1,offset=0;else offset=ProbLimSupSeg(indx-1);end
%     plot([0 0]+maxdif*LimVelSeg(indx)+minval,[0 1])   
%     plot(maxdif*xx+minval,offset+fac(indx)*polyval(pol(indx,:),xx));
% end
% pause

% Now, we compute the random halfspace velocity

% Identify segment and velocity limits
iseg=find(P<=ProbLimSupSeg,1,'first');
FinalPol=fac(iseg)*pol(iseg,:);
xminseg=LimVelSeg(iseg);
xmaxseg=LimVelSeg(iseg+1);

% Solve P = Pacum(x) = Prob[min. vel. < halfspace_velocity <x]
if xminseg==xmaxseg
    x=xminseg;
else
    if iseg==1
        x=fzero(@(y)polyval([FinalPol(1:end-1),FinalPol(end)-P],y),[xminseg,xmaxseg]);
    else
        x=fzero(@(y)polyval([FinalPol(1:end-1),FinalPol(end)-P+ProbLimSupSeg(iseg-1)],y),[xminseg,xmaxseg]);
    end
end
x=x*maxdif+minval;% Undo normalization
return