
% Compute and show mean model, std deviation, covariance matrix and 
% and correlation matrix
%
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
%

[DATAERROR,ISDEPTH]=GetV('DATAERROR','ISDEPTH'); % Load DATAERROR structure with a list of models and misfits
if isfield(DATAERROR, 'Mis')
    
    % Settings for representation
    linewidth=5;
    meanP='k';col_errbar_P=[0 0 0];
    meanS='k';col_errbar_S=[0 0 0];

    % Take number of models and layers
    ndatos=length(DATAERROR.Mis); % number of models
    ncapas=size(DATAERROR.mod(:,:,1),1);% number of layers, including halfspace
    ncapasX4=4*ncapas;    
    if ISDEPTH
        DATAERROR.mod(1:end-1,1,:)=cumsum(DATAERROR.mod(1:end-1,1,:),1);% No operation with thicknesses
        simbolos{1}='z';simbolos{2}='Vp';simbolos{3}='Vs';simbolos{4}='Dens';% for labels
    else
        % We use thickness for the correlation matrix. We will change to depths for the model later 
        simbolos{1}='h';simbolos{2}='Vp';simbolos{3}='Vs';simbolos{4}='Dens';% for labels
    end
    % Select model parameters for statistics (from left side to right side and from up to down)
    % E.G.: icoords=[1:4:4*(ncapas-1) 2:4:4*ncapas 3:4:4*ncapas ] for all thicknesses, P and S wave velocities
    icoords=[find(max(DATAERROR.mod(:,1,:),[],3)~=min(DATAERROR.mod(:,1,:),[],3))*4-3;...
             find(max(DATAERROR.mod(:,2,:),[],3)~=min(DATAERROR.mod(:,2,:),[],3))*4-2;...
             find(max(DATAERROR.mod(:,3,:),[],3)~=min(DATAERROR.mod(:,3,:),[],3))*4-1;...
             find(max(DATAERROR.mod(:,4,:),[],3)~=min(DATAERROR.mod(:,4,:),[],3))*4].'; % all changing parameters
    dim=length(icoords); 
    icoordsFil=ceil(icoords/4);
    icoordsCol=mod(icoords-1,4)+1;
    icoordstrans= icoordsFil+(icoordsCol-1)*ncapas;
    
    % Find best model and compute mean model
    [kk ,pos]=min(DATAERROR.Mis);
    best_model=DATAERROR.mod(:,:,pos).';
    modelo_medio_completo=reshape(DATAERROR.mod(:,:,pos).',1,[]);% Initializing.Row vector
    modelo_promedio=zeros(1,dim);% Initializing mean model.Row vector
    for index=0:ndatos-1,
        modelo_promedio=modelo_promedio+DATAERROR.mod(icoordstrans+index*ncapasX4);
    end
    modelo_promedio=modelo_promedio/ndatos;
    modelo_medio_completo(icoords)=modelo_promedio;

    if dim>0
        % Covariance matrix (matriz_cov)
        matriz_cov=zeros(dim);
        for index=0:ndatos-1,
            matriz_cov=matriz_cov+(DATAERROR.mod(icoordstrans+index*ncapasX4)-modelo_promedio).'*...
                                  (DATAERROR.mod(icoordstrans+index*ncapasX4)-modelo_promedio);
        end
        matriz_cov=matriz_cov/ndatos;

        % Compute and show correlation matrix
        matriz_corrs=matriz_cov(1:dim+1:end);% take diagonal from matriz_cov as a row [C11,C22,...]
        matriz_corrs=sqrt(matriz_corrs.' * matriz_corrs);% transform into matrix sqrt(Cii*Cjj)
        matriz_corrs=abs(matriz_cov)./matriz_corrs;% correlation matrix
        matriz_corrs_ampliada=zeros(size(matriz_corrs)+1);% a trick that avoid missing a column and a row in the representation
        matriz_corrs_ampliada(1:end-1,1:end-1) = matriz_corrs;
        matriz_corrs_ampliada(end,1:end-1)     = matriz_corrs_ampliada(end-1,1:end-1);
        matriz_corrs_ampliada(1:end-1,end)     = matriz_corrs_ampliada(1:end-1,end-1);
        matriz_corrs_ampliada(end,end)         = matriz_corrs_ampliada(end,end-1);
        figure;eje=axes;surf(matriz_corrs_ampliada);view([0 -90]);
        tickcell={dim};
        for index=1:dim,tickcell{index}=[simbolos{rem(icoords(index)-1,4)+1},num2str(ceil(icoords(index)/4))];end
        axis equal;axis tight;caxis([0 1]);colormap(1-gray);
        colorbar;
        set(eje,'xtickmode','manual','xticklabelmode','manual','xtick',(1:dim)+0.5,'xticklabel',tickcell);
        set(eje,'ytickmode','manual','yticklabelmode','manual','ytick',(1:dim)+0.5,'yticklabel',tickcell,'fontsize',14);
        set(eje,'ytickmode','manual','yticklabelmode','manual','ytick',(1:dim)+0.5,'yticklabel',tickcell,'fontsize',14);
        set(eje,'XAxisLocation','top');
    else
        errordlg(sprintf('Too few low-misfit models to estimate uncertainties. Try to:\n \n -Decrease perturbation range in Global optimizations\n -Check that experimental uncertainties are not too small\n -Increase Number of Iterations for Monte Carlo Sampling\n -Increase Last Màrkov Chain length for Simulated Annealing'));
    end

    if ~ISDEPTH
        % Replace computations using thickness with computations using depths
        DATAERROR.mod(1:end-1,1,:)=cumsum(DATAERROR.mod(1:end-1,1,:),1);% thicknesses to depths
        icoordsDEPTH=transpose(find(max(DATAERROR.mod(:,1,:),[],3)~=min(DATAERROR.mod(:,1,:),[],3))*4-3);
        dimDEPTH=length(icoordsDEPTH);
        dimTOTAL=dim+dimDEPTH;
        matriz_cov(dim+1:dimTOTAL,dim+1:dimTOTAL)=0;% Enlarge matriz_cov
        % cancel depths as valid variable:
        icoords(icoordsCol==1)=0;
        icoordsFil(icoordsCol==1)=0;
        icoordsCol(icoordsCol==1)=0;
        icoordstrans(icoordsCol==1)=0;
        % Compute the coordinate indexes for the moving interfaces:
        icoordsFilDEPTH=ceil(icoordsDEPTH/4);        
        icoordsColDEPTH=mod(icoordsDEPTH-1,4)+1;
        icoordstransDEPTH= icoordsFilDEPTH+(icoordsColDEPTH-1)*ncapas;
        % Append at the end of the complete list
        icoords=[icoords, icoordsDEPTH];
        icoordsFil=[icoordsFil,icoordsFilDEPTH];
        icoordsCol=[icoordsCol,icoordsColDEPTH];
        icoordstrans=[icoordstrans,icoordstransDEPTH];        
        % Find best model and compute mean model
%        [kk ,pos]=min(DATAERROR.Mis);
        modelo_medio_completo=reshape(DATAERROR.mod(:,:,pos).',1,[]);% Initializing with depths.Row vector
%        modelo_promedio=zeros(1,dim);% Initializing mean model.Row vector
        modelo_promedio=[modelo_promedio,zeros(1,dimDEPTH)];% Append depths of moving layers
        % Compute mean depths: 
        for index=0:ndatos-1,
            modelo_promedio(dim+1:dimTOTAL)=modelo_promedio(dim+1:dimTOTAL)+...
                                            DATAERROR.mod(icoordstransDEPTH+index*ncapasX4);
        end
        modelo_promedio(dim+1:dimTOTAL)=modelo_promedio(dim+1:dimTOTAL)/ndatos;
        % Obtain the full mean model in terms of depths Vp Vs density        
        modelo_medio_completo(icoords(icoords>0))=modelo_promedio(icoords>0);
        if dimTOTAL>0
            % Covariance matrix (matriz_cov)
            for index=0:ndatos-1,
    %             size(matriz_cov(dim+1:dimTOTAL,dim+1:dimTOTAL))
    % 
    %             size(DATAERROR.mod(icoordstransDEPTH+index*ncapasX4))
    %             size(modelo_promedio(dim+1:dimTOTAL))
    %             size((DATAERROR.mod(icoordstransDEPTH+index*ncapasX4)-modelo_promedio(dim+1:dimTOTAL)).'*...
    %                 (DATAERROR.mod(icoordstransDEPTH+index*ncapasX4)-modelo_promedio(dim+1:dimTOTAL)))


                matriz_cov(dim+1:dimTOTAL,dim+1:dimTOTAL)=matriz_cov(dim+1:dimTOTAL,dim+1:dimTOTAL)+...
                    (DATAERROR.mod(icoordstransDEPTH+index*ncapasX4)-modelo_promedio(dim+1:dimTOTAL)).'*...
                    (DATAERROR.mod(icoordstransDEPTH+index*ncapasX4)-modelo_promedio(dim+1:dimTOTAL));
            end
            matriz_cov(dim+1:dimTOTAL,dim+1:dimTOTAL)=matriz_cov(dim+1:dimTOTAL,dim+1:dimTOTAL)/ndatos;
        end
    end

    % Plot mean Vs profile and std deviation
    %d    = modelo_medio_completo(1:4:(ncapas-1)*4);
    %d(2:end)=diff(d(1:end-1));% Depth to thickness
    z=[0 modelo_medio_completo(1:4:(ncapas-1)*4)];% superficie libre incluida excluida
    bta  = modelo_medio_completo(3:4:ncapas*4);
    figure(5);subplot(1,3,2);hold on
    hg=evalin('base','hg');
    DATAERROR=evalin('base','DATAERROR');
    %y(2:2:2*length(d)+2)=[cumsum(d) DATAERROR.lim];
    %y(3:2:2*length(d)+1)=y(2:2:2*length(d));
    y(2:2:2*ncapas)=[z(2:end) DATAERROR.lim];
    y(3:2:2*ncapas-1)=y(2:2:2*ncapas-2);    
    x(2:2:2*ncapas)=bta;
    x(1:2:2*ncapas-1)=bta;
    dibujo=plot(x,y,'color',[.0 .0 .0],'LineWidth',linewidth);
    % Axes, limits and legend
    set(gca,'YDir','reverse');
    ylim([0 y(end)]);
    l4=legend([hg(2) dibujo],'Best model V_s','Mean model V_s','Location','southwest');
    set(l4,'FontSize',10);    
    % Error bars
    for index=1:length(icoords),
        if icoordsCol(index)==1
            oid=line([1 1]*mean(bta(icoordsFil(index):icoordsFil(index)+1)),[-1 1]*sqrt(matriz_cov(index,index))+z(icoordsFil(index)+1)); % plot standard deviation of layer depth
            set(oid,'linewidth',2,'color',col_errbar_S);% set color col_errbar_S
        elseif icoordsCol(index)==3
            if icoordsFil(index)<ncapas
                oid=line([-1 1]*sqrt(matriz_cov(index,index))+bta(icoordsFil(index)),[1 1]*(mean(z(icoordsFil(index):icoordsFil(index)+1))));% plot standard deviation of Vs
            else
                oid=line([-1 1]*sqrt(matriz_cov(index,index))+bta(icoordsFil(index)),[1 1]*(z(icoordsFil(index))+(max(y)-z(icoordsFil(index)))/2));% Same for halfspace            
            end
            set(oid,'linewidth',2,'color',col_errbar_S);% set color col_errbar_S
        end
    end

    % Plot mean Vp profile and std deviation
    alfa = modelo_medio_completo(2:4:ncapas*4);
    figure(5);subplot(1,3,1);hold on;
    x(2:2:2*ncapas)=alfa;
    x(1:2:2*ncapas-1)=alfa;
    dibujo=plot(x,y,'color',[.0 .0 .0],'LineWidth',linewidth);    
    % Axes, limits and legend
    set(gca,'YDir','reverse');
    ylim([0 y(end)]);
    l4=legend([hg(1) dibujo],'Best model V_p','Mean model V_p','Location','southwest');
    set(l4,'FontSize',10)
    % Error bars
    for index=1:length(icoords),
        if icoordsCol(index)==1
            oid=line([1 1]*mean(alfa(icoordsFil(index):icoordsFil(index)+1)),[-1 1]*sqrt(matriz_cov(index,index))+z(icoordsFil(index)+1)); % plot standard deviation of layer interface
            set(oid,'linewidth',2,'color',col_errbar_P);% set color col_errbar_P
        elseif icoordsCol(index)==2
            if icoordsFil(index)<ncapas
                oid=line([-1 1]*sqrt(matriz_cov(index,index))+alfa(icoordsFil(index)),[1 1]*(mean(z(icoordsFil(index):icoordsFil(index)+1)))); % plot standard deviation of layer Vs                
            else
                oid=line([-1 1]*sqrt(matriz_cov(index,index))+alfa(icoordsFil(index)),[1 1]*(z(icoordsFil(index))+(max(y)-z(icoordsFil(index)))/2));% Same for halfspace            
            end
            set(oid,'linewidth',2,'color',col_errbar_P);% set color col_errbar_P
        end
    end

    % Plot mean Rho profile and std deviation
    rho = modelo_medio_completo(4:4:ncapas*4);
    figure(5);subplot(1,3,3);hold on;
    x(2:2:2*ncapas)=rho;
    x(1:2:2*ncapas-1)=rho;
    dibujo=plot(x,y,'color',[.0 .0 .0],'LineWidth',linewidth);    
    % Axes, limits and legend
    set(gca,'YDir','reverse');
    ylim([0 y(end)]);
    l4=legend([hg(3) dibujo],'Best model \rho','Mean model \rho','Location','southwest');% 'best'
    set(l4,'FontSize',10)
    % Error bars
    for index=1:length(icoords),
        if icoordsCol(index)==1  
            oid=line([1 1]*mean(rho(icoordsFil(index):icoordsFil(index)+1)),[-1 1]*sqrt(matriz_cov(index,index))+z(icoordsFil(index)+1)); % plot standard deviation of layer interface
            set(oid,'linewidth',2,'color',col_errbar_P);% set color col_errbar_P
        elseif icoordsCol(index)==4
            if icoordsFil(index)<ncapas
                oid=line([-1 1]*sqrt(matriz_cov(index,index))+rho(icoordsFil(index)),[1 1]*(mean(z(icoordsFil(index):icoordsFil(index)+1)))); % plot standard deviation of layer Vs
            else
                oid=line([-1 1]*sqrt(matriz_cov(index,index))+rho(icoordsFil(index)),[1 1]*(z(icoordsFil(index))+(max(y)-z(icoordsFil(index)))/2));% Same for halfspace            
            end
            set(oid,'linewidth',2,'color',col_errbar_P);% set color col_errbar_P
        end
    end
    
    %figure(5);
    %subplot(1,3,3);
    %ylim([0 DATAERROR.lim]);%set depth limit for density subplot
end


