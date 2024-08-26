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

function VisualFit_v8(MDLid,handles)

global ground_par keep_poi keep_hv LockAxis change1thickness dataFHV pol_backup
%
% INPUTS
% ground_par= 2 for Vp
%             3 for Vs
%             4 for density
% MDLid = handle to aux model line
% keep_poi = true /false
% keep_hv = true / false
% LockAxis = 1 locks X (ground model parameter); 2 locks Y (thickness)
% change1thickness=0 moves interface; =1 changes thickness of the layer above the vertex
%
% OUTPUTS
% dataFHV
% pol_backup
% Update of handles.uitable1,'Data'
%
%

MIN_PAR=[1 2  1   100];% Minima for Thickness,Vp,Vs,density
MAX_PAR=[inf 90000 30000 10000];% Maxima for Thickness,Vp,Vs,density

%%%%%%%%%%%%%%% Obtain X - Y %%%%%%%%%%%%%%%%%%%
[B,A,R,E]=MODELO(dataFHV(:,1:4),0);
E(end)=sum(dataFHV(:,1))*1.25;
ProfilePolys=[E;A;B;R].';

pol_backup=[ProfilePolys(:,ground_par) ProfilePolys(:,1) ];

%%%%%%%%%%%% Make de poligon and save handles to base workspace
[MDLid.XData,MDLid.YData,MDLid.Parent.YDir]=deal(pol_backup(:,1),pol_backup(:,2),'reverse');
h = impoly(MDLid.Parent,pol_backup,'Closed',false,'PositionConstraintFcn',@modify);
g=get(h);
LoadV('h',h,'g',g);
    function pol=modify(pol)
            
        if keep_hv && (ground_par==2 || ground_par==3)
            % Aborts left-down and right-up motion of vertexes
            ifact=find(pol(3:2:end-1,2)~=pol_backup(3:2:end-1,2),1);
            if ~isempty(ifact),
                if sign(pol(ifact*2+1,2)-pol_backup(ifact*2+1,2))*sign(pol(ifact*2+1,1)-pol_backup(ifact*2+1,1))<0
                    pol=pol_backup;
                end                
            end
        end
        % Prohibe cambios en el eje bloqueado
        if LockAxis==1 && ~keep_hv             % keep_hv es incompatible con LockAxis==1
            pol(:,1)=pol_backup(:,1);
        end
        if LockAxis==2 || keep_hv              % Por simplicidad, keep_hv aplica un bloqueo temporal al eje Y
            pol(:,2)=pol_backup(:,2);
        else
            % Conserva el arranque del perfil en la superficie libre
            if pol(1,2)~=0,pol(1,2)=0;end
            if ~isempty(find(pol(3:2:end-1,2)<0,1)),pol(:,2)=pol_backup(:,2);end% Restauración total de las posiciones de las interfaces;        
            if change1thickness
                % Mover solidariamente las interfaces, conservando espesores
                pol(1:2:end-1,2)=pol_backup(1:2:end-1,2)+cumsum(pol(1:2:end-1,2)-pol_backup(1:2:end-1,2));
            end
            % Verifica que los espesores/profundidades son consistentes y conserva la parte horizontal del escalón, horizontal
            kk=find(pol(3:2:end-1,2)~=pol(2:2:end-1,2));% vértices móviles que se han movido
            aux=MIN_PAR(1)-(pol(2*kk+1,2)-pol(2*kk-1,2));
            ibad=find(aux>=0,1);% Al menos el kk(ibad)-ésimo vértice movil se ha movido demasiado arriba, hay que bajarlo aux(ibad)
            if ~isempty(ibad)
                if change1thickness
                    % Llevo el (primer) vertice problemático a la profundidad mínima,
                    % moviendo solidariamente hacia abajo todo el modelo tentativo
                    %pol(2*kk(ibad)+1:2:end,2)=pol(2*kk(ibad)+1:2:end,2)+aux(ibad);
                    pol(2*kk(ibad)+1:2:end,2)=pol(2*kk(ibad)-1,2)+MIN_PAR(1)+pol_backup(2*kk(ibad)+1:2:end,2)-pol_backup(2*kk(ibad)+1,2);

                else
                    % Llevo el (primer) vertice problemático a la profundidad mínima,
                    pol(2*kk(ibad)+1,2)=pol(2*kk(ibad)-1,2)+MIN_PAR(1);
                end
                %pol(2*kk(ibad(1))+1:2:end,2)=pol(2*kk(ibad(1)):2:end-1,2);% Prohibo el movimiendo desde ese vértice hacia abajo, copiando de los fijos
            end
            %if ~isempty(ibad),fprintf(1,'arreglo del tentativo (1)\n');pol,end

            if ~change1thickness
                % En este modo, tiene sentido comprobar que el vértice movido no
                % está demasiado cerca del vértice movil inferior
                ibad=(pol(2*kk(2*kk+3<=size(pol,1))+3,2)-pol(2*kk(2*kk+3<=size(pol,1))+1,2))<MIN_PAR(1);
                pol(2*kk(ibad)+1,2)=pol(2*kk(ibad)+3,2)-MIN_PAR(1);% El vértice movil al menos MIN_PAR(1) mas superficial que la capa siguiente
                %if ~isempty(ibad),fprintf(1,'arreglo del tentativo (2)\n');pol,end
            end
            [pol(2*kk,2),pol(end,2)]=deal(pol(2*kk+1,2),pol(end-1,2)*1.25);% enderezamos el escalón
            if pol(end,2)==0,pol(end,2)=20;end;% 20m for representing a halfspace            
        end
        
        % Keep the minimum of the ground property
        kk=find(pol(1:2:end,1)<MIN_PAR(ground_par));
        pol(2*kk-1,1)=MIN_PAR(ground_par);
        
        % Keep the maximum of the ground property
        kk=find(pol(1:2:end,1)>MAX_PAR(ground_par));
        pol(2*kk-1,1)=MAX_PAR(ground_par);
        
        if keep_hv
            ifact=find(pol(1:2:end-1,1)~=pol_backup(1:2:end-1,1),2);
            if length(ifact)==2
                pol=pol_backup;
            elseif ~isempty(ifact),
                % compute the factor introduced for the ground property
                fact=(pol(2*ifact-1,1)/pol_backup(2*ifact-1,1));
                % Apply this factor to the rest of the odd vertexes
                pol(1:2:2*ifact-3,1)=pol_backup(1:2:2*ifact-3,1)*fact;
                pol(2*ifact+1:2:end,1)=pol(2*ifact+1:2:end,1)*fact;
                if ~isempty(find(pol(1:2:end-1,1)>MAX_PAR(ground_par),1))||...
                   ~isempty(find(pol(1:2:end-1,1)<MIN_PAR(ground_par),1))
                   % Abort if the property is out of range for any layer 
                   pol=pol_backup;
                elseif ground_par~=4
                    % Change thicknesses to preserve travel time
                    pol(:,2)=pol_backup(:,2)*fact;
                    % Abort if the thickness is out of range for any layer
                    if find(pol(2:2:end-2,2)-pol(1:2:end-3,2)<MIN_PAR(1),1)
                        pol=pol_backup;
                    end
                end
            end
        end
        
        % Conserva la parte vertical del escalón, vertical
        kk=find(pol(1:2:end,1)~=pol(2:2:end,1));
        pol(2*kk,1)=pol(2*kk-1,1);

        if (keep_poi||keep_hv)&&(ground_par==2||ground_par==3)
            % Hay que cambiar la otra velocidad para conservar el Poisson
            if ground_par==2
                TentativeVel=real(pol(1:2:end-1,1).*...
                    (sqrt((1-(2.*dataFHV(:,5)))./...
                    (2*(1-dataFHV(:,5))))));
                if find(TentativeVel>MAX_PAR(3)|TentativeVel<MIN_PAR(3),1),
                    pol=pol_backup;
                else
                    % Changing the changed velocity (current ground
                    % parameter) and the other velocity, stored in TentativeVel
                    [pol_backup,dataFHV(:,1),dataFHV(:,2)]=deal(pol,[diff(pol(1:2:end-1,2));0],pol(1:2:end-1,1));
                    dataFHV(:,3)=TentativeVel;
                end
            elseif ground_par==3
                TentativeVel=real(pol(1:2:end-1,1).*...
                    (sqrt((2*(1-dataFHV(:,5)))./...
                    (1-(2.*dataFHV(:,5))))));
                if find(TentativeVel>MAX_PAR(2)|TentativeVel<MIN_PAR(2),1),
                    pol=pol_backup;
                else
                    % Changing the changed velocity (current ground
                    % parameter) and the other velocity, stored in TentativeVel
                    [pol_backup,dataFHV(:,1),dataFHV(:,3)]=deal(pol,[diff(pol(1:2:end-1,2));0],pol(1:2:end-1,1));
                    dataFHV(:,2)=TentativeVel;
                end                
            end
        elseif ground_par==4
            % Changing density
            % Save model pol as pol_backup and replace thickness and densities in dataFHV
            [pol_backup,dataFHV(:,1),dataFHV(:,4)]=deal(pol,[diff(pol(1:2:end-1,2));0],pol(1:2:end,1));
        else
            % Changing one velocity and updating Poisson ratio
            % Save model pol as pol_backup and replace the velociy in dataFHV
            [pol_backup,dataFHV(:,1),dataFHV(:,ground_par)]=deal(pol,[diff(pol(1:2:end-1,2));0],pol(1:2:end-1,1));            
            % Save update Poisson ratio in dataFHV
            dataFHV(:,5)=(0.5-(dataFHV(:,3)./dataFHV(:,2)).^2)./...
                         (1-(dataFHV(:,3)./dataFHV(:,2)).^2);
        end
        
        [bmodel,limddFHV,ZoomOn]=GetV('bmodel','limddFHV','ZoomOn');

%         get(VScbarFHV.VerticalScrollBar,'Value')

%         get(VScbarFHV.VerticalScrollBar,'Value')
        pol(end,2)=max(limddFHV*1.25,pol(end,2));
        % Actualizando rango de los ejes dependiendo el modelo
        if ~ZoomOn
            if length(bmodel)==1;
                handles.uipanel14.Children.XLim=[min(dataFHV(:,ground_par))-min(dataFHV(:,ground_par))*0.05...
                    max(dataFHV(:,ground_par))+max(dataFHV(:,ground_par))*0.05  ];
            else
                handles.uipanel14.Children.XLim=[min(min(bmodel(:,ground_par))-min(bmodel(:,ground_par))*0.05 , min(dataFHV(:,ground_par))-min(dataFHV(:,ground_par))*0.05)...
                    max(max(bmodel(:,ground_par))+max(bmodel(:,ground_par))*0.05 , max(dataFHV(:,ground_par))+max(dataFHV(:,ground_par))*0.05)];
            end
            handles.uipanel14.Children.YLim=[0 max(limddFHV*1.25,pol(end,2)) ];
        end
        if length(bmodel)>=1;
        MDLid.YData(end)=max(limddFHV*1.25,pol(end,2));
        end  
%        set(handles.uitable1,'Data',num2cell(dataFHV));
% handles.uitable1.Enable='off';
        handles.uitable1.Data=num2cell(dataFHV);
%            [VScbarFHV]=GetV('VScbarFHV');
%          set(VScbarFHV.VerticalScrollBar,'Value',50);
%         set(VScbarFHV.VerticalScrollBar,'ValueIsAdjusting',10);
% handles.uitable1.Enable='on';
% set(VScbarFHV.VerticalScrollBar,'Enable',0)
    end % of modifica

% Actualizando rango de los ejes dependiendo el modelo
[bmodel,limddFHV]=GetV('bmodel','limddFHV');
if length(bmodel)==1;
    handles.uipanel14.Children.XLim=[min(dataFHV(:,ground_par))-min(dataFHV(:,ground_par))*0.05...
                                     max(dataFHV(:,ground_par))+max(dataFHV(:,ground_par))*0.05  ];
else
    handles.uipanel14.Children.XLim=[min(min(bmodel(:,ground_par))-min(bmodel(:,ground_par))*0.05 , min(dataFHV(:,ground_par))-min(dataFHV(:,ground_par))*0.05)...
                                     max(max(bmodel(:,ground_par))+max(bmodel(:,ground_par))*0.05 , max(dataFHV(:,ground_par))+max(dataFHV(:,ground_par))*0.05)];
end
handles.uipanel14.Children.XLim=[min(dataFHV(:,ground_par))-min(dataFHV(:,ground_par))*0.05...
                                 max(dataFHV(:,ground_par))+max(dataFHV(:,ground_par))*0.05  ];
handles.uipanel14.Children.YLim=[0 E(end) ];
handles.uitable1.Data=num2cell(dataFHV);


end % of VisualFit_v8
            
