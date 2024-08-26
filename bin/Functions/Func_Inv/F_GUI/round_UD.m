function [rx] = round_UD(x,pos,dir)
% Rounds up or down the pos-th decimal position of x
% dir = 'u' for for rounded x >= x
% dir = 'd' for for rounded x <= x
% pos = 0 for unity, 1 for tenths, 2 for hundredths 
% x can be a vector o a matrix
% pos is scalar
rx=round(x,pos);
sg=sign(x-rx);
switch dir
case {'u','U'}
    corregir=(sg==1);
    rx(corregir)=round(rx(corregir)+10^-pos,pos);
case {'d','D'}
    corregir=(sg==-1);
    rx(corregir)=round(rx(corregir)-10^-pos,pos);
end        
end

