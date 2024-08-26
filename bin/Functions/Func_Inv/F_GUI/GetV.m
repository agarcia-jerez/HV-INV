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

%% Function for get variables from workspace
function varargout=GetV(varargin)
varargout=cell(nargin,1);
for i=1:nargin
    varargout{i,1}=evalin('base',varargin{i});
end