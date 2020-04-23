function varargout = db2lin(varargin)
% DESCRIPTION y = db2lin(x)
%  Converts dB to linear.
%  Takes any number of arguments -Inf will map to 0.
%  Often used to go from dBm to mW.
% INPUT
%  Any number of any real matrix.
% OUTPUT
%  y -- The input arguments converted.  
% TRY
%  db2lin(3), db2lin([0 -inf]), db2lin([0 0], -3+zeros(2,2,2)) 
% SEE ALSO
%  lin2db

% by Magnus Almgren 990301

for i = 1:nargin
 varargout{i} = 10 .^ (varargin{i}/10);
end
