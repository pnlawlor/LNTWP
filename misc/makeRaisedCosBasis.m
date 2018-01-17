function [iht, ihbas, ihbasis] = makeRaisedCosBasis(nb, dt, endpoints, ...
						  b, zflag);
                      
% This code is an earlier version of, or is based upon code produced by the
% Pillow laboratory, made available under the MIT license. I received it
% from a colleague, Ian Stevenson, and have tested it and found it
% appropriate for my own use. 
% The most updated version of this script can be found here: 
% https://github.com/pillowlab/raisedCosineBasis. 
% It was originally implemented in Pillow et al, Nature 2008. 
                      
% [iht, ihbas, ihbasis] = makeRaisedCosBasis(nb, dt, pks, b, zflag);
%
% Make nonlinearly stretched basis consisting of raised cosines
% Inputs:  nb = # of basis vectors
%          dt = time bin separation for representing basis
%          endpoints = 2-vector containg [1st_peak  last_peak], the peak 
%                  (i.e. center) of the last raised cosine basis vectors
%          b = offset for nonlinear stretching of x axis:  y = log(x+b) 
%              (larger b -> more nearly linear stretching)
%          zflag = flag for making (if = 1) finest-timescale basis
%                  vector constant below its peak
%
%  Outputs:  iht = time lattice on which basis is defined
%            ihbas = orthogonalized basis
%            ihbasis = basis itself
%
%  Example call
%  [iht, ihbas, ihbasis] = makeRaisedCosBasis(10, .01, [0 10], .1);

if nargin < 3
    if nargin == 2
	zflag = dt;
    else
	zflag = 0;
    end
    ihprs = nb;
    nb = ihprs.nh;
    dt = ihprs.hdt;
    if isfield(ihprs, 'endpoints')
	endpoints = ihprs.endpoints;
    elseif isfield(ihprs, 'hspan')
	endpoints = ihprs.hspan;
    else
	error('missing field');
    end
    b = ihprs.b;

elseif (nargin == 4)
    zflag = 0;
end

% nonlinearity for stretching x axis (and its inverse)
nlin = inline('log(x+1e-20)');
invnl = inline('exp(x)-1e-20');

if b <= 0
    error('b must be greater than 0');
end

yrnge = nlin(endpoints+b);  
db = diff(yrnge)/(nb-1);      % spacing between raised cosine peaks
ctrs = yrnge(1):db:yrnge(2);  % centers for basis vectors
mxt = invnl(yrnge(2)+2*db)-b; % maximum time bin
iht = [0:dt:mxt]';
nt = length(iht);        % number of points in iht
ff = inline('(cos(max(-pi,min(pi,(x-c)*pi/dc/2)))+1)/2', 'x', 'c', 'dc');  % raised cosine basis vector
ihbasis = ff(repmat(feval(nlin, iht+b), 1, nb), repmat(ctrs, nt, 1), db);

if zflag  % set first basis vector bins (before 1st peak) to 1
    ii = find(iht<=endpoints(1));
    ihbasis(ii,1) = 1;
end

ihbas = orth(ihbasis);  % use orthogonalized basis
