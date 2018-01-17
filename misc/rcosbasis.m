% This code was originally written by Pavan Ramkumar, a former postdoc in
% the Kording lab at Northwestern University. I received his written
% permission to use and publish this script on github. 

% -------------------------------------------
% Make a raised cosine basis
% -------------------------------------------
% function [t, basis] = rcosbasis(t, c, w);
%
% t is a vector of time points in samples
% c is a vector of raised cosine centers
% w is a vector of raised cosine widths
%
% (c) Pavan Ramkumar 2013
%--------------------------------------------

function [t, basis] = rcosbasis(t, c, w)

basis = zeros(length(t)+1, length(c));

% make t into a column vector
if size(t,2) > size(t,1) t = t'; end
    
for i=1:length(c)
  basis(2:end,i) = (cos(max(-pi,min(pi,(t-c(i))*pi/w(i)/2)))+1)/2;
end

basis(end,:) = [];
%basis = orth(basis);