function y = medianoutnan(x,dim)
%MEDIAN Median value.
%   For vectors, MEDIANOUTNAN(X) is the median value of the finite 
%   elements in X. For matrices, MEDIANOUTNAN(X) is a row vector 
%   containing the median value of each column.  For N-D arrays, 
%   MEDIANOUTNAN(X) is the median value of the elements along the 
%   first non-singleton dimension of X.
%
%   MEDIANOUTNAN(X,DIM) takes the median along the dimension DIM of X.
%
%   Example: If X = [0 1 2   NaN
%                    3 4 NaN NaN]
%
%   then medianoutnan(X,1) is [1.5 2.5 2 NaN] 
%   and medianoutnan(X,2) is [1
%                             3.5]
%
%   See also MEANOUTNAN, STDOUTNAN, MIN, MAX, COV.

%   P. Lherminier, 21/03/2002. From modifications of median.m

if nargin==1, 
  dim = min(find(size(x)~=1)); 
  if isempty(dim), dim = 1; end
end
if isempty(x), y = []; return, end

siz = [size(x) ones(1,dim-ndims(x))];
n = size(x,dim);

% Permute and reshape so that DIM becomes the row dimension of a 2-D array
perm = [dim:max(length(size(x)),dim) 1:dim-1];
x = reshape(permute(x,perm),n,prod(siz)/n);

% Sort along first dimension
x = sort(x,1);
sizx=size(x);

[nnan,ncol]=find(diff(isnan(x)));
nn=nan*ones(prod(siz)/n,1);
y=nn';
nn(ncol)=nnan;
nn(~isnan(x(end,:)))=n;

ieven=(rem(nn,2)==0);           % Even number of elements along DIM
iodd=find(~ieven & ~isnan(nn)); % Odd number of elements along DIM
ieven=find(ieven==1);
x=x(:);
if ~isempty(iodd), 
    y(iodd)=x(sub2ind(sizx,(nn(iodd)+1)/2,iodd)); 
end;
if ~isempty(ieven), 
    y(ieven)=(x(sub2ind(sizx,nn(ieven)/2,ieven))+ x(sub2ind(sizx,nn(ieven)/2+1,ieven)))/2;
end;

% Permute and reshape back
siz(dim) = 1;
y = ipermute(reshape(y,siz(perm)),perm);
