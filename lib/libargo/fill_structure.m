function [S,Dim]=fill_structure(var,varname,dimname,S,Dim)
% -========================================================
%   USAGE : [S,Dim]=fill_structure(var,varname,dimname,S,Dim)
%   PURPOSE : trouve une chaine de caractère dans un tableau de caractere
%   ou cell array
% -----------------------------------
%   INPUT :
%     IN1   (class)  -comments-
%             additional description
%     IN2   (class)  -comments-
%
%   OPTIONNAL INPUT :
%    OPTION1  (class)  -comments-
% -----------------------------------
%   OUTPUT :
%     OUT1   (class)  -comments-
%             additional description
%     OUT2   (class)  -comments-
%             additional description
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%            : modified (yyyy) byxxx 
%   CALLED SUBROUTINES: none
% ========================================================

if nargin<=3
    S=struct;
    Dim=struct;
end

% analyse des dimensions de var
m=size(var);
existdimname=fieldnames(Dim);
for l=1:length(dimname)
        ischar=[];
        if isempty(existdimname)==0
        h=strfind(existdimname,lower(dimname{l}));
        ischar=~cellfun('isempty',h);
        end
        if sum(ischar)==0
        % rajoute des dimensions à Dim
        Dim.(lower(dimname{l})).name=dimname{l};
        Dim.(lower(dimname{l})).dimlength=m(l);
        end
end

S.(lower(varname)).data=var;
S.(lower(varname)).name=upper(varname);
S.(lower(varname)).dim=dimname;
if isa(var,'double')
   S.(lower(varname)).type=6;
   S.(lower(varname)).FillValue_=99999;
end
if isa(var,'single')
    S.(lower(varname)).type=5;
    S.(lower(varname)).FillValue_=single(99999);
end
if isa(var,'integer')
    S.(lower(varname)).type=4;
    S.(lower(varname)).FillValue_=int16(99999);
end
if isa(var,'char')
   S.(lower(varname)).type=2;
   S.(lower(varname)).FillValue_=' ';
end
if isa(var,'logical')
   S.(lower(varname)).type=1;
end
S.(lower(varname))=orderfields(S.(lower(varname)),{'name','dim','data','FillValue_','type'});
    
S.obj='ObsInSitu';
