function [num_ligne,varargout]=get_txtfile_col(ficname,D,select_ligne)
% -========================================================
%   USAGE : [num_ligne,col1,col2,col3,...coln]=get_txtfile_col(ficname,D,select_ligne)
%   PURPOSE : read a txt file with nbcol columns, one variable per colunm (column are delimited by  "white space" (default) or one of the characters in D. Any leading delimiter characters are ignored => use of the matlab function strtock)                               
% -----------------------------------
%   INPUT :
%     ficname   (string)  file(txt)  name
%            
%   OPTIONNAL INPUT :
%    D            (string)  column delimiter
%    select_ligne (string)  a character in the firt column of the file that is used to identify lignes that should be read 
% -----------------------------------
%   OUTPUT :
%     num_ligne (array of integer)   line number as it is in the file (all lines)
%     col1   (cell array of string)  size(1,nbligne) (nbligne :only the selected lines, nbligne <= num_ligne)
%             variable in the first column
%     col2   (cell array of string)  size(1,nbligne)
%            variable in the second column
%     ...
% -----------------------------------
%   EXAMPLE :
%   [temp,sal,wmoid]=get_txtfile_col('TSArgo.txt',',')       
%   (ex lig1) 22,35,1900265    => temp{1}='22' sal{1}='35'  wmoid{1}='1900265'
%   (ex lig2) 12,35.2,1900264  => temp{2}='12' sal{2}='35.1'  wmoid{2}='1900264'
%
%   [ficname,temp,sal,wmoid]=get_txtfile_col('TSArgo_new.txt',' ','.nc')
%   (ex lig1)  filename temp sal wmoid    => nothing is read
%   (ex lig2)  truc1.nc  22  35 1900265   => ficname{1}='truc1.nc' temp{1}='22' sal{1}='35'  wmoid{1}='1900265'
%   (ex lig3)  truc2.nc  12 35.2 1900264  => ficname{2}='truc2.nc' temp{2}='12' sal{2}='35.1'  wmoid{2}='1900264'
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%   CALLED SUBROUTINES: none
% ========================================================

nout=nargout;  % nombre de colonnes à lire (peut etre inferieur au nombre total de colonne)

if nargin==1
   D= ' ';
   select_ligne='';
end
if nargin==2
   select_ligne='';
end

except_ligne='#'; % ne lit pas les lignes debutant par '#'

icompt = 1; % determine la fin du fichier quand ==0
ilin=0;
nblin=0;
if exist(ficname,'file')==2
    fid=fopen(ficname);
    while icompt
	tline=fgetl(fid); % tline est de type char
	
	if ~ischar(tline) 
	    icompt=0;     % arrive a la fin du fichier
	else
	    if isspace(D)==0  % si jamais separateur de champ, remplace les separateurs accolés par deux separateur espacés
		
		tline= strrep(tline,[D D],[D ' ' D]);
		tline= strrep(tline,[D D],[D ' ' D]);
	    end
	    nblin=nblin+1;
	    
	    [part1,r]=strtok(tline,D); % coupe tline au niveau du premier espace (1ere col)  ' '=> part1 avant, r apres
	    if isempty(select_ligne)==1 % on lit toutes les lignes sans distinction 
		if strcmp(tline(1),except_ligne)==0 
		    ilin=ilin+1;            % compteur de ligne
		    num_ligne(ilin)=nblin;
		    col.l1{ilin}=part1; 
		    for nbcol=2:nout
			dynf=['l' num2str(nbcol)];
			[part1,r]=strtok(r,D);
			col.(dynf){ilin}=part1;
		    end
		end
	    else
		if findstr(part1,select_ligne)  % verifie que l'on est bien sur une ligne de modif de flag (qui commence par le nom du fichier netcdf)
		    if strcmp(tline(1),except_ligne)==0 
			ilin=ilin+1;         % compteur de ligne
			col.l1{ilin}=part1; 
			num_ligne(ilin)=nblin;
			for nbcol=2:nout
			    dynf=['l' num2str(nbcol)];
			    [part1,r]=strtok(r,D);
			    col.(dynf){ilin}=part1;
			end
		    end
		end
	    end
	end
    end
    fclose(fid);
else
    error(['le fichier ' ficname ' n''existe pas'])
end
for nbcol=1:nout
    dynf=['l' num2str(nbcol)];
    varargout{nbcol,:} = {col.(dynf){:}};
end