function [file_list]=select_float_files_on_ftp(floatname,dacname,DIRFTP,FileType,IncludeDescProf)
% -========================================================
%   USAGE : [file_list]=select_float_files_on_ftp(floatname,dacname,DIR_FTP,FileType)
%   PURPOSE : selectionne les fichiers mono profils sur le ftp argo en fonction du type de fichier souhaité:
%  Selections des fichiers mono-profil a lire: 
%  -  fichiers Merged ('M','MR' ou 'MD')  donnees PTS + param bio finaux
%  -  fichiers Bio    ('B', 'BR', ou 'BD') données P + bio (intermediaire + param finaux)
%  -  fichiers Core   ('C', 'CR', ou 'CD'); données PTS
%  -  IncludeDescProf (0 si on n'inclue pas les profils descendant, 1 sinon (defaut)
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
if nargin<5
IncludeDescProf=1;
end

%if isempty(dacname)==0
repnc = [DIRFTP dacname '/' floatname '/profiles/'];
%else
%repnc = [DIRFTP dacname '/' floatname '/'];
%end
if isempty(dacname)&isempty(floatname)
repnc = [DIRFTP];
end
list=dir([repnc '*.nc']);

%  Selections des fichiers mono-profil a lire: 
%  -  fichiers Merged ('M','MR' ou 'MD')  donnees PTS + param bio finaux
%  -  fichiers Bio    ('B', 'BR', ou 'BD') données P + bio (intermediaire + param finaux)
%  -  fichiers Core   ('C', 'CR', ou 'CD'); données PTS

  
%keyboard

filenames={list.name};
d=filenames;

d=strcat('C',d);

ischarcellMR = strfind(d,'CMR');
ischarlogMR  = ~cellfun(@isempty, ischarcellMR);
ischarcellMD = strfind(d,'CMD');
ischarlogMD  = ~cellfun(@isempty, ischarcellMD);
ischarcellBR = strfind(d,'CBR');
ischarlogBR  = ~cellfun(@isempty, ischarcellBR);
ischarcellBD = strfind(d,'CBD');
ischarlogBD  = ~cellfun(@isempty, ischarcellBD);
ischarcellD = strfind(d,'CD');
ischarlogD  = ~cellfun(@isempty, ischarcellD);
ischarcellR = strfind(d,'CR');
ischarlogR  = ~cellfun(@isempty, ischarcellR);


if FileType=='M'
   the_files=filenames(ischarlogMR|ischarlogMD);
end
if FileType=='B'
   the_files=filenames(ischarlogBR|ischarlogBD);
end
if FileType=='C'
   the_files=filenames(ischarlogR|ischarlogD);
end
if FileType=='MR'
   the_files=filenames(ischarlogMR);
end
if FileType=='BR'
   the_files=filenames(ischarlogBR);
end
if FileType=='CR'
   the_files=filenames(ischarlogR);
end
if FileType=='MD'
   the_files=filenames(ischarlogMD);
end
if FileType=='BD'
   the_files=filenames(ischarlogBD);
end
if FileType=='CD'
   the_files=filenames(ischarlogD);
end


% selection ou non des profils descendants
%repD = input('Prise en compte du profil D ? (o/n) ','s');

if IncludeDescProf==1
    Asc_files = the_files;
else   
% ne prend pas les profils descendants
    ischarcellD = strfind(the_files,'D.nc');

    if ~isempty(ischarcellD)% enleve descending profiles
         Asc_files = the_files(cellfun(@isempty,ischarcellD));
    end
end

% reordonne rep.name avec les profils descendants qui precedent les profils ascendants pour un meme cycle
% c'est plus simple pour les plots
d=Asc_files;
[t,r]=strtok(d,'_');
r=strrep(r,'.nc','Z.nc');
[rsorted,isortr]=sort(r);% range dans l'ordre alphabetique ASCII '10DZ.nc' precede '10Z.nc'
d=d(isortr);
file_list=d;   
if isempty(file_list)
    warning('******** No file selected, check your folder name')
end

