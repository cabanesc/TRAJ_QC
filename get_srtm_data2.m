% ------------------------------------------------------------------------------
% Lecture des �l�vations (et des positions associ�es) dans un fichier STRM30.
%
% SYNTAX :
%   function [o_elev, o_lon , o_lat] = ...
%      get_srtm_data(DIR_BATHY, a_fileNumber, a_lonMin, a_lonMax, a_latMin, a_latMax)
%
% INPUT PARAMETERS :
%   a_fileNumber : num�ro du fichier SRTM concern�
%   a_lonMin     : longitude minimale de la zone concern�e
%   a_lonMax     : longitude maximale de la zone concern�e
%   a_latMin     : latitude minimale de la zone concern�e
%   a_latMax     : latitude maximale de la zone concern�e
%
% OUTPUT PARAMETERS :
%   o_elev : �l�vations collect�es
%   o_lon  : longitudes des �l�vations collect�es
%   o_lat  : latitudes des �l�vations collect�es
%
% EXAMPLES :
%
% SEE ALSO : 
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   20/11/2008 - RNU - creation
% ------------------------------------------------------------------------------
function [o_elev, o_lon , o_lat] = ...
   get_srtm_data2(DIR_BATHY, a_fileNumber, a_lonMin, a_lonMax, a_latMin, a_latMax)

global g_elevDef;

% initialisation des valeurs par d�faut
init_valdef;

o_elev = [];
o_lon = [];
o_lat = [];

% fichiers binaires relatifs aux donn�es des diff�rentes tuiles
fileName{1} = 'w180n90';
fileName{2} = 'w140n90';
fileName{3} = 'w100n90';
fileName{4} = 'w060n90';
fileName{5} = 'w020n90';
fileName{6} = 'e020n90';
fileName{7} = 'e060n90';
fileName{8} = 'e100n90';
fileName{9} = 'e140n90';
fileName{10} = 'w180n40';
fileName{11} = 'w140n40';
fileName{12} = 'w100n40';
fileName{13} = 'w060n40';
fileName{14} = 'w020n40';
fileName{15} = 'e020n40';
fileName{16} = 'e060n40';
fileName{17} = 'e100n40';
fileName{18} = 'e140n40';
fileName{19} = 'w180s10';
fileName{20} = 'w140s10';
fileName{21} = 'w100s10';
fileName{22} = 'w060s10';
fileName{23} = 'w020s10';
fileName{24} = 'e020s10';
fileName{25} = 'e060s10';
fileName{26} = 'e100s10';
fileName{27} = 'e140s10';
fileName{28} = 'w180s60';
fileName{29} = 'w120s60';
fileName{30} = 'w060s60';
fileName{31} = 'w000s60';
fileName{32} = 'e060s60';
fileName{33} = 'e120s60';
currentFileName  = fileName{a_fileNumber};

test=ismember(fileName,currentFileName);
a_fileNumber = find(test==1);

% chemin d'acc�s au fichier demand�
%if exist('config_perso.txt','file')==0
%copyfile('config.txt','config_perso.txt')
%end
%C_FILE=load_configuration('config_perso.txt');

%srtmDirName = C_FILE.DIR_SRTM
srtmFileName = [DIR_BATHY currentFileName '.Bathymetry.srtm'];
%srtmFileName = currentFileName;

% ouverture du fichier (big-endian format)
fId = fopen(srtmFileName, 'r', 'b');
if (fId == -1)
   fprintf('Erreur ouverture fichier : %s/n', srtmFileName);
   return;
end

% nombre de points par ligne du fichier
nbPtsLine = 4800;
if (a_fileNumber > 27)
   nbPtsLine = 7200;
end

% coordonn�es du premier point du fichier
factLon = 1;
if (currentFileName(1) == 'w')
   factLon = -1;
end
ficLonMin = str2num(currentFileName(2:4))*factLon;

factLat = 1;
if (currentFileName(5) == 's')
   factLat = -1;
end
ficLatMax = str2num(currentFileName(6:7))*factLat;

% lecture des donn�es dans le fichier demand�
nbLon = round((a_lonMax-a_lonMin)*120+1);
nbLat = round((a_latMax-a_latMin)*120+1);
elev = ones(nbLat, nbLon)*g_elevDef;
%keyboard
idStart = round((ficLatMax-a_latMax)*120*nbPtsLine + (a_lonMin-ficLonMin)*120);
for id = 1:nbLat
   Val=[];
   fseek(fId, (idStart+(id-1)*nbPtsLine)*2, 'bof');   %% déplacement du curseur dans le fichier (bof: "begin of file")
   Val = fread(fId, [1 nbLon], 'int16');
      if(isempty(Val)==0)
        elev(id, 1:length(Val)) = Val;
      %else
      %  elev(id,:) = g_elevDef;
      end
end

lon = [round(a_lonMin*120):round(a_lonMax*120)]/120;
lat = [round(a_latMax*120):-1:round(a_latMin*120)]/120;

o_elev = elev;
o_lon = lon;
o_lat = lat;

%keyboard
fclose(fId);

return;
