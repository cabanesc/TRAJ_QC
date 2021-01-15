% ------------------------------------------------------------------------------
% Test de vérification que les poisitions Argos ne sont pas trop agglomérées
%
% SYNTAX :
% o_idOk = check_argos_positions_aglo2(a_longitude, a_latitude, a_maxDist)
%
% INPUT PARAMETERS :
%   a_longitude : longitude des positions Argos observées
%   a_latitude  : latitude des positions Argos observées
%   a_posQc     : classe de localisation des positions Argos observées
%   a_maxDist   : distance maximale acceptable
%
% OUTPUT PARAMETERS :
%   o_idOk      : vaut 1 si au moins une position Agos est à plus de a_maxDist
%                 du point de référence (vaut 0 dans le cas contraire)
%   o_lonRef    : longitude de la position de référence (barycentre)
%   o_latRef    : longitude de la position de référence (barycentre)
%   o_lonCercle : longitudes des points du cercle englobant
%   o_latCercle : latitudes des points du cercle englobant
%
% EXAMPLES :
%
% SEE ALSO : dist
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   16/09/2008 - RNU - creation
% ------------------------------------------------------------------------------
function [o_idOk, o_lonRef, o_latRef, o_lonCercle, o_latCercle] = ...
   check_argos_positions_aglo2(a_longitude, a_latitude, a_posQc, a_maxDist)

o_idOk = [];
o_lonRef = [];
o_latRef = [];
o_lonCercle = [];
o_latCercle = [];

global g_latDef g_lonDef;

% initialisation des valeurs par défaut
init_valdef;

% on associe à chaque classe de localisation, la précision estimée en latitude
% et longitude fournie (en m) par Argos
precision(1) = 1000;
precision(2) = 350;
precision(3) = 150;

% precision affectée à une position GPS
precision(4) = 30;
a_posQc(find(a_posQc == 'G')) = '4';

% distance maximale acceptable
qcFirst = str2num(a_posQc(1));
qcLast = str2num(a_posQc(end));
if (qcFirst == 4)
   maxDist = a_maxDist(4);
elseif ((qcFirst == 1) && (qcLast == 1))
   maxDist = a_maxDist(1);
elseif ((qcFirst == 1) || (qcLast == 1))
   maxDist = a_maxDist(2);
else
   maxDist = a_maxDist(3);
end

nbPos = length(a_longitude);
tablon = ones(nbPos*2-1, 1)*g_lonDef;
tablat = ones(nbPos*2-1, 1)*g_latDef;

% détermination du point de référence choisi comme le  barycentre des positions
% Argos
% passage en repère cartésien
xPosOri = a_longitude(1);
yPosOri = a_latitude(1);

valCos = cosd(a_latitude(1))*1.852*60;
xPos = (a_longitude - xPosOri)*valCos;
yPos = (a_latitude - yPosOri)*1.852*60;

xPos = xPos./(precision(str2num(a_posQc'))/1000);
yPos = yPos./(precision(str2num(a_posQc'))/1000);

invPrec = precision(str2num(a_posQc'))/1000;
invPrec = 1./invPrec;

xRef = sum(xPos)/sum(invPrec);
yRef = sum(yPos)/sum(invPrec);

% on repasse en dégrés décimaux
lonRef = xRef/valCos;
lonRef = lonRef + xPosOri;
latRef = yRef/(1.852*60);
latRef = latRef + yPosOri;

o_lonRef = lonRef;
o_latRef = latRef;

% tableau de calcul des distances
tablon(1:2:length(tablon)) = a_longitude;
tablat(1:2:length(tablat)) = a_latitude;

tablon(2:2:length(tablon)) = ones(nbPos-1, 1)*lonRef;
tablat(2:2:length(tablat)) = ones(nbPos-1, 1)*latRef;

% tous les points sont-ils à moins de a_maxDist du point de référence ?
ranges = dist(tablat, tablon);
if (~isempty(find(ranges > maxDist, 1)))
   o_idOk = 1;
else
   o_idOk = 0;
end

% calcul des positions du cercle de recherche
xPoint = [];
yPoint = [];
for id = 0:1:360
   xPoint = [xPoint cosd(id)];
   yPoint = [yPoint sind(id)];
end

xPoint = xPoint*maxDist/1000;
yPoint = yPoint*maxDist/1000;

xPoint = xPoint/(cosd(o_latRef)*1.852*60);
xPoint = xPoint + o_lonRef;
yPoint = yPoint/(1.852*60);
yPoint = yPoint + o_latRef;

o_lonCercle = xPoint;
o_latCercle = yPoint;

return;

