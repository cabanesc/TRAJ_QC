% ------------------------------------------------------------------------------
% Suppression des doublons de dates sur un ensemble de localisations Argos.
%
% SYNTAX :
% [o_posDate, o_posLon, o_posLat, o_posQc, o_posTemp] = ...
%    dep_delete_doubles(a_posDate, a_posLon, a_posLat, a_posQc, a_posTemp)
%
% INPUT PARAMETERS :
%   a_posDate  : dates des localisations Argos
%   a_posLon   : longitudes des localisations Argos
%   a_posLat   : latitudes des localisations Argos
%   a_posQc    : classes des localisations Argos
%   a_posTemp  : températures associées aux localisations Argos
%
% OUTPUT PARAMETERS :
%
% EXAMPLES :
%
% SEE ALSO : 
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   13/05/2008 - RNU - creation
%  modif G.Herbert
% ------------------------------------------------------------------------------
function [isdouble_toremove, idDoublon_date, idDoublon_lonlat] = ...
   find_doubles_datelatlon(a_posDate, a_posLon, a_posLat, a_posQc, a_posTemp)
global PARAM

o_posDate = a_posDate;
o_posLon = a_posLon;
o_posLat = a_posLat;
o_posQc = a_posQc;
o_posTemp = a_posTemp;
% on supprime les éventuels doublons de dates
idDelete = [];
idDoublon_date=[];
idDoublon_lonlat = [];
argosDates = unique(o_posDate);
% correction cc 24/09/2020 : attention les dates ne sont pas forcement
% ordonnees: on peut trouver un double plus loin dans la serie
[o_posDate_sorted,isort]=sort(o_posDate);
o_posLon_sorted=o_posLon(isort);
o_posLat_sorted=o_posLat(isort);
o_posQc_sorted=o_posQc(isort);
diffDates = diff(o_posDate_sorted);
%diffDates = diff(o_posDate);

idDoublon_d = find(diffDates<= PARAM.TIME_DIFF_DOUBLE_DATE_LOC);   %%considère un doublon lorsque la différence entre les deux positions est <=30s
 
   if (length(idDoublon_d) > 1)
     
     for id = 1:length(idDoublon_d)
          idDoublon=[idDoublon_d(id) idDoublon_d(id)+1];
          idDoublon_date=[idDoublon_date idDoublon];
          %% (récupère id quand doublons date associées à doublons lon-lat) --> non 
          if(diff(o_posLon_sorted(idDoublon))==0 & diff(o_posLat_sorted(idDoublon))==0)
              idDoublon_lonlat = [idDoublon_lonlat; idDoublon];
              % idDoublon_date=[idDoublon_date idDoublon];
          end
          % parmi les doublons, on conserve la localisation de meilleure classe Argos
          if (sum(isstrprop(o_posQc_sorted(idDoublon), 'digit')) ~= 0)
              % il n'y a que des classes 1 2 3
             [qcMax id] = max(o_posQc_sorted(idDoublon));
             idChoix = id(1);
          elseif (sum(isstrprop(o_posQc_sorted(idDoublon), 'alpha')) ~= 0)
              % il n'y a que des classes A B Z
             [qcMin id] = min(o_posQc_sorted(idDoublon));
             idChoix = id(1);
          else
             idDigit = find(isstrprop(o_posQc_sorted(idDoublon), 'digit') == 1);
             [qcMax id] = max(o_posQc_sorted(idDoublon(idDigit)));
             if (~isempty(id))
               idChoix = idDigit(id(1));
             else
               idChoix = 1; % cas ou on a aucune classe de loc pour décider
             end
          end
          idDoublon(idChoix) = [];
          if (size(idDoublon, 1) < size(idDoublon, 2))
             idDoublon = idDoublon';
          end
          idDelete = [idDelete; idDoublon];
     end
   end
%idDoublon_lonlat = idDoublon_lonlat(:);
idDoublon_lonlat = isort(idDoublon_lonlat(:)); % cc correction 24/09/2020
%idDoublon_date = idDoublon_date(:);
idDoublon_date = isort(idDoublon_date(:));
%isdouble_toremove=idDelete;
isdouble_toremove=isort(idDelete);

return;
