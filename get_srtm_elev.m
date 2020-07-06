% ------------------------------------------------------------------------------
% R�cup�ration des �l�vations SRTM (et des positions associ�es) pour une zone
% donn�e. Les longitudes doivent �tre dans l'intervalle [-180; 180[ avec
% �ventuellement la ligne de changement de date dans la zone concern�e.
%
% SYNTAX :
%   function [o_elev, o_lon , o_lat] = get_srtm_elev(a_lonMin, a_lonMax,
%                                                    a_latMin, a_latMax)
%
% INPUT PARAMETERS :
%   a_lonMin     : longitude minimale de la zone concern�e
%   a_lonMax     : longitude maximale de la zone concern�e
%   a_latMin     : latitude minimale de la zone concern�e
%   a_latMax     : latitude maximale de la zone concern�e
%
% OUTPUT PARAMETERS :
%   o_elev : �l�vations r�cup�r�es
%   o_lon  : longitudes des �l�vations r�cup�r�es
%   o_lat  : latitudes des �l�vations r�cup�r�es
%
% EXAMPLES :
%
% SEE ALSO : init_valdef, get_srtm_data
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   20/11/2008 - RNU - creation
% ------------------------------------------------------------------------------
function [o_elev, o_lon , o_lat] = get_srtm_elev(a_lonMin, a_lonMax, a_latMin, a_latMax)

global g_elevDef;
DIR_BATHY = '/home/gherbert/andro/andro/data/bathy/SRTM30plus/data/';
% initialisation des valeurs par d�faut
init_valdef;

o_elev = [];
o_lon = [];
o_lat = [];

% contr�le des positions d�finissant la zone g�ographique demand�e
if (a_latMin > a_latMax)
   fprintf('*** ATTENTION *** get_srtm_elev: zone g�ographique incoh�rente!\n');
   return;
else
   if ((a_latMin < -90) || (a_latMax > 90))
      fprintf('*** ATTENTION *** get_srtm_elev: zone g�ographique incoh�rente!\n');
      return;
   end
end

if (a_lonMin >= 180)
   a_lonMin = a_lonMin - 360;
   a_lonMax = a_lonMax - 360;
end

if (a_lonMax < a_lonMin)
   a_lonMax = a_lonMax + 360;
end

% zone g�ographique minimale englobant la zone demand�e et correspondant
% aux points de grille (1/120 degr�)
lonMin = floor(a_lonMin*120)/120;
lonMax = ceil(a_lonMax*120)/120;
latMin = floor(a_latMin*120)/120;
latMax = ceil(a_latMax*120)/120;

% on ne g�re qu'un chevauchement de tuiles (soit une utilisation simultan�e d'au
% plus quatre tuiles)
if (abs(lonMax - lonMin) >= 20) || (abs(latMax - latMin) >= 30)
   fprintf('*** ATTENTION *** get_srtm_elev: zone g�ographique trop �tendue!\n');
   return;
end

% d�termination du cas de figure en latitude
if (latMin > 40); bandeLat = 1;
elseif (latMax > 40); bandeLat = 2;
elseif (latMin > -10); bandeLat = 3;
elseif (latMax > -10); bandeLat = 4;
elseif (latMin > -60); bandeLat = 5;
elseif (latMax > -60); bandeLat = 6;
else
   bandeLat = 7;
end

% d�termination du cas de figure en longitude
if (bandeLat <= 5)
   if (lonMax >= 180); bandeLon = 1;
   elseif (lonMin >= 140); bandeLon = 18;
   elseif (lonMax >= 140); bandeLon = 17;
   elseif (lonMin >= 100); bandeLon = 16;
   elseif (lonMax >= 100); bandeLon = 15;
   elseif (lonMin >= 60); bandeLon = 14;
   elseif (lonMax >= 60); bandeLon = 13;
   elseif (lonMin >= 20); bandeLon = 12;
   elseif (lonMax >= 20); bandeLon = 11;
   elseif (lonMin >= -20); bandeLon = 10;
   elseif (lonMax >= -20); bandeLon = 9;
   elseif (lonMin >= -60); bandeLon = 8;
   elseif (lonMax >= -60); bandeLon = 7;
   elseif (lonMin >= -100); bandeLon = 6;
   elseif (lonMax >= -100); bandeLon = 5;
   elseif (lonMin >= -140); bandeLon = 4;
   elseif (lonMax >= -140); bandeLon = 3;
   else
      bandeLon = 2;
   end
elseif (bandeLat == 6)
   if (lonMax >= 180); bandeLon = 1;
   elseif (lonMin >= 140); bandeLon = 24;
   elseif (lonMax >= 140); bandeLon = 23;
   elseif (lonMin >= 120); bandeLon = 22;
   elseif (lonMax >= 120); bandeLon = 21;
   elseif (lonMin >= 100); bandeLon = 20;
   elseif (lonMax >= 100); bandeLon = 19;
   elseif (lonMin >= 60); bandeLon = 18;
   elseif (lonMax >= 60); bandeLon = 17;
   elseif (lonMin >= 20); bandeLon = 16;
   elseif (lonMax >= 20); bandeLon = 15;
   elseif (lonMin >= 0); bandeLon = 14;
   elseif (lonMax >= 0); bandeLon = 13;
   elseif (lonMin >= -20); bandeLon = 12;
   elseif (lonMax >= -20); bandeLon = 11;
   elseif (lonMin >= -60); bandeLon = 10;
   elseif (lonMax >= -60); bandeLon = 9;
   elseif (lonMin >= -100); bandeLon = 8;
   elseif (lonMax >= -100); bandeLon = 7;
   elseif (lonMin >= -120); bandeLon = 6;
   elseif (lonMax >= -120); bandeLon = 5;
   elseif (lonMin >= -140); bandeLon = 4;
   elseif (lonMax >= -140); bandeLon = 3;
   else
      bandeLon = 2;
   end
else
   if (lonMax >= 180); bandeLon = 1;
   elseif (lonMin >= 120); bandeLon = 12;
   elseif (lonMax >= 120); bandeLon = 11;
   elseif (lonMin >= 60); bandeLon = 10;
   elseif (lonMax >= 60); bandeLon = 9;
   elseif (lonMin >= 0); bandeLon = 8;
   elseif (lonMax >= 0); bandeLon = 7;
   elseif (lonMin >= -60); bandeLon = 6;
   elseif (lonMax >= -60); bandeLon = 5;
   elseif (lonMin >= -120); bandeLon = 4;
   elseif (lonMax >= -120); bandeLon = 3;
   else
      bandeLon = 2;
   end
end

% r�cup�ration des donn�es
if (bandeLon > 1)
   if (bandeLat <= 5)
      if (rem(bandeLat, 2) == 1)
         if (rem(bandeLon, 2) == 0)
            fileNum = fix(bandeLat/2)*9 + fix(bandeLon/2);
           % [o_elev, o_lon , o_lat] = get_srtm_data(fileNum, lonMin, lonMax, latMin, latMax);
            [o_elev, o_lon , o_lat] = get_srtm_data2(DIR_BATHY,fileNum, lonMin, lonMax, latMin, latMax);
         else
            lonCut = -140 + 40*(fix(bandeLon/2)-1);
            fileNum1 = fix(bandeLat/2)*9 + fix(bandeLon/2);
            %[elev1, lon1, lat1] = get_srtm_data(fileNum1, lonMin, lonCut-1/120, latMin, latMax);
            [elev1, lon1 , lat1] = get_srtm_data2(DIR_BATHY,fileNum1, lonMin, lonCut-1/120, latMin, latMax);
            fileNum2 = fileNum1 + 1;
            %[elev2, lon2, lat2] = get_srtm_data(fileNum2, lonCut, lonMax, latMin, latMax);
            [elev2, lon2 , lat2] = get_srtm_data2(DIR_BATHY,fileNum2, lonCut, lonMax, latMin, latMax);

            o_elev = [elev1 elev2];
            o_lon = [lon1 lon2];
            o_lat = lat1;
         end
      else
         if (rem(bandeLon, 2) == 0)
            latCut = 40 - 50*(fix(bandeLat/2)-1);
            fileNum1 = (fix(bandeLat/2)-1)*9 + fix(bandeLon/2);
            %[elev1, lon1, lat1] = get_srtm_data(fileNum1, lonMin, lonMax, latCut+1/120, latMax);
            [elev1, lon1, lat1] = get_srtm_data2(DIR_BATHY,fileNum1, lonMin, lonMax, latCut+1/120, latMax);
            fileNum2 = fileNum1 + 9;
            %[elev2, lon2, lat2] = get_srtm_data(fileNum2, lonMin, lonMax, latMin, latCut);
            [elev2, lon2, lat2] = get_srtm_data2(DIR_BATHY,fileNum2, lonMin, lonMax, latMin, latCut);
            o_elev = [elev1; elev2];
            o_lon = lon1;
            o_lat = [lat1 lat2];
         else
            lonCut = -140 + 40*(fix(bandeLon/2)-1);
            latCut = 40 - 50*(fix(bandeLat/2)-1);
            fileNum1 = (fix(bandeLat/2)-1)*9 + fix(bandeLon/2);
            %[elev1, lon1, lat1] = get_srtm_data(fileNum1, lonMin, lonCut-1/120, latCut+1/120, latMax);
            [elev1, lon1, lat1] = get_srtm_data2(DIR_BATHY,fileNum1, lonMin, lonCut-1/120, latCut+1/120, latMax);
            fileNum2 = fileNum1 + 1;
            %[elev2, lon2, lat2] = get_srtm_data(fileNum2, lonCut, lonMax, latCut+1/120, latMax);
            [elev2, lon2, lat2] = get_srtm_data2(DIR_BATHY,fileNum2, lonCut, lonMax, latCut+1/120, latMax);
            fileNum3 = fileNum1 + 9;
            %[elev3, lon3, lat3] = get_srtm_data(fileNum3, lonMin, lonCut-1/120, latMin, latCut);
            [elev3, lon3, lat3] = get_srtm_data2(DIR_BATHY,fileNum3, lonMin, lonCut-1/120, latMin, latCut);
            fileNum4 = fileNum2 + 9;
           % [elev4, lon4, lat4] = get_srtm_data(fileNum4, lonCut, lonMax, latMin, latCut);
            [elev4, lon4, lat4] = get_srtm_data2(DIR_BATHY,fileNum4, lonCut, lonMax, latMin, latCut);
            o_elev = [elev1 elev2; elev3 elev4];
            o_lon = [lon1 lon2];
            o_lat = [lat1 lat3];
         end
      end
   elseif (bandeLat == 6)
      if (rem(bandeLon, 2) == 0)
         tabFileNum = [19 20 20 21 22 23 23 24 25 26 26 27];
         latCut = -60;
         fileNum1 = tabFileNum(round(bandeLon/2));
         %[elev1, lon1, lat1] = get_srtm_data(fileNum1, lonMin, lonMax, latCut+1/120, latMax);
         [elev1, lon1, lat1] = get_srtm_data2(DIR_BATHY,fileNum1, lonMin, lonMax, latCut+1/120, latMax);
         tabBandeLon = [2 4 6 8 10 12 14 16 18 20 22 24];
         tabFileNum = [28 28 29 29 30 30 31 31 32 32 33 33];
         fileNum2 = tabFileNum(find(tabBandeLon == bandeLon));
         %[elev2, lon2, lat2] = get_srtm_data(fileNum2, lonMin, lonMax, latMin, latCut);
         [elev2, lon2, lat2] = get_srtm_data2(DIR_BATHY,fileNum2, lonMin, lonMax, latMin, latCut);
         o_elev = [elev1; elev2];
         o_lon = lon1;
         o_lat = [lat1 lat2];
      else
         tabBandeLon = [3 5 7 9 11 13 15 17 19 21 23];
         tabFileNum1 = [19 20 20 21 22 23 23 24 25 26 26];
         tabFileNum2 = [20  0 21 22 23  0 24 25 26  0 27];
         tabFileNum3 = [28 28 29 29 30 30 31 31 32 32 33];
         tabFileNum4 = [ 0 29  0 30  0 31  0 32  0 33  0];
         fileNum1 = tabFileNum1(find(tabBandeLon == bandeLon));
         fileNum2 = tabFileNum2(find(tabBandeLon == bandeLon));
         fileNum3 = tabFileNum3(find(tabBandeLon == bandeLon));
         fileNum4 = tabFileNum4(find(tabBandeLon == bandeLon));

         if (fileNum2 == 0)
            tabBandeLon = [5 13 21];
            tabLonCut = [-120 0 120];
            lonCut = tabLonCut(find(tabBandeLon == bandeLon));
            latCut = -60;
            %[elev1, lon1, lat1] = get_srtm_data(fileNum1, lonMin, lonMax, latCut+1/120, latMax);
            [elev1, lon1, lat1] = get_srtm_data2(DIR_BATHY,fileNum1, lonMin, lonMax, latCut+1/120, latMax);
            %[elev3, lon3, lat3] = get_srtm_data(fileNum3, lonMin, lonCut-1/120, latMin, latCut);
            [elev3, lon3, lat3] = get_srtm_data2(DIR_BATHY,fileNum3, lonMin, lonCut-1/120, latMin, latCut);
            %[elev4, lon4, lat4] = get_srtm_data(fileNum4, lonCut, lonMax, latMin, latCut);
            [elev4, lon4, lat4] = get_srtm_data2(DIR_BATHY,fileNum4, lonCut, lonMax, latMin, latCut);
            o_elev = [elev1; elev3 elev4];
            o_lon = lon1;
            o_lat = [lat1 lat3];
         elseif (fileNum4 == 0)
            tabBandeLon = [3 7 11 15 19 23];
            tabLonCut = [-140 -100 -20 20 100 140];
            lonCut = tabLonCut(find(tabBandeLon == bandeLon));
            latCut = -60;
            %[elev1, lon1, lat1] = get_srtm_data(fileNum1, lonMin, lonCut-1/120, latCut+1/120, latMax);
            [elev1, lon1, lat1] = get_srtm_data2(DIR_BATHY,fileNum1, lonMin, lonCut-1/120, latCut+1/120, latMax);
            %[elev2, lon2, lat2] = get_srtm_data(fileNum2, lonCut, lonMax, latCut+1/120, latMax);
            [elev2, lon2, lat2] = get_srtm_data2(DIR_BATHY,fileNum2, lonCut, lonMax, latCut+1/120, latMax);
            %[elev3, lon3, lat3] = get_srtm_data(fileNum3, lonMin, lonMax, latMin, latCut);
            [elev3, lon3, lat3] = get_srtm_data2(DIR_BATHY,fileNum3, lonMin, lonMax, latMin, latCut);
            o_elev = [elev1 elev2; elev3];
            o_lon = lon3;
            o_lat = [lat1 lat3];
         else
            tabBandeLon = [9 17];
            tabLonCut = [-60 60];
            lonCut = tabLonCut(find(tabBandeLon == bandeLon));
            latCut = -60;
            %[elev1, lon1, lat1] = get_srtm_data(fileNum1, lonMin, lonCut-1/120, latCut+1/120, latMax);
            [elev1, lon1, lat1] = get_srtm_data2(DIR_BATHY,fileNum1, lonMin, lonCut-1/120, latCut+1/120, latMax);
            %[elev2, lon2, lat2] = get_srtm_data(fileNum2, lonCut, lonMax, latCut+1/120, latMax);
            [elev2, lon2, lat2] = get_srtm_data2(DIR_BATHY,fileNum2, lonCut, lonMax, latCut+1/120, latMax);
            %[elev3, lon3, lat3] = get_srtm_data(fileNum3, lonMin, lonCut-1/120, latMin, latCut);
            [elev3, lon3, lat3] = get_srtm_data2(DIR_BATHY,fileNum3, lonMin, lonCut-1/120, latMin, latCut);
            %[elev4, lon4, lat4] = get_srtm_data(fileNum4, lonCut, lonMax, latMin, latCut);
            [elev4, lon4, lat4] = get_srtm_data2(DIR_BATHY,fileNum4, lonCut, lonMax, latMin, latCut);

            o_elev = [elev1 elev2; elev3 elev4];
            o_lon = [lon1 lon2];
            o_lat = [lat1 lat3];
         end
      end
   else
      if (rem(bandeLon, 2) == 0)
         fileNum = fix(bandeLon/2) + 27;
         %[o_elev, o_lon , o_lat] = get_srtm_data(fileNum, lonMin, lonMax, latMin, latMax);
         [o_elev, o_lon , o_lat] = get_srtm_data2(DIR_BATHY,fileNum, lonMin, lonMax, latMin, latMax);
      else
         lonCut = -120 + 60*(fix(bandeLon/2)-1);
         fileNum1 = fix(bandeLon/2) + 27;
         %[elev1, lon1, lat1] = get_srtm_data(fileNum1, lonMin, lonCut-1/120, latMin, latMax);
         [elev1, lon1, lat1] = get_srtm_data2(DIR_BATHY,fileNum1, lonMin, lonCut-1/120, latMin, latMax);
         fileNum2 = fileNum1 + 1;
        % [elev2, lon2, lat2] = get_srtm_data(fileNum2, lonCut, lonMax, latMin, latMax);(
         [elev2, lon2, lat2] = get_srtm_data2(DIR_BATHY,fileNum2, lonCut, lonMax, latMin, latMax);
         o_elev = [elev1 elev2];
         o_lon = [lon1 lon2];
         o_lat = lat1;
      end
   end
else
   % cas o� la ligne de changement de date fait partie de la zone g�ographique
   % demand�e
   if (bandeLat <= 5)
      if (rem(bandeLat, 2) == 1)
         lonCut = -140 + 40*(fix(bandeLon/2)-1);
         fileNum2 = fix(bandeLat/2)*9 + 1;
         %[elev2, lon2, lat2] = get_srtm_data(fileNum2, lonCut, lonMax-360, latMin, latMax);
         [elev2, lon2, lat2] = get_srtm_data2(DIR_BATHY,fileNum2, lonCut, lonMax-360, latMin, latMax);
         fileNum1 = fileNum2 + 8;
         %[elev1, lon1, lat1] = get_srtm_data(fileNum1, lonMin, 360+lonCut-1/120, latMin, latMax);
         [elev1, lon1, lat1] = get_srtm_data2(DIR_BATHY,fileNum1, lonMin, 360+lonCut-1/120, latMin, latMax);
         o_elev = [elev1 elev2];
         o_lon = [lon1 lon2+360];
         o_lat = lat1;
      else
         lonCut = -140 + 40*(fix(bandeLon/2)-1);
         latCut = 40 - 50*(fix(bandeLat/2)-1);
         fileNum2 = (fix(bandeLat/2)-1)*9 + 1;
        % [elev2, lon2, lat2] = get_srtm_data(fileNum2, lonCut, lonMax-360, latCut+1/120, latMax);
         [elev2, lon2, lat2] = get_srtm_data2(DIR_BATHY,fileNum2, lonCut, lonMax-360, latCut+1/120, latMax);
         fileNum1 = fileNum2 + 8;
         %[elev1, lon1, lat1] = get_srtm_data(fileNum1, lonMin, 360+lonCut-1/120, latCut+1/120, latMax);
         [elev1, lon1, lat1] = get_srtm_data2(DIR_BATHY,fileNum1, lonMin, 360+lonCut-1/120, latCut+1/120, latMax);
         fileNum4 = fileNum2 + 9;
         %[elev4, lon4, lat4] = get_srtm_data(fileNum4, lonCut, lonMax-360, latMin, latCut);
         [elev4, lon4, lat4] = get_srtm_data2(DIR_BATHY,fileNum4, lonCut, lonMax-360, latMin, latCut);
         fileNum3 = fileNum1 + 9;
         %[elev3, lon3, lat3] = get_srtm_data(fileNum3, lonMin, 360+lonCut-1/120, latMin, latCut);
         [elev3, lon3, lat3] = get_srtm_data2(DIR_BATHY,fileNum3, lonMin, 360+lonCut-1/120, latMin, latCut);

         o_elev = [elev1 elev2; elev3 elev4];
         o_lon = [lon1 lon2+360];
         o_lat = [lat1 lat3];
      end
   elseif (bandeLat == 6)
      lonCut = -140 + 40*(fix(bandeLon/2)-1);
      latCut = -60;
      fileNum2 = 19;
      %[elev2, lon2, lat2] = get_srtm_data(fileNum2, lonCut, lonMax-360, latCut+1/120, latMax);
      [elev2, lon2, lat2] = get_srtm_data2(DIR_BATHY,fileNum2, lonCut, lonMax-360, latCut+1/120, latMax);
      fileNum1 = 27;
      %[elev1, lon1, lat1] = get_srtm_data(fileNum1, lonMin, 360+lonCut-1/120, latCut+1/120, latMax);
      [elev1, lon1, lat1] = get_srtm_data2(DIR_BATHY,fileNum1, lonMin, 360+lonCut-1/120, latCut+1/120, latMax);
      fileNum4 = 28;
      %[elev4, lon4, lat4] = get_srtm_data(fileNum4, lonCut, lonMax-360, latMin, latCut);
      [elev4, lon4, lat4] = get_srtm_data2(DIR_BATHY,fileNum4, lonCut, lonMax-360, latMin, latCut);
      fileNum3 = 33;
     % [elev3, lon3, lat3] = get_srtm_data(fileNum3, lonMin, 360+lonCut-1/120, latMin, latCut);
      [elev3, lon3, lat3] = get_srtm_data2(DIR_BATHY,fileNum3, lonMin, 360+lonCut-1/120, latMin, latCut);

      o_elev = [elev1 elev2; elev3 elev4];
      o_lon = [lon1 lon2+360];
      o_lat = [lat1 lat3];
   else
      lonCut = -140 + 40*(fix(bandeLon/2)-1);
      fileNum2 = 28;
      %[elev2, lon2, lat2] = get_srtm_data(fileNum2, lonCut, lonMax-360, latMin, latMax);
      [elev2, lon2, lat2] = get_srtm_data2(DIR_BATHY,fileNum2, lonCut, lonMax-360, latMin, latMax);
      fileNum1 = 33;
      %[elev1, lon1, lat1] = get_srtm_data(fileNum1, lonMin, 360+lonCut-1/120, latMin, latMax);
      [elev1, lon1, lat1] = get_srtm_data2(DIR_BATHY,fileNum1, lonMin, 360+lonCut-1/120, latMin, latMax);

      o_elev = [elev1 elev2];
      o_lon = [lon1 lon2+360];
      o_lat = lat1;
   end
end

% on v�rifie que le tableau des �l�vations est enti�rement rempli
if (~isempty(find(o_elev == g_elevDef, 1)))
   fprintf('*** ATTENTION *** get_srtm_elev: tableau des �l�vations incomplet!\n');
end

return;
