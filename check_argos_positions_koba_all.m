% ------------------------------------------------------------------------------
% Passage du critère JAMSTEC sur un ensemble de positions Argos de surface
% (avec prise en compte du flag à 4).
%
% SYNTAX :
%  [o_date, o_longitude, o_latitude, o_posAcc, o_posQcIn, o_posQcOut ...
%    o_idBadPos] = check_argos_positions_koba_all( ...
%    a_date, a_longitude, a_latitude, a_posAcc, a_posQc, ...
%    a_juldLast, a_longLast, a_latLast, a_fidOut)
%
% INPUT PARAMETERS :
%   a_date      : date des positions Argos observées
%   a_longitude : longitude des positions Argos observées
%   a_latitude  : latitude des positions Argos observées
%   a_posAcc    : classe des positions Argos observées
%   a_posQc     : flag des positions Argos observées
%   a_juldLast  : date de la dernière bonne positions Argos du cycle
%                 précédent
%   a_longLast  : longitude de la dernière bonne positions Argos du cycle
%                 précédent
%   a_latLast   : latitude de la dernière bonne positions Argos du cycle
%                 précédent
%   a_fidOut    : Id du fichier log
%
% OUTPUT PARAMETERS :
%   o_date      : date des positions Argos (éventuellement reclassées)
%   o_longitude : longitude des positions Argos (éventuellement reclassées)
%   o_latitude  : latitude des positions Argos (éventuellement reclassées)
%   o_posAcc    : classe des positions Argos (éventuellement reclassées)
%   o_posQcIn   : flag original des positions Argos (éventuellement reclassées)
%   o_posQcOut  : flag final des positions Argos (éventuellement reclassées)
%   o_idBadPos  : indices des positions Argos rejetées (flag 3 ou 4)
%
% EXAMPLES :
%
% SEE ALSO : dist, julian_2_gregorian
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   12/11/2012 - RNU - creation
% ------------------------------------------------------------------------------
function [o_date, o_longitude, o_latitude, o_posAcc, o_posQcIn, o_posQcOut ...
   o_idBadPos] = check_argos_positions_koba_all( ...
   a_date, a_longitude, a_latitude, a_posAcc, a_posQc, ...
   a_juldLast, a_longLast, a_latLast, a_fidOut)

global g_dateDef;

% initialisation des valeurs par défaut
init_valdef;

global g_verbose;
global g_cycleNumber;

o_idBadPos = [];

% vitesse maximale autorisée (m/s)
MAX_VEL = 3;

% for idP = 1:length(a_date)
%    a_date(idP) = juld_jamstec(a_date(idP));
% end

% on convertit les classes Argos en numérique
posAccNum = ones(length(a_posAcc), 1)*-1;
posAccNum(find(a_posAcc == 'G')) = 1;
posAccNum(find(a_posAcc == 'D')) = 1;
posAccNum(find(a_posAcc == 'E')) = 1;
posAccNum(find(a_posAcc == 'F')) = 1;
posAccNum(find(a_posAcc == 'H')) = 1;
posAccNum(find(a_posAcc == '3')) = 2;
posAccNum(find(a_posAcc == '2')) = 3;
posAccNum(find(a_posAcc == '1')) = 4;
posAccNum(find(a_posAcc == '0')) = 5;
posAccNum(find(a_posAcc == 'A')) = 6;
posAccNum(find(a_posAcc == 'B')) = 6;
posAccNum(find(a_posAcc == 'Z')) = 6;
posAccNum(find(a_posAcc == 'I')) = 6;  % iridium
posAccNum(find(a_posAcc == 'R')) = 4;  % rafos
posAccNum(find(a_posAcc == 'U')) = 6;
posAccNum(find(a_posAcc == ' ')) = 6;  % ajout cc 29/04/2022
% on associe à chaque classe de localisation, la précision estimée en latitude
% et longitude fournie (en m) par Argos
precision(4) = 1000;
precision(3) = 350;
precision(2) = 150;

% precision affectée à une position GPS
precision(1) = 30;

% precision affectée aux positions "0", "A", "B" ou "Z"
precision(5) = 1500;
precision(6) = 5000;

% classement temporel des données
[date, idSort] = sort(a_date);
longitude = a_longitude(idSort);
latitude = a_latitude(idSort);
posAccStr = a_posAcc(idSort);
posAccNum = posAccNum(idSort);
posQc = [];
if (~isempty(a_posQc))
   posQc = a_posQc(idSort);
end
flag = ones(length(date), 1);
order = 1:length(date);

driftSpeed = ones(length(date), 1)*-1;
driftDst = ones(length(date), 1)*-1;
distErr = ones(length(date), 1)*-1;

n = 1;
m = 1;
%keyboard

while (1)

   flagged = -9;
   exeptid = 0;
   maxSpeed = 0;
   maxId = 0;
   distSum = 0;

   distance = ones(length(date), 1)*-1;
   deltaTime = ones(length(date), 1)*-1;
   speed = ones(length(date), 1)*-1;

   % calculation of subsurface drift speed
   if ((a_juldLast ~= g_dateDef) && (length(order) > 0))
      distance(order(1)) = distance_lpo([a_latLast latitude(order(1))], [a_longLast longitude(order(1))]);
      deltaTime(order(1)) = abs(date(order(1)) - a_juldLast);
      speed(order(1)) = distance(order(1))/(deltaTime(order(1))*86400);
      if (m == 1)
         driftSpeed(order(1)) = speed(order(1));
         driftDst(order(1)) = distance(order(1))/1000;
      elseif (speed(order(1)) > MAX_VEL)
         flag(order(1)) = 4;
         order(1) = [];
         m = m + 1;
         continue;
      end
   elseif (isempty(order))
      if (~isempty(a_fidOut))
         fprintf(a_fidOut, ' %d ---- ---- .. no data\n', n);
      end
      break
   end

   % calculation of surface drift speed
   for ii = 2:length(order)
      im1 = order(ii-1);
      i0 = order(ii);
      distance(i0) = distance_lpo([latitude(im1) latitude(i0)], [longitude(im1) longitude(i0)]);
      distSum = distSum + distance(i0);
      deltaTime(i0) = date(i0) - date(im1);
      if (deltaTime(i0) < 1/86400)
         deltaTime(i0) = 1/86400;
      end
      speed(i0) = distance(i0)/(deltaTime(i0)*86400);

      if (m == 1)
         driftSpeed(i0) = speed(i0);
         driftDst(i0) = distance(i0)/1000;
      else
         % deletion of duplicated data
         if ((latitude(im1) == latitude(i0)) && (longitude(im1) == longitude(i0)) && ...
               ((date(i0) == date(im1))))
            flag(i0) = 4;
            exeptid = exeptid + 1;
         elseif (deltaTime(i0) > 1)
            % deletion of another cycle data
            flag(i0) = 4;
            exeptid = exeptid + 1;
         end
      end

      if (flag(i0) == 1)
         speed(i0) = abs(speed(i0));
         if (speed(i0) > maxSpeed)
            maxSpeed = speed(i0);
            maxId = ii - exeptid;
         end
      end
   end

   if ((exeptid > 0) || (m == 1))
      idDel = find(flag(order) == 4);
      order(idDel) = [];
      if (isempty(order))
         if (~isempty(a_fidOut))
            fprintf(a_fidOut, ' %d ---- ---- .. no data\n', n);
         end
         break
      end
      m = m + 1;
      continue
   end

   if (length(order) == 1)
      maxSpeed = speed(order(1));
      avgSpeed = speed(order(1));
   else
      avgSpeed = distSum/((date(order(end)) - date(order(1)))*86400);
   end

   if (~isempty(a_fidOut))
      fprintf(a_fidOut, '%2d %4.2f %4.2f ..%3dpoints\n', n, avgSpeed, maxSpeed, length(order));
   end

   % discrimination of abnormal position
   if (maxSpeed > MAX_VEL)
      posAcc0 = posAccNum(order(maxId-1));
      posAcc1 = posAccNum(order(maxId));
      prec0 = precision(posAcc0);
      prec1 = precision(posAcc1);

      % case of another Argos class
      if (posAcc0 ~= posAcc1)
         if (posAcc0 > posAcc1)
            flagged = maxId - 1;
         else
            flagged = maxId;
         end
      else
         % case of same Argos class

         if (length(order) == 2)
            % 2 points at same cycle

         elseif (length(order) == 3)
            % 3 points at same cycle

            if (maxId == 2)
               ip1 = order(maxId+1);
               im1 = order(maxId-1);
               i0 = order(maxId);
               distanceF = distance_lpo([latitude(i0) latitude(ip1)], [longitude(i0) longitude(ip1)]);
               deltaTimeF = date(ip1) - date(i0);
               distanceL = distance_lpo([latitude(im1) latitude(ip1)], [longitude(im1) longitude(ip1)]);
               deltaTimeL = date(ip1) - date(im1);
            else % maxId == 3
               im1 = order(maxId-1);
               im2 = order(maxId-2);
               i0 = order(maxId);
               distanceF = distance_lpo([latitude(im2) latitude(i0)], [longitude(im2) longitude(i0)]);
               deltaTimeF = date(i0) - date(im2);
               distanceL = distance_lpo([latitude(im2) latitude(im1)], [longitude(im2) longitude(im1)]);
               deltaTimeL = date(im1) - date(im2);
            end

            spdF = distanceF/(deltaTimeF*86400);
            spdL = distanceL/(deltaTimeL*86400);
            if (spdF > spdL)
               flagged = maxId;
            else
               flagged = maxId - 1;
            end

         elseif (length(order) > 3)
            % 4 points or more at same cycle

            if (maxId == 2) % maxspeed is the first section
               ip2 = order(maxId+2);
               ip1 = order(maxId+1);
               distance1F = distance_lpo([latitude(ip1) latitude(ip2)], [longitude(ip1) longitude(ip2)]);
               deltaTime1F = date(ip2) - date(ip1);
               distance1L = distance1F;
               deltaTime1L = deltaTime1F;
            elseif (maxId > 2)
               im1 = order(maxId-1);
               im2 = order(maxId-2);
               i0 = order(maxId);
               distance1F = distance_lpo([latitude(im2) latitude(i0)], [longitude(im2) longitude(i0)]);
               deltaTime1F = date(i0) - date(im2);
               distance1L = distance_lpo([latitude(im2) latitude(im1)], [longitude(im2) longitude(im1)]);
               deltaTime1L = date(im1) - date(im2);
            else
               fprintf('ERROR: maxId error! maxId = %d\n', maxId);
               return;
            end

            if (maxId == length(order)) % maxspeed is the end section
               im3 = order(maxId-3);
               im2 = order(maxId-2);
               distance2F = distance_lpo([latitude(im3) latitude(im2)], [longitude(im3) longitude(im2)]);
               deltaTime2F = date(im2) - date(im3);
               distance2L = distance2F;
               deltaTime2L = deltaTime2F;
            elseif (maxId < length(order))
               im1 = order(maxId-1);
               i0 = order(maxId);
               ip1 = order(maxId+1);
               distance2F = distance_lpo([latitude(i0) latitude(ip1)], [longitude(i0) longitude(ip1)]);
               deltaTime2F = date(ip1) - date(i0);
               distance2L = distance_lpo([latitude(im1) latitude(ip1)], [longitude(im1) longitude(ip1)]);
               deltaTime2L = date(ip1) - date(im1);
            else
               fprintf('ERROR: maxId error! maxId = %d\n', maxId);
               return;
            end

            % comparing speeds of 2 routes
            spdF = (distance1F + distance2F)/((deltaTime1F + deltaTime2F)*86400);
            spdL = (distance1L + distance2L)/((deltaTime1L + deltaTime2L)*86400);
            if (spdF > spdL)
               flagged = maxId;
            else
               flagged = maxId - 1;
            end
         end
      end

      % calculation of distance error
      dderr = 1*sqrt(prec0*prec0 + prec1*prec1);
      distErr(order(maxId)) = dderr/1000;

      % flagging at abnormal point
      if (distance(order(maxId)) >= dderr)
         if (length(order) == 2)
            flag(order(maxId-1)) = 3;
            flag(order(maxId)) = 3;
         else
            flag(order(flagged)) = 3;
         end
      end
      if (length(order) == 2)
         break
      end

      order(flagged) = [];
   end

   if (flagged < 0)
      break
   end

   n = n + 1;
   m = m + 1;
end

if (~isempty(a_fidOut))
   if (~isempty(date))
      fprintf(a_fidOut, '   Lat    Lon     Date       Time    cls org new    speed  distance_lpo(km) distErr\n');
      for idP = 1:length(date)
         driftSpeedStr = [];
         if (driftSpeed(idP) ~= -1)
            driftSpeedStr = sprintf('%7.2f', driftSpeed(idP));
         end
         driftDstStr = [];
         if (driftDst(idP) ~= -1)
            driftDstStr = sprintf('%7.3f', driftDst(idP));
         end
         distErrStr = [];
         if (distErr(idP) ~= -1)
            distErrStr = sprintf('%5.3f', distErr(idP));
         end

         fprintf(a_fidOut, ' %5.3f  %6.3f %19s %3s %3s %3s  %8s %8s %6s\n', ...
            latitude(idP), longitude(idP), ...
            julian_2_gregorian(date(idP)), ...
            posAccStr(idP), posQc(idP), num2str(flag(idP)), ...
            driftSpeedStr, driftDstStr, distErrStr);
      end
   end
end

% output parameters
o_date = date;
o_longitude = longitude;
o_latitude = latitude;
o_posAcc = posAccStr;
o_posQcIn = posQc;
o_posQcOut = flag;
o_idBadPos = find(o_posQcOut ~= 1);

if (~isempty(o_idBadPos))
   if (g_verbose == 1)
      fprintf('KOBAYASHI cycle #%d exclusion des positions Argos de surface:', g_cycleNumber);
      fprintf(' #%d', o_idBadPos);
      fprintf('\n');
   end
end
%keyboard
return;

% % calcul de distance par JAMSTEC (pour effectuer les comparaisons)
% function [o_distance] = distance_jamstec(a_lat, a_lon)
% 
% R_EARTH = 6378.140;
% 
% lat1 = a_lat(1);
% lat2 = a_lat(2);
% lon1 = a_lon(1);
% lon2 = a_lon(2);
% 
% if ((lat1 == lat2) && (lon1 == lon2))
%    o_distance = 0;
% else
%    while (lon1 < 0)
%       lon1 = lon1 + 360;
%    end
%    while (lon2 < 0)
%       lon2 = lon2 + 360;
%    end
%    lat1_r = lat1 / 180.0 * pi;
%    lon1_r = lon1 / 180.0 * pi;
%    lat2_r = lat2 / 180.0 * pi;
%    lon2_r = lon2 / 180.0 * pi;
% 
%    x1 = R_EARTH * cos(lat1_r) * cos(lon1_r);
%    y1 = R_EARTH * cos(lat1_r) * sin(lon1_r);
%    z1 = R_EARTH * sin(lat1_r);
% 
%    x2 = R_EARTH * cos(lat2_r) * cos(lon2_r);
%    y2 = R_EARTH * cos(lat2_r) * sin(lon2_r);
%    z2 = R_EARTH * sin(lat2_r);
% 
%    innerproduct = x1 * x2 + y1 * y2 + z1 * z2;
%    norm_1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
%    norm_2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
%    pre_angle = innerproduct / (norm_1 * norm_2);
%    angle = acos(pre_angle);
% 
%    o_distance = R_EARTH * angle * 1000; % km -> m
% end
% 
% return;
% 
% % conversion jour julien -> jour grégorien par JAMSTEC (pour effectuer les
% % comparaisons)
% function [o_juld] = juld_jamstec(a_juld)
%  
% mmonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
% juldate = fix(a_juld);
% nyear = fix(juldate / 365.2425);
% nday = fix(juldate - (nyear * 365.2425));
% nyear = nyear + 1950;
% 
% if (mod(nyear, 4) == 0)
%    mmonth(2) = 29;
% end
% 
% j = 1;
% while (j < 13)
%    if (nday > mmonth(j))
%       nday = nday - mmonth(j);
%    else
%       nmonth = j;
%       break
%    end
%    j = j + 1;
% end
% if (j == 13)
%    nmonth = 12;
% end
% 
% jultime = fix((a_juld - juldate) * 86400);
% nhour = fix(jultime / 3600);
% jultime = jultime - nhour * 3600;
% nmin = fix(jultime / 60);
% nsec = jultime - nmin * 60;
% 
% o_juld = gregorian_2_julian(sprintf('%04d/%02d/%02d %02d:%02d:%02d', ...
%    nyear, nmonth, nday, nhour, nmin, nsec));
% 
% return;
