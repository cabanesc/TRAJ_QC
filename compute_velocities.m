% ------------------------------------------------------------------------------
% Addaptation du code dep_export_to_yo.m Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% Pour calcul et stockage des vitesses de surface et en profondeur dans les fichiers netcdf
%
% SYNTAX :
%   compute_velocities((T,alertes_cycle)
%
% INPUT PARAMETERS :
%   T :structure fichier Traj (ajout des variables vitesses)
%   alertes_cycle : numero de cycle en alerte
%
% OUTPUT PARAMETERS :
%
% EXAMPLES :
%
% SEE ALSO : init_valdef, init_valflag, get_config, get_float_info,
%            compute_surf_vel_start_end
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   07/03/2012 - RNU - creation
%   22/09/2015 - BRN - modification (configVar)
%   19/10/2016 - BRN - modification
%                      PSAL sur RPP mis à FillValue (format ANDRO)
% ------------------------------------------------------------------------------
function T=compute_velocities(T,P,alertes_cycle);


% durée maximale de prise ne compte des positions Argos utilisées pour
% déterminer les vitesses en surface (à la remontée et juste avant la plongée)
SURF_VEL_MAX_DURATION = 6; % en heures

% pas d'affichage des informations intermédiaires
g_verbose = 0;
VERBOSE_LOCAL = 0;
g_typeArgosLoc=703;
% offset de convertion des dates juliennes
JUL_CONVERT = 0;  % on reste par rapport au 01/01/1950 dans les fichiers netcdf (dans atlas c'est par rapport au 01/01/2000)

% mise en oeuvre de critères supplémentaires
KOBA_CHECK = 1;
GROUNDED_CHECK = 1;

% intervalle pour lequel les valeurs de pression en parking sont forcées à 0
% intervalle : ]MAX_OFFSET_FOR_PARK_PRESSURE, 0[
MAX_OFFSET_FOR_PARK_PRESSURE = -1.0;

% on associe à chaque classe de localisation, la précision estimée en latitude
% et longitude fournie (en m) par Argos
precision(1) = 1000;
precision(2) = 350;
precision(3) = 150;

% precision affectée à une position GPS
precision(4) = 30;
[yoCycle, id_sorted] = sort(unique(T.cycle_number.data(find(T.cycle_number.data>=0))));
nbCycles = length(yoCycle);

GroundedNum = ones(length(T.grounded.data), 1)*-1;
GroundedNum(find(T.grounded.data == 'B')) = 1;


floatNum = str2num(T.platform_number.data');
floatNumStr = num2str(floatNum);

yoLonLastLocPrev = ones(nbCycles, 1)*NaN;
yoLatLastLocPrev = ones(nbCycles, 1)*NaN;
yoJuldLastLocPrev = ones(nbCycles, 1)*NaN;

yoLonFirstLocCur = ones(nbCycles, 1)*NaN;
yoLatFirstLocCur = ones(nbCycles, 1)*NaN;
yoJuldFirstLocCur = ones(nbCycles, 1)*NaN;

yoLonLastLocCur = ones(nbCycles, 1)*NaN;
yoLatLastLocCur = ones(nbCycles, 1)*NaN;
yoJuldLastLocCur = ones(nbCycles, 1)*NaN;

yoNbLoc = NaN*zeros(nbCycles, 1);
yoWmo = ones(nbCycles, 1)*(floatNum);
yoCycleNum = ones(nbCycles, 1)*-1;
yoProfNum = ones(nbCycles, 1)*NaN;

% consigne le décalage entre les first/last loc et les locs utilisées pour le calcul des vitesses, en raison de la prise en compte des flags
% 0 indique que les locs utilisées correspondent bien aux first/last loc
% +3 indique que c'est la 4eme loc qui est utilisé à la remontée
% -1 indique que c'est l'avant dernière loc qui est utilisé à la descente 
first_good_loc = ones(nbCycles, 1)*NaN;
last_good_loc =ones(nbCycles, 1)*NaN;


cyNumPrev = -1;
juldSurfPrev = NaN;
longSurfPrev = NaN;
latSurfPrev = NaN;



logfile=[P.DIR_HOME '/logs/compute_velocities.log'];

flog=fopen(logfile,'a');
% récupération des localisations Argos et détermination des vitesses de
% surface

for idCy_sorted = 1:length(yoCycle)
    idCy = id_sorted(idCy_sorted);
    numCycle = yoCycle(idCy_sorted);
    g_cycleNumber = numCycle;
    idCycle = find(T.cycle_number.data == numCycle);
    
    % on ne prend pas en compte les cycles #0 des flotteurs ni les cycles avec Alerte
    
    if (numCycle == 0)| (ismember(numCycle,alertes_cycle))
        
		continue;
    end
    % % contrôle de l'information GROUNDED
    % if (GROUNDED_CHECK == 1)
        % if GroundedNum==1
            % continue;
        % end
    % end
    
    % VALEURS DES LOCALISATIONS ARGOS DE SURFACE
    if isfield(T,'measurement_code')==0
        idLoc =   find(~isnan(T.longitude.data(idCycle))& ~isnan(T.latitude.data(idCycle)));
    else
        idLoc =   find(T.measurement_code.data(idCycle)== g_typeArgosLoc);     %%%ne prend pas la première loca en surface car launch
    end
    
    
    cycleDate = T.juld.data(idCycle(idLoc));
    cycleDateQc = T.juld_qc.data(idCycle(idLoc));
    cycleLon = T.longitude.data(idCycle(idLoc));
    cycleLat = T.latitude.data(idCycle(idLoc));
    cycleAccuracy = T.position_accuracy.data(idCycle(idLoc));
    cycleQc = T.position_qc.data(idCycle(idLoc));
    cycleQc_koba = T.position_qc_koba.data(idCycle(idLoc));
	
	
    if (KOBA_CHECK == 1)
        idGood = find((cycleQc ~= 6) & (cycleQc ~= 4) & (cycleQc_koba ~= 6) & (cycleDateQc ~= 4)& (cycleDateQc ~= 6));
    else
        idGood = find((cycleQc ~= 6) & (cycleQc ~= 4)  & (cycleDateQc ~= 4)& (cycleDateQc ~= 6));
    end
	
	first_good_loc(idCy_sorted)=idGood(1)-1;
	last_good_loc(idCy_sorted)=idGood(end)-length(cycleDate);
    
    cycleDate = cycleDate(idGood);
    cycleLon = cycleLon(idGood);
    cycleLat = cycleLat(idGood);
    cycleQc = cycleQc(idGood);
    cycleAccuracy = cycleAccuracy(idGood);
    cycleQc_koba = cycleQc_koba(idGood);
	
	
	[cycleDate,idloc_sorted] = sort(cycleDate);
    cycleLon = cycleLon(idloc_sorted);
    cycleLat = cycleLat(idloc_sorted);
    cycleQc = cycleQc(idloc_sorted);
    cycleAccuracy = cycleAccuracy(idloc_sorted);
    cycleQc_koba = cycleQc_koba(idloc_sorted);
	
    
    if (numCycle ~= cyNumPrev + 1)
        cyNumPrev = -1;
        juldSurfPrev = NaN;
        longSurfPrev = NaN;
        latSurfPrev = NaN;
    end
    
    if (~isempty(cycleDate))
        cyNumPrev = numCycle;
        juldSurfPrev = cycleDate(end);
        longSurfPrev = cycleLon(end);
        latSurfPrev = cycleLat(end);
    end
    
    
    
    %if (GroundedNum(idCy) == -1)
        if (~isempty(cycleLon))
             % if numCycle==252
		   % keyboard
		   % end
            % estimation de la vitesse de surface
            [T.juld_first_surface_velocity.data(idCy), T.longitude_first_surface_velocity.data(idCy), T.latitude_first_surface_velocity.data(idCy), ...
                T.u_first_surface_velocity.data(idCy), T.v_first_surface_velocity.data(idCy), ...
                T.uerr_first_surface_velocity.data(idCy), T.verr_first_surface_velocity.data(idCy), ...
                T.juld_last_surface_velocity.data(idCy), T.longitude_last_surface_velocity.data(idCy), T.latitude_last_surface_velocity.data(idCy), ...
                T.u_last_surface_velocity.data(idCy), T.v_last_surface_velocity.data(idCy), ...
                T.uerr_last_surface_velocity.data(idCy), T.verr_last_surface_velocity.data(idCy)] = ...
                compute_surf_vel_start_end(cycleDate-JUL_CONVERT, cycleLon, cycleLat, SURF_VEL_MAX_DURATION);
            
            % test sur les positions agglomérées
            % on estime que, pour une trajectoire Argos de surface de plus de 3h
            % dont toutes les localisations peuvent être contenues dans un
            % cercle de R m de rayon centré sur le barycentre des
            % localisations, la vitesse de surface peut être considérée comme
            % nulle et l'erreur associée fournie par les précisions Argos
            % relatives aux première et dernière localisations
            % pour les flotteurs Argos on choisit, selon la classe de loc des
            % première et dernière position de surface:
            % - classe 1 vs classe 1: R = 500 m
            % - classe 1 vs classe 2 ou 3: R = 375 m
            % - classe 2 ou 3 vs classe 2 ou 3: R = 175 m
            % pour les flotteurs Iridium/GPS on choisit R = 50 m
            radius = [500; 375; 175; 50];
            
            if ((cycleDate(end)-cycleDate(1))*24 > 3)
                % on vérifie que les localisations ne sont pas agglomérées
                [testAgloOk, lonBary, latBary, lonCirclePts, latCirclePts] = ...
                    check_argos_positions_aglo2(cycleLon', cycleLat', cycleAccuracy', radius);
                if (testAgloOk == 0)
                    T.u_first_surface_velocity.data(idCy) = 0;
                    T.v_first_surface_velocity.data(idCy) = 0;
                    T.u_last_surface_velocity.data(idCy) = 0;
                    T.v_last_surface_velocity.data(idCy) = 0;
                    
                    cycleAccuracy(find(cycleAccuracy == 'G')) = '4';
                    firstLocErr = precision(str2num(cycleAccuracy(1)))*100;
                    lastLocErr = precision(str2num(cycleAccuracy(end)))*100;
                    surfErr = sqrt(firstLocErr*firstLocErr + lastLocErr*lastLocErr);
                    surfErr = surfErr/((cycleDate(end)-cycleDate(1))*86400);
                    
                    T.uerr_first_surface_velocity.data(idCy) = surfErr;
                    T.verr_first_surface_velocity.data(idCy) = surfErr;
                    T.uerr_last_surface_velocity.data(idCy) = surfErr;
                    T.verr_last_surface_velocity.data(idCy) = surfErr;
                    
                    fprintf(flog,'%d #%d: positions AGGLOMEREES, vitesse de surface forcée à 0, SurfStartUErr=SurfStartVErr=SurfEndUErr=SurfEndVErr=%.1f cm/s\n', ...
                        floatNum, numCycle, surfErr);
                end
            end
            
            % contrôle de vitesses aberrantes (module de la vitesse > 3 m/s)
            if (~isnan(T.u_first_surface_velocity.data(idCy) )) && (~isnan(T.v_first_surface_velocity.data(idCy) ))
                surfStartVel = sqrt(T.u_first_surface_velocity.data(idCy).^2 + T.v_first_surface_velocity.data(idCy).^2);
				%[numCycle,surfStartVel]
				
                if (surfStartVel > 300)
                    T.juld_first_surface_velocity.data(idCy) = NaN;
                    T.longitude_first_surface_velocity.data(idCy) = NaN;
                    T.latitude_first_surface_velocity.data(idCy) = NaN;
                    T.u_first_surface_velocity.data(idCy) = NaN;
                    T.v_first_surface_velocity.data(idCy) = NaN;
                    T.uerr_first_surface_velocity.data(idCy) = NaN;
                    T.verr_first_surface_velocity.data(idCy) = NaN;
                    
                    fprintf(flog,'%d #%d: vitesse de surface à la remontée ABERRANTE ||v||=%.1f > 3 m.s (forcée à la valeur par défaut)\n', ...
                        floatNum, numCycle, surfStartVel);
						%keyboard
                end
            end
            if (~isnan(T.u_last_surface_velocity.data(idCy) )) && (~isnan(T.v_last_surface_velocity.data(idCy) ))
                surfEndVel = sqrt(T.u_last_surface_velocity.data(idCy).^2 + T.v_last_surface_velocity.data(idCy).^2);
                if (surfEndVel > 300)
                    T.juld_last_surface_velocity.data(idCy) = NaN;
                    T.longitude_last_surface_velocity.data(idCy) = NaN;
                    T.latitude_last_surface_velocity.data(idCy) = NaN;
                    T.u_last_surface_velocity.data(idCy) = NaN;
                    T.v_last_surface_velocity.data(idCy) = NaN;
                    T.uerr_last_surface_velocity.data(idCy) = NaN;
                    T.verr_last_surface_velocity.data(idCy) = NaN;
                    
                    fprintf(flog,'%d #%d: vitesse de surface avant la descente ABERRANTE ||v||=%.1f > 3 m.s (forcée à la valeur par défaut)\n', ...
                        floatNum, numCycle, surfEndVel);
                end
            end
            
            % contrôle des valeurs aberrantes (hors format)
            if (abs(T.u_first_surface_velocity.data(idCy)) >  9999.99) || ...
                    (abs(T.v_first_surface_velocity.data(idCy)) >  9999.99) 
                
                fprintf(flog,'%d #%d: valeur (USurfStart, VSurfStart)=(%f, %f) HORS FORMAT (forcée à la valeur par défaut)\n', ...
                    floatNum, numCycle, T.u_first_surface_velocity.data(idCy), T.v_first_surface_velocity.data(idCy));
                
                    T.juld_first_surface_velocity.data(idCy) = NaN;
                    T.longitude_first_surface_velocity.data(idCy) = NaN;
                    T.latitude_first_surface_velocity.data(idCy) = NaN;
                    T.u_first_surface_velocity.data(idCy) = NaN;
                    T.v_first_surface_velocity.data(idCy) = NaN;
                    T.uerr_first_surface_velocity.data(idCy) = NaN;
                    T.verr_first_surface_velocity.data(idCy) = NaN;
            end
            if (abs(T.u_last_surface_velocity.data(idCy)) >  9999.99) || ...
                    (abs(T.v_last_surface_velocity.data(idCy)) >  9999.99) 
                
                fprintf(flog,'%d #%d: valeur (USurfEnd, VSurfEnd)=(%f, %f) HORS FORMAT (forcée à la valeur par défaut)\n', ...
                    floatNum, numCycle,  T.u_last_surface_velocity.data(idCy), T.v_last_surface_velocity.data(idCy));
                
                 T.juld_last_surface_velocity.data(idCy) = NaN;
                    T.longitude_last_surface_velocity.data(idCy) = NaN;
                    T.latitude_last_surface_velocity.data(idCy) = NaN;
                    T.u_last_surface_velocity.data(idCy) = NaN;
                    T.v_last_surface_velocity.data(idCy) = NaN;
                    T.uerr_last_surface_velocity.data(idCy) = NaN;
                    T.verr_last_surface_velocity.data(idCy) = NaN;
            end
            
            % contrôle des valeurs aberrantes (hors format)
            if (abs(T.uerr_first_surface_velocity.data(idCy)) >  9999.99) || ...
                    (abs(T.verr_first_surface_velocity.data(idCy)) >  9999.99) 
                
                fprintf(flog,'%d #%d: valeur (USurfStartErr, VSurfStartErr)=(%f, %f) HORS FORMAT (forcée à la valeur par défaut)\n', ...
                    floatNum, numCycle, T.uerr_first_surface_velocity.data(idCy), T.verr_first_surface_velocity.data(idCy));
                
                T.uerr_first_surface_velocity.data(idCy) = NaN;
                T.verr_first_surface_velocity.data(idCy) = NaN;
            end
            if (abs(T.uerr_last_surface_velocity.data(idCy)) >  9999.99) || ...
                    (abs(T.verr_last_surface_velocity.data(idCy)) >  9999.99) 
                
                fprintf(flog,'%d #%d: valeur (USurfEndErr, VSurfEndErr)=(%f, %f) HORS FORMAT (forcée à la valeur par défaut)\n', ...
                    floatNum, numCycle, T.uerr_last_surface_velocity.data(idCy), T.verr_last_surface_velocity.data(idCy));
                
                T.uerr_last_surface_velocity.data(idCy) = NaN;
                T.verr_last_surface_velocity.data(idCy) = NaN;
            end
            
            % première loc du cycle courant
            yoLonFirstLocCur(idCy) = cycleLon(1);
            yoLatFirstLocCur(idCy) = cycleLat(1);
            yoJuldFirstLocCur(idCy) = cycleDate(1)-JUL_CONVERT;
            
            % dernière loc du cycle courant
            yoLonLastLocCur(idCy) = cycleLon(end);
            yoLatLastLocCur(idCy) = cycleLat(end);
            yoJuldLastLocCur(idCy) = cycleDate(end)-JUL_CONVERT;
            
            % nombre de (bonnes) localisations Argos (ayant passé avec succès le
            % critère de Kobayashi)
            yoNbLoc(idCy) = length(idGood) ;
			
		    diffTimeLoc = mean(diff(cycleDate));
			if isempty(diffTimeLoc)
			diffTimeLoc=-999.99
			end

        end
        
        % numéro de cycle
        yoCycleNum(idCy) = numCycle;
            
        %yoParkPres(idCy) = parkPres;
        %yoParkTemp(idCy) = depTemp(idCycle(idMeanParkMesCycle));
        %yoParkSal(idCy) = depSal(idCycle(idMeanParkMesCycle));
        
        
        % numéro du profil NetCDF
        %yoProfNum(idCy) = unique(depProfNum(idCycle));
    %end
    
    if (~isempty(cycleLon))
        % même lorsque le cycle est GROUNDED, il faut stocker la dernière
        % localisation Argos
        idNextCycle = find(yoCycle == numCycle+1);
        if (~isempty(idNextCycle))
            yoLonLastLocPrev(idNextCycle) = cycleLon(end);
            yoLatLastLocPrev(idNextCycle) = cycleLat(end);
            yoJuldLastLocPrev(idNextCycle) = cycleDate(end)-JUL_CONVERT;
        end
    end
    
end



% % changement de valeur par défaut entre DEP et ANDRO
% yoParkPres(find(yoParkPres == g_presDef)) = g_yoPresDef;
% yoParkTemp(find(yoParkTemp == g_tempDef)) = g_yoTempDef;
% %yoParkSal(find(yoParkSal == g_salDef)) = g_yoSalDef;
% yoParkSal(1:end) = g_yoSalDef;

% yoProfNum(find(yoProfNum == g_profNumDef)) = g_yoProfNumDef;

% estimation des vitesses en profondeur
for idCy_sorted = 1:length(yoCycle)
    idCy = id_sorted(idCy_sorted);
    if (~isnan(yoLonLastLocPrev(idCy)) && ...
		   ~isnan(yoLonFirstLocCur(idCy)) && ...
		   (GroundedNum(idCy) == -1))
        
        % gestion du passage de la ligne de changement de date
        yoLonStart = yoLonLastLocPrev(idCy);
        yoLonEnd = yoLonFirstLocCur(idCy);
        if ((abs(yoLonStart-yoLonEnd) > 180))
            if (yoLonStart < 0)
                yoLonStart = yoLonStart + 360;
            else
                yoLonEnd = yoLonEnd + 360;
            end
        end
        
        T.longitude_deep_velocity.data(idCy) = yoLonStart + (yoLonEnd-yoLonStart)/2;
        if (T.longitude_deep_velocity.data(idCy) >= 180)
            T.longitude_deep_velocity.data(idCy) = T.longitude_deep_velocity.data(idCy) - 360;
        end
        T.latitude_deep_velocity.data(idCy) = yoLatLastLocPrev(idCy) + (yoLatFirstLocCur(idCy)-yoLatLastLocPrev(idCy))/2;
        T.juld_deep_velocity.data(idCy) = yoJuldLastLocPrev(idCy) + (yoJuldFirstLocCur(idCy)-yoJuldLastLocPrev(idCy))/2;
        
        % composante en U de la vitesse en profondeur
        locLat(1) = T.latitude_deep_velocity.data(idCy);
        locLat(2) = locLat(1);
        locLon(1) = yoLonStart;
        locLon(2) = yoLonEnd;
        rangeLon = dist(locLat, locLon);
        T.u_deep_velocity.data(idCy) = (rangeLon*100)/((yoJuldFirstLocCur(idCy)-yoJuldLastLocPrev(idCy))*24*3600);
        if (yoLonEnd < yoLonStart)
            T.u_deep_velocity.data(idCy) = T.u_deep_velocity.data(idCy)*-1;
        end
        
        % estimation de l'erreur sur U
        numCycle = yoCycle(idCy_sorted);
        idPrevCycle = find(yoCycle == numCycle-1);
        if (~isempty(idPrevCycle))
            if (~isnan( T.u_last_surface_velocity.data(idPrevCycle)) &&~isnan(T.u_first_surface_velocity.data(idCy) )&& ...
                    ~isnan(T.u_deep_velocity.data(idCy))  && ~isnan(T.representative_park_pressure.data(idCy))  && ...
                    ~isnan(yoJuldLastLocPrev(idCy) ) && ~isnan(yoJuldFirstLocCur(idCy)))
                [T.uerr_deep_velocity.data(idCy)] = ...
                    compute_deep_vel_err(T.u_last_surface_velocity.data(idPrevCycle), T.u_first_surface_velocity.data(idCy), ...
                    T.u_deep_velocity.data(idCy), T.representative_park_pressure.data(idCy), ...
                    yoJuldLastLocPrev(idCy), yoJuldFirstLocCur(idCy));
            end
        end
        
        % composante en V de la vitesse en profondeur
        locLat(1) = yoLatLastLocPrev(idCy);
        locLat(2) = yoLatFirstLocCur(idCy);
        locLon(1) = yoLonStart;
        locLon(2) = locLon(1);
        rangeLat = dist(locLat, locLon);
        T.v_deep_velocity.data(idCy) = (rangeLat*100)/((yoJuldFirstLocCur(idCy)-yoJuldLastLocPrev(idCy))*24*3600);
        if (yoLatFirstLocCur(idCy) < yoLatLastLocPrev(idCy))
            T.v_deep_velocity.data(idCy) = T.v_deep_velocity.data(idCy)*-1;
        end
        
        % estimation de l'erreur sur V
        if (~isempty(idPrevCycle))
            if (~isnan(T.v_last_surface_velocity.data(idPrevCycle)) && ~isnan(T.v_first_surface_velocity.data(idCy))  && ...
                    ~isnan(T.v_deep_velocity.data(idCy))  && ~isnan(T.representative_park_pressure.data(idCy)) && ...
                    ~isnan(yoJuldLastLocPrev(idCy))  && ~isnan(yoJuldFirstLocCur(idCy)) )
                [T.verr_deep_velocity.data(idCy)] = ...
                    compute_deep_vel_err(T.v_last_surface_velocity.data(idPrevCycle), T.v_first_surface_velocity.data(idCy), ...
                    T.v_deep_velocity.data(idCy), T.representative_park_pressure.data(idCy), ...
                    yoJuldLastLocPrev(idCy), yoJuldFirstLocCur(idCy));
            end
        end
        
        % contrôle des valeurs aberrantes
        if ~isnan(T.uerr_deep_velocity.data(idCy)) && ~isnan(T.verr_deep_velocity.data(idCy))
            if (abs(T.uerr_deep_velocity.data(idCy))  > 9999.99|| ...
                   abs(T.verr_deep_velocity.data(idCy)) > 9999.99)
                
                fprintf(flog,'%d #%d: valeur (uerr_deep_velocity, verr_deep_velocity)=(%f, %f) aberrante (forcée à la valeur par défaut)\n', ...
                    floatNum, yoCycle(idCy), T.uerr_deep_velocity.data(idCy), T.verr_deep_velocity.data(idCy));
                
                T.uerr_deep_velocity.data(idCy) = NaN;
                T.verr_deep_velocity.data(idCy) = NaN;
            end
        end
		
		% consigne les vitesses qui sont calculés avec des positions éloignées des first/last loc
		if first_good_loc(idCy)~=0
		fprintf(flog,'%d #%d: Loc utilisee a la remontee ne correspond pas a la premiere loc argos. First Loc+%d Vitesse de surface: %f Temps moyen(h) entre deux locs  %f \n', ...
                    floatNum, yoCycle(idCy), first_good_loc(idCy), T.v_first_surface_velocity.data(idCy), diffTimeLoc*24);
	    end
		if last_good_loc(idCy)~=0
		fprintf(flog,'%d #%d: Loc utilisee a la descente ne correspond pas a la derniere loc argos. Last Loc%d Vitesse de surface: %f Temps moyen(h) entre deux locs  %f \n', ...
                    floatNum, yoCycle(idCy), last_good_loc(idCy), T.v_last_surface_velocity.data(idCy), diffTimeLoc*24);
	    end
		
    end

end

fclose(flog)

return;

% ------------------------------------------------------------------------------
% Détermination de la date, la position et les composantes U et V des vitesses
% de surface ainsi que de l'erreur associée.
% Deux vitesses de surface sont déterminées:
%    - la première, à la remontée, est déterminée avec les localisations Argos
%      disponibles jusqu'à a_surfMaxDuration heures après la première.
%    - la seconde, avant la plongée, est déterminée avec les localisations Argos
%      disponibles jusqu'à a_surfMaxDuration heures avant la dernière.
%
% SYNTAX :
%  [o_dateSurfStart, o_lonSurfStart, o_latSurfStart, ...
%    o_surfStartU, o_surfStartV, o_surfStartUErr, o_surfStartVErr, ...
%    o_dateSurfEnd, o_lonSurfEnd, o_latSurfEnd, ...
%    o_surfEndU, o_surfEndV, o_surfEndUErr, o_surfEndVErr] = ...
%    compute_surf_vel_start_end(a_locDate, a_locLon, a_locLat, a_surfMaxDuration)
%
% INPUT PARAMETERS :
%   a_locDate         : date des localisations Argos
%   a_locLon          : longitude des localisations Argos
%   a_locLat          : latitude des localisations Argos
%   a_surfMaxDuration : durée maximale de prise ne compte des positions Argos
%                       utilisées pour déterminer les vitesses en surface
%
% OUTPUT PARAMETERS :
%   o_dateSurfStart : date de la vitesse de surface à la remontée
%   o_lonSurfStart  : longitude de la vitesse de surface à la remontée
%   o_latSurfStart  : latitude de la vitesse de surface à la remontée
%   o_surfStartU    : composante U de la vitesse de surface à la remontée
%   o_surfStartV    : composante V de la vitesse de surface à la remontée
%   o_surfStartUErr : erreur sur la composante U de la vitesse de surface à la
%                     remontée
%   o_surfStartVErr : erreur sur la composante V de la vitesse de surface à la
%                     remontée
%   o_dateSurfEnd   : date de la vitesse de surface avant la plongée
%   o_lonSurfEnd    : longitude de la vitesse de surface avant la plongée
%   o_latSurfEnd    : latitude de la vitesse de surface avant la plongée
%   o_surfEndU      : composante U de la vitesse de surface avant la plongée
%   o_surfEndV      : composante V de la vitesse de surface avant la plongée
%   o_surfEndUErr   : erreur sur la composante U de la vitesse de surface avant
%                     la plongée
%   o_surfEndVErr   : erreur sur la composante V de la vitesse de surface avant
%                     la plongée
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   07/03/2012 - RNU - creation
% ------------------------------------------------------------------------------
function [o_dateSurfStart, o_lonSurfStart, o_latSurfStart, ...
    o_surfStartU, o_surfStartV, o_surfStartUErr, o_surfStartVErr, ...
    o_dateSurfEnd, o_lonSurfEnd, o_latSurfEnd, ...
    o_surfEndU, o_surfEndV, o_surfEndUErr, o_surfEndVErr] = ...
    compute_surf_vel_start_end(a_locDate, a_locLon, a_locLat, a_surfMaxDuration)

% vitesse de surface à la remontée
idLocStart = find((a_locDate - a_locDate(1)) <= a_surfMaxDuration/24);

[o_dateSurfStart, o_lonSurfStart, o_latSurfStart, ...
    o_surfStartU, o_surfStartV, o_surfStartUErr, o_surfStartVErr] = ...
    compute_surf_vel(a_locDate(idLocStart), a_locLon(idLocStart), a_locLat(idLocStart));

% vitesse de surface avant la plongée
idLocEnd = find((a_locDate(end) - a_locDate) <= a_surfMaxDuration/24);

[o_dateSurfEnd, o_lonSurfEnd, o_latSurfEnd, ...
    o_surfEndU, o_surfEndV, o_surfEndUErr, o_surfEndVErr] = ...
    compute_surf_vel(a_locDate(idLocEnd), a_locLon(idLocEnd), a_locLat(idLocEnd));

return;

% ------------------------------------------------------------------------------
% Détermination de la date, la position et les composantes d'une vitesse de
% surface ainsi que de l'erreur associée.
%
% SYNTAX :
%   [o_dateSurf, o_lonSurf, o_latSurf, ...
%      o_surfU, o_surfV, o_surfUErr, o_surfVErr] = compute_surf_vel(a_locDate, a_locLon, a_locLat)
%
% INPUT PARAMETERS :
%   a_locDate : date des localisations Argos
%   a_locLon  : longitude des localisations Argos
%   a_locLat  : latitude des localisations Argos
%
% OUTPUT PARAMETERS :
%   o_dateSurf : date de la vitesse de surface
%   o_lonSurf  : longitude de la vitesse de surface
%   o_latSurf  : latitude de la vitesse de surface
%   o_surfU    : composante U de la vitesse de surface
%   o_surfV    : composante V de la vitesse de surface
%   o_surfUErr : erreur sur la composante U de la vitesse de surface
%   o_surfVErr : erreur sur la composante V de la vitesse de surface
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   28/11/2008 - RNU - creation
% ------------------------------------------------------------------------------
function [o_dateSurf, o_lonSurf, o_latSurf, ...
    o_surfU, o_surfV, o_surfUErr, o_surfVErr] = compute_surf_vel(a_locDate, a_locLon, a_locLat)

global g_yoLonDef g_yoLatDef g_yoJuldDef g_yoUVDef g_yoDeepUVErrDef;

% initialisation des valeurs par défaut
%init_valdef;

o_dateSurf = NaN;
o_lonSurf = NaN;
o_latSurf = NaN;
o_surfU = NaN;
o_surfV = NaN;
o_surfUErr = NaN;
o_surfVErr = NaN;

% contrôle du passage de la ligne de changement de date
deltaLon = abs(min(a_locLon) - max(a_locLon));
if (deltaLon > 180)
    id = find(a_locLon < 0);
    a_locLon(id) = a_locLon(id) + 360;
    %    fprintf('COUPE LA LIGNE\n');
end

% la date milieu est prise comme date de référence
dateRef = (a_locDate(end)-a_locDate(1))/2;
locDate = a_locDate - dateRef;

% termes communs
nbLoc = length(a_locDate);
dateMoy = mean(locDate);
lonMoy = mean(a_locLon);
latMoy = mean(a_locLat);
s2 = sum(locDate.*locDate);
denum = (s2 - nbLoc*dateMoy*dateMoy);

% il faut au moins deux localisations Argos pour effectuer un calcul de vitesse
if (nbLoc == 1)
    return;
end

% date et localisation de la vitesse de surface
o_dateSurf = mean(a_locDate);
o_lonSurf = lonMoy;
if (o_lonSurf >= 180)
    o_lonSurf = o_lonSurf - 360;
end
o_latSurf = latMoy;

% calcul de la vitesse en U
s1U = sum(a_locLon.*locDate);
vitLon = (s1U - nbLoc*lonMoy*dateMoy)/(s2 - nbLoc*dateMoy*dateMoy);
o_surfU = vitLon*cosd(a_locLat(1))*(60*1852/864);

% calcul de la vitesse en V
s1V = sum(a_locLat.*locDate);
vitLat = (s1V - nbLoc*latMoy*dateMoy)/(s2 - nbLoc*dateMoy*dateMoy);
o_surfV = vitLat*(60*1852/864);

% il faut au moins trois localisations Argos pour effectuer un calcul d'erreur
% sur la vitesse déterminée
if (nbLoc == 2)
    return;
end

% calcul de l'erreur sur la vitesse en U
s2U = sum(a_locLon.*a_locLon);
num = s2U - nbLoc*lonMoy*lonMoy - (vitLon*(s1U - nbLoc*lonMoy*dateMoy));
o_surfUErr = sqrt(num/((nbLoc-2)*denum));
o_surfUErr = o_surfUErr*cosd(a_locLat(1))*(60*1852/864);

% calcul de l'erreur sur la vitesse en V
s2V = sum(a_locLat.*a_locLat);
num = s2V - nbLoc*latMoy*latMoy - (vitLat*(s1V - nbLoc*latMoy*dateMoy));
o_surfVErr = sqrt(num/((nbLoc-2)*denum));
o_surfVErr = o_surfVErr*(60*1852/864);

return;

% ------------------------------------------------------------------------------
% Estimation de l'erreur sur une composante de la vitesse en profondeur
%
% SYNTAX :
%   [o_deepUVErr] = compute_deep_vel_err( ...
%      a_surfUVPrev, a_surfUVCur, a_deepUV, a_parkPres, a_juldLastLocPrev, a_juldFirstLocCur)
%
% INPUT PARAMETERS :
%   a_surfUVPrev      : composante de la vitesse en surface du cyle précédent
%   a_surfUVCur       : composante de la vitesse en surface du cyle courant
%   a_deepUV          : composante de la vitesse en profondeur
%   a_parkPres        : immersion de dérive
%   a_juldLastLocPrev : date de la dernière localisation Argos du cycle
%                       précédent
%   a_juldFirstLocCur : date de la première localisation Argos du cycle
%                       courant
%
% OUTPUT PARAMETERS :
%   o_deepUVErr : erreur sur la composante de la vitesse en profondeur
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   28/11/2008 - RNU - creation
% ------------------------------------------------------------------------------
function [o_deepUVErr] = compute_deep_vel_err( ...
    a_surfUVPrev, a_surfUVCur, a_deepUV, a_parkPres, a_juldLastLocPrev, a_juldFirstLocCur)

global g_yoDeepUVErrDef;

% initialisation des valeurs par défaut
%init_valdef;

o_deepUVErr = NaN;

% on fait l'hypothèse d'une vitesse ascensionnelle de 10 cm/s
ASC_SPEED = 10;

% estimation de l'erreur sur la vitesse en profondeur
velUVSurfMoy = (a_surfUVPrev + a_surfUVCur)/2;
alpha = (a_parkPres*100)/ASC_SPEED/((a_juldFirstLocCur-a_juldLastLocPrev)*86400);

o_deepUVErr = abs((velUVSurfMoy-a_deepUV)*alpha/(1-alpha));

% % pour test, calcul tel que préconisé dans YoMaHa
% o_deepUVErr = abs((a_surfUVCur-a_deepUV)*alpha/(1-alpha));

return;
