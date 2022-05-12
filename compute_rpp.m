% cc 29/09/2020 add a function that compute pres_drift_mes from
% pressure at park (take flag into account).
% from JP RANNOU code dep_compute_park_mes.m
% cc 21/01/2021 ajout de tabFinalMaxParkPres (valeur max de la pression en derive)
function [tabFinalParkPres,tabFinalParkTemp,tabFinalParkPsal,tabFinalParkEtat,tabFinalMaxParkPres]= compute_rpp(T,M,idCyc_drift,idCyc_loc,prof_fileName,verbose)

global file_alerte;
global floatname;
global PARAM

depPres=T.pres.data(idCyc_drift);
depPres_qc=T.pres_qc.data(idCyc_drift);
if isfield(T,'temp')
    depTemp=T.temp.data(idCyc_drift);
    depTemp_qc=T.temp_qc.data(idCyc_drift);
else
    depTemp=NaN*depPres;
    depTemp_qc=NaN*depPres_qc;
end
if isfield(T,'psal')
    depPsal=T.psal.data(idCyc_drift);
    depPsal_qc=T.psal_qc.data(idCyc_drift);
else
    depPsal=NaN*depPres;
    depPsal_qc=NaN*depPres_qc;
end
depType=T.measurement_code.data(idCyc_drift);
depType=T.measurement_code.data(idCyc_drift);
idPresKo = (depPres_qc==6|depPres_qc==4|depPres_qc==3|isnan(depPres_qc));
idTempKo = (depTemp_qc==6|depTemp_qc==4|depTemp_qc==3|isnan(depTemp_qc));
idPsalKo = (depPsal_qc==6|depPsal_qc==4|depPsal_qc==3|isnan(depPsal_qc));
idPresNaN = (isnan(depPres_qc));
idTempNaN = (isnan(depTemp_qc));
idPsalNaN = (isnan(depPsal_qc));

depPres(idPresKo)=NaN;
depTemp(idTempKo)=NaN;
depPsal(idPsalKo)=NaN;

g_typeAdjPark=189;
g_typeDriftMes=290;
g_typeMeanParkMes=296;
g_typeMedianParkMes=295;
g_typeParkMes=299;
g_typeParkMes2=300;
g_typeMinPresAtParkPres=297;
g_typeMaxPresAtParkPres=298;

g_typeArgosLoc=703;

g_etatFromNcMeta = 6;
g_etatFromDec = 4;
g_etatFromTxt = 5;
finalOk = 0;
tabFinalParkPres=NaN;
tabFinalMaxParkPres=NaN;
tabFinalParkTemp=NaN;
tabFinalParkPsal=NaN;
tabFinalParkEtat=' ';
% Priorite 1: cas ou l'on dispose de mesures en phase de
% stabilisation et des mesures en derive
% JMA A3 et A6: mesures en stabilisation toutes les 1.5h puis mesures en
% derive toutes les 6h => moyenne ponderee de ces mesures
% Priorite 2: cas ou l'on dispose de mesures en derive
%keyboard

if (finalOk == 0)
    
    tabPresAdjPark = [];
    tabTempAdjPark = [];
    tabPresDriftMes = [];
    tabTempDriftMes = [];
    tabPsalAdjPark =[];
    tabPsalDriftMes = [];
    
    idCycleAdjPark = find(depType ==g_typeAdjPark );
    if (~isempty(idCycleAdjPark))
        tabPresAdjPark = depPres(idCycleAdjPark);
        tabTempAdjPark = depTemp(idCycleAdjPark);
        tabPsalAdjPark = depPsal(idCycleAdjPark);
        
        tabPresAdjPark(isnan(tabPresAdjPark)) = [];
        tabTempAdjPark(isnan(tabTempAdjPark)) = [];
        tabPsalAdjPark(isnan(tabPsalAdjPark)) = [];
        
    end
    
    idCycleDriftMes = find(depType == g_typeDriftMes);
    if (~isempty(idCycleDriftMes))
        tabPresDriftMes = depPres(idCycleDriftMes);
        tabTempDriftMes = depTemp(idCycleDriftMes);
        tabPsalDriftMes = depPsal(idCycleDriftMes);
        
        tabPresDriftMes(isnan(tabPresDriftMes)) = [];
        tabTempDriftMes(isnan(tabTempDriftMes)) = [];
        tabPsalDriftMes(isnan(tabPsalDriftMes)) = [];
        % cas des APEX IR 1004 de CSIRO: la derni�re mesure en
        % d�rive est parfois r�alis�e � une immersion proche de
        % celle de d�but de profil, on supprime cette mesure de la
        % moyenne
        %                if (length(tabPresDriftMes) > 5)
        %                   idDel = find(abs(tabPresDriftMes-mean(tabPresDriftMes)) > 7*std(tabPresDriftMes));
        %                   if (~isempty(idDel))
        %                      fprintf('INFO : cycle %d non prise en compte mesure %d/%d\n', ...
        %                         numCycle, idDel, length(tabPresDriftMes));
        %                      tabPresDriftMes(idDel) = [];
        %                      tabTempDriftMes(idDel) = [];
        %                   end
        %                elseif (length(tabPresDriftMes) > 2)
        %                   idDel = [];
        %                   for id = 1:length(tabPresDriftMes)
        %                      tabPresTmp = tabPresDriftMes;
        %                      tabPresTmp(id) = [];
        %                      if (abs(tabPresDriftMes(id)-mean(tabPresTmp)) > 7*std(tabPresTmp))
        %                         idDel = [idDel; id];
        %                         fprintf('INFO : cycle %d non prise en compte mesure %d/%d\n', ...
        %                            numCycle, id, length(tabPresDriftMes));
        %                      end
        %                   end
        %                   tabPresDriftMes(idDel) = [];
        %                   tabTempDriftMes(idDel) = [];
        %                end
    end
    
    nbAdjPark = length(tabPresAdjPark);
    nbDriftMes = length(tabPresDriftMes);
    
    if ((nbAdjPark > 0) || (nbDriftMes > 0))
        tabFinalParkPres = sum(tabPresAdjPark)*1.5 + sum(tabPresDriftMes)*6;
        tabFinalParkPres = tabFinalParkPres/(nbAdjPark*1.5+nbDriftMes*6);
        tabFinalMaxParkPres = max([max(tabPresDriftMes),max(tabPresAdjPark)]);
        if ((length(tabTempAdjPark) > 0) || (length(tabTempDriftMes) > 0))
            tabFinalParkTemp = sum(tabTempAdjPark)*1.5 + sum(tabTempDriftMes)*6;
            tabFinalParkTemp = tabFinalParkTemp/(length(tabTempAdjPark)*1.5+length(tabTempDriftMes)*6);
        end
        if ((length(tabPsalAdjPark) > 0) || (length(tabPsalDriftMes) > 0))
            tabFinalParkPsal = sum(tabPsalAdjPark)*1.5 + sum(tabPsalDriftMes)*6;
            tabFinalParkPsal = tabFinalParkPsal/(length(tabPsalAdjPark)*1.5+length(tabPsalDriftMes)*6);
        end
        tabFinalParkEtat = '1';
        finalOk = 1;
        
    end
    
    
    % Priorite 2: cas ou l'on dispose de la moyenne des mesures en
    % derive
    if (finalOk == 0)
        idCycleMeanParkMes = find(depType == g_typeMeanParkMes & idPresKo==0);
        if (~isempty(idCycleMeanParkMes))
            tabFinalParkPres = mean(depPres((idCycleMeanParkMes)));
            tabFinalParkTemp = mean(depTemp((idCycleMeanParkMes)));
            tabFinalParkPsal = mean(depPsal((idCycleMeanParkMes)));
            tabFinalParkEtat = '2';
            finalOk = 1;
            
            idCycleMaxPresAtParkPres = find(depType == g_typeMaxPresAtParkPres & idPresKo==0);
            if (~isempty(idCycleMaxPresAtParkPres))
                tabFinalMaxParkPres =depPres((idCycleMaxPresAtParkPres));
            else
                tabFinalMaxParkPres =tabFinalParkPres;
            end
        end
    end
    
    % Priorite3: cas ou l'on dispose de la mediane des mesures en
    % derive
    if (finalOk == 0)
        idCycleMedianParkMes = find(depType == g_typeMedianParkMes & idPresKo==0);
        if (~isempty(idCycleMedianParkMes))
            tabFinalParkPres = depPres(dCycle(idCycleMedianParkMes));
            tabFinalParkTemp = NaN;
            tabFinalParkPsal = NaN;
            tabFinalParkEtat = '3';
            finalOk = 1;
            idCycleMaxPresAtParkPres = find(depType == g_typeMaxPresAtParkPres & idPresKo==0);
            if (~isempty(idCycleMaxPresAtParkPres))
                tabFinalMaxParkPres =depPres((idCycleMaxPresAtParkPres));
            else
                tabFinalMaxParkPres =tabFinalParkPres;
            end
        end
    end
    
    % Priorite 4: cas ou l'on dispose d'une mesure en fin de phase de
    % derive
    if (finalOk == 0)
        idCycleParkMes = find(depType == g_typeParkMes & idPresKo==0);
        if (~isempty(idCycleParkMes))
            tabFinalParkPres = depPres((idCycleParkMes));
            tabFinalParkTemp = depTemp((idCycleParkMes));
            tabFinalParkPsal = depPsal((idCycleParkMes));
            tabFinalParkEtat = '4';
            finalOk = 1;
            idCycleMaxPresAtParkPres = find(depType == g_typeMaxPresAtParkPres & idPresKo==0);
            if (~isempty(idCycleMaxPresAtParkPres))
                tabFinalMaxParkPres =depPres((idCycleMaxPresAtParkPres));
            else
                tabFinalMaxParkPres =tabFinalParkPres;
            end
        end
        
    end
    if (finalOk == 0)
        idCycleParkMes = find(depType == g_typeParkMes2 & idPresKo==0);
        if (~isempty(idCycleParkMes))
            tabFinalParkPres = depPres((idCycleParkMes));
            tabFinalParkTemp = depTemp((idCycleParkMes));
            tabFinalParkPsal = depPsal((idCycleParkMes));
            tabFinalParkEtat = '4';
            finalOk = 1;
            idCycleMaxPresAtParkPres = find(depType == g_typeMaxPresAtParkPres & idPresKo==0);
            if (~isempty(idCycleMaxPresAtParkPres))
                tabFinalMaxParkPres =depPres((idCycleMaxPresAtParkPres));
            else
                tabFinalMaxParkPres =tabFinalParkPres;
            end
        end
    end
    
    % Priorite 5: cas ou l'on dispose des mesures de pression
    % extremales de la phase de derive
    if (finalOk == 0)
        idCycleMinPresAtParkPres = find(depType == g_typeMinPresAtParkPres & idPresKo==0);
        idCycleMaxPresAtParkPres = find(depType == g_typeMaxPresAtParkPres & idPresKo==0);
        if ((~isempty(idCycleMinPresAtParkPres)) && (~isempty(idCycleMaxPresAtParkPres)))
            presMin = depPres((idCycleMinPresAtParkPres));
            presMax = depPres((idCycleMaxPresAtParkPres));
            
            
            tabFinalParkPres = mean([presMin presMax]);
            tabFinalMaxParkPres =presMax;
            tabFinalParkTemp = NaN;
            tabFinalParkPsal = NaN;
            tabFinalParkEtat = '5';
            finalOk = 1;
            
        else
            if (~isempty(idCycleMaxPresAtParkPres))
                presMax = depPres((idCycleMaxPresAtParkPres));
                tabFinalMaxParkPres =presMax;
            end
        end
    end
    
    %     % Priorite 6: cas d'un cycle BOUNCE
    %     if (finalOk == 0)
    %         idCycleProfAsc = find(depType(idCycle) == g_typeProfAsc);
    %         if (~isempty(idCycleProfAsc))
    %             tabOrdre = depOrdre(idCycle(idCycleProfAsc));
    %             tabPres = depPres(idCycle(idCycleProfAsc));
    %             tabTemp = depTemp(idCycle(idCycleProfAsc));
    %             if (min(tabOrdre) > 10000)
    %                 % il s'agit de profils bounce
    %                 tabPresProfBMax = [];
    %                 tabPresProfBMean = [];
    %                 tabTempProfBMax = [];
    %                 tabTempProfBMean = [];
    %                 for idProfB = 1:g_nbProfInBounceCycle
    %                     idMesProfB = find((tabOrdre > idProfB*10000) & (tabOrdre < (idProfB+1)*10000));
    %                     if (~isempty(idMesProfB))
    %                         tabPresProfBMes = tabPres(idMesProfB);
    %                         tabTempProfBMes = tabTemp(idMesProfB);
    %
    %                         idKo = find((tabPresProfBMes == g_presDef) | (tabTempProfBMes == g_tempDef));
    %                         tabPresProfBMes(idKo) = [];
    %                         tabTempProfBMes(idKo) = [];
    %
    %                         if (~isempty(tabPresProfBMes))
    %                             tabPresProfBMax = [tabPresProfBMax; tabPresProfBMes(1)];
    %                             tabPresProfBMean = [tabPresProfBMean; mean(tabPresProfBMes)];
    %                             tabTempProfBMax = [tabTempProfBMax; tabTempProfBMes(1)];
    %                             tabTempProfBMean = [tabTempProfBMean; mean(tabTempProfBMes)];
    %                         end
    %                     end
    %                 end
    %
    %                 if (~isempty(tabPresProfBMax))
    %                     tabFinalParkPres(idCy) = (sum(tabPresProfBMax) + sum(tabPresProfBMean))/(2*length(tabPresProfBMax));
    %                     tabFinalParkTemp(idCy) = (sum(tabTempProfBMax) + sum(tabTempProfBMean))/(2*length(tabTempProfBMax));
    %                     tabFinalParkEtat(idCy) = g_etatFromDec;
    %                     finalOk = 1;
    %                 end
    %             end
    %         end
    %     end
    
    % Priorite 7: a defaut, on met la valeur de derive fournie par le
    % meta-donnee
    
    if (finalOk == 0)
        % on ne le fait que si le cycle a au moins une position Argos et si
        % il n'y a pas de mesure en derive disponible
        if (~isempty(idCyc_loc))&(isempty(idCyc_drift)||sum(idPresNaN)==length(idPresNaN))
            % on verifie en plus que le flotteur a bien fait un profil  qui a atteint au moins  ParkPressure pour ce cycle
            if exist(prof_fileName)
                
                P = read_netcdf_allthefile(prof_fileName);
                P = replace_fill_bynan(P);
                P = format_flags_char2num(P);
                idcycleProf = find(P.cycle_number.data==unique(T.cycle_number.data(idCyc_loc)));
                pres_prof = P.pres.data(idcycleProf,:);
                pres_prof_qc = P.pres_qc.data(idcycleProf,:);
                isok_pres=find(~isnan(pres_prof)&pres_prof_qc<3);
                if ~isempty (isok_pres)
                    max_pres_prof = max(pres_prof(isok_pres));
                else
                    max_pres_prof =0;
                end
                
                theidMis = find(M.config_mission_number==T.config_mission_number.data(T.cycle_number_index.data==unique(T.cycle_number.data(idCyc_loc))));
                meta_park_pressure = M.ParkPressure(theidMis);
                if isempty(theidMis)||length(theidMis)>1    % add cc 29/04/2022 on perd beaucoup de donnees lorsque T.config_mission_number est n'importe quoi eg csio
                    if length(unique(M.ParkPressure))==1
                        theidMis =1;
                        meta_park_pressure = M.ParkPressure(theidMis);
                        if verbose==1
                            fid_alerte=fopen(file_alerte,'a');
                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(unique(T.cycle_number.data(idCyc_loc))) ',warning, Mission number in traj does not allow to retrieve PARK META, but we used the UNIQUE VALUE of PARK META (same for all missions)']);
                            fclose(fid_alerte);
                            fprintf('%s\n',[ floatname ', cycle ' num2str(unique(T.cycle_number.data(idCyc_loc))) ',warning, Mission number in traj does not allow to retrieve PARK META, but we used the UNIQUE VALUE of PARK META (same for all missions)']);
                        end
                    end
                end
                if length(theidMis)==1
                    if meta_park_pressure<=max_pres_prof+PARAM.PRESS_PARK_DIFF_BATH
                        
                        if verbose==1
                            fid_alerte=fopen(file_alerte,'a');
                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(unique(T.cycle_number.data(idCyc_loc))) ',warning, RPP == META, : ' num2str(meta_park_pressure)]);
                            fclose(fid_alerte);
                            fprintf('%s\n',[ floatname ', cycle ' num2str(unique(T.cycle_number.data(idCyc_loc))) ',warning, RPP == META, : ' num2str(meta_park_pressure)]);
                        end
                        tabFinalParkPres = M.ParkPressure(theidMis);
                        tabFinalMaxParkPres =NaN;
                        tabFinalParkEtat = '6';
                        
                        
                        
                    else
                        if verbose==1&unique(T.cycle_number.data(idCyc_loc))~=0
                            
                            fid_alerte=fopen(file_alerte,'a');
                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(unique(T.cycle_number.data(idCyc_loc))) ',warning, RPP == FILLVALUE, NO PROFILE']);
                            fclose(fid_alerte);
                            fprintf('%s\n',[ floatname ', cycle ' num2str(unique(T.cycle_number.data(idCyc_loc))) ',warning, RPP == FILLVALUE, NO PROFILE']);
                        end
                    end
                end
                
                
            else
                if verbose==1&unique(T.cycle_number.data(idCyc_loc))~=0
                    fid_alerte=fopen(file_alerte,'a');
                    fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(unique(T.cycle_number.data(idCyc_loc))) ',warning, RPP == FILLVALUE, NO pROFILE']);
                    fclose(fid_alerte);
                    fprintf('%s\n',[ floatname ', cycle ' num2str(unique(T.cycle_number.data(idCyc_loc))) ',warning, RPP == FILLVALUE, NO PROFILE']);
                end
            end
            
            
        end
    end
end
