% cc 29/09/2020 add a function that compute pres_drift_mes from
% pressure at park (take flag into account).
% from JP RANNOU code dep_compute_park_mes.m
function [tabFinalParkPres,tabFinalParkTemp,tabFinalParkPsal,tabFinalParkEtat]= compute_rpp(T,M,idCyc_drift,idCyc_loc)

depPres=T.pres.data(idCyc_drift);
depTemp=T.temp.data(idCyc_drift);
depPsal=T.psal.data(idCyc_drift);
depType=T.measurement_code.data(idCyc_drift);
depPres_qc=T.pres_qc.data(idCyc_drift);
depTemp_qc=T.temp_qc.data(idCyc_drift);
depPsal_qc=T.psal_qc.data(idCyc_drift);
depType=T.measurement_code.data(idCyc_drift);
idPresKo = (depPres_qc==6|depPres_qc==4|isnan(depPres_qc));
idTempKo = (depTemp_qc==6|depTemp_qc==4|isnan(depTemp_qc));
idPsalKo = (depPsal_qc==6|depPsal_qc==4|isnan(depPsal_qc));

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
tabFinalParkTemp=NaN;
tabFinalParkPsal=NaN;
tabFinalParkEtat=0;
% Priorite 1: cas ou l'on dispose de mesures en phase de
% stabilisation et des mesures en derive
% JMA A3 et A6: mesures en stabilisation toutes les 1.5h puis mesures en
% derive toutes les 6h => moyenne ponderee de ces mesures
% Priorite 2: cas ou l'on dispose de mesures en derive
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
        if ((length(tabTempAdjPark) > 0) || (length(tabTempDriftMes) > 0))
            tabFinalParkTemp = sum(tabTempAdjPark)*1.5 + sum(tabTempDriftMes)*6;
            tabFinalParkTemp = tabFinalParkTemp/(length(tabTempAdjPark)*1.5+length(tabTempDriftMes)*6);
        end
        if ((length(tabPsalAdjPark) > 0) || (length(tabPsalDriftMes) > 0))
            tabFinalParkPsal = sum(tabPsalAdjPark)*1.5 + sum(tabPsalDriftMes)*6;
            tabFinalParkPsal = tabFinalParkPsal/(length(tabPsalAdjPark)*1.5+length(tabPsalDriftMes)*6);
        end
        tabFinalParkEtat = g_etatFromDec;
        finalOk = 1;
    end
    
    
    % Priorite 2: cas ou l'on dispose de la moyenne des mesures en
    % derive
    if (finalOk == 0)
        idCycleMeanParkMes = find(depType == g_typeMeanParkMes & idPresKo==0);
        if (~isempty(idCycleMeanParkMes))
            tabFinalParkPres = depPres((idCycleMeanParkMes));
            tabFinalParkTemp = depTemp((idCycleMeanParkMes));
            tabFinalParkPsal = depPsal((idCycleMeanParkMes));
            tabFinalParkEtat = g_etatFromDec;
            finalOk = 1;
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
            tabFinalParkEtat = g_etatFromTxt;
            finalOk = 1;
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
            tabFinalParkEtat = g_etatFromDec;
            finalOk = 1;
        end
    end
    if (finalOk == 0)
        idCycleParkMes = find(depType == g_typeParkMes2 & idPresKo==0);
        if (~isempty(idCycleParkMes))
            tabFinalParkPres = depPres((idCycleParkMes));
            tabFinalParkTemp = depTemp((idCycleParkMes));
            tabFinalParkPsal = depPsal((idCycleParkMes));
            tabFinalParkEtat = g_etatFromDec;
            finalOk = 1;
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
            
            if ((presMin ~= g_presDef) && (presMax ~= g_presDef))
                tabFinalParkPres = mean([presMin presMax]);
                tabFinalParkTemp = NaN;
                tabFinalParkPsal = NaN;
                tabFinalParkEtat = g_etatFromDec;
                finalOk = 1;
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
        % il n'y a pas de mesure en derive dispo
        if (~isempty(idCyc_loc))&isempty(idCyc_drift)
            %keyboard
            theidMis=find(M.config_mission_number==T.config_mission_number.data(T.cycle_number_index.data==unique(T.cycle_number.data(idCyc_drift))));
            if isempty(theidMis)
                tabFinalParkPres=M.ParkPressure(theidMis);
            end
            
            
            tabFinalParkTemp = NaN;
            tabFinalParkPsal = NaN;
            tabFinalParkEtat = g_etatFromNcMeta;
        end
    end
end
