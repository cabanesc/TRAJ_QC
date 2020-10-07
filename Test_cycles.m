

function [o_alerte3, o_alerte4, o_alerte5, o_alerte6, o_alertCyc_e4, idCycprec] = ...
    Test_cycles(a_idLoc, a_cycles, a_cycles_1, a_cycles_sorted, a_duree_cycle,a_mission, a_dureeMedianCycle)


%%------------------------------------------------------------------------------------------------------------
% BOUCLE SUR CHAQUE CYCLE POUR VERIFIER:
% Vérifie: - les double de cycles
%          - pb de numéros de cycle non croissants
%          - mauvaise duree du cycle
%%INPUT:
%% a_idLoc: indice des locs de surface
%% a_cycles: tableau des numéros de cycles
%% a_cycles_1: tableau des cycles théoriques
%% a_cycles_sorted: tableau des cycles tries
%% a_duree_cycle: durée des cycles
%% a_mission: numéro de mission donné
%% a_dureeMedianCycle: mediane de la duree de cycle calculee pour le numero
%% de mission donné.
%% k: indice du flotteur
%%OUTPUT:
%%% o_alerte3,o_alerte4,o_alerte5,o_alerte6: recuperation des flotteurs avec alertes sur les cycles
%%% o_alertCyc_e4: recuperations des cycles concernes par l'alerte 'cycle'.
%%%idCycprec: indices du cycle n-1
%%---------------------------------------------------------------------------------------------------------------

global PARAM;
global T;
global M;
global file_alerte;
global floatname;

cyc=(1);
idCyc=(1);
idCycprec{1}=[];
o_alertCyc_e4 = [];
o_alerte3=[]; o_alerte4=[]; o_alerte5=[];o_alerte6=[];

for id=1:length(a_cycles)
    %idCyc_prec=idCyc;
    numCycle = a_cycles(id);
    isCyc = (T.cycle_number.data(a_idLoc) == numCycle);
    idCyc = a_idLoc(isCyc);
    numMis = T.config_mission_number.data(id);
    isMis = find(M.config_mission_number == numMis);
    if(id>1)
        isCycprec = (T.cycle_number.data(a_idLoc)==(numCycle-1));
        idCycprec{id} = [a_idLoc(isCycprec)];
    end
    
    
    if (id>1 & sum(idCyc)>0 & sum(idCycprec{id})>0)
        
        cyc(id,1:length(T.juld.data(idCyc)))=T.juld.data(idCyc);
        
        
        for i = 1:size(cyc,1)
            a = cyc(i,:);
            a(find(ismember(cyc(i,:),0)==1))=NaN;    %%%remplace les 0 en Nan pour utilisation ismember.
            cyc(i,:) = a;
        end
        
        if(isempty(find(ismember(cyc(1:id-1,:),cyc(id,:))==1))==0 & id>1)
            %if ismember(cyc(1:id-1,:),cyc(id,:),'rows')>1 & id >1 % NC
            o_alerte3(a_cycles(id)+1)=str2double(floatname);
            fid_alerte=fopen(file_alerte,'a');
            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(a_cycles(id)) ', ? PB DOUBLE CYCLE ?']);
            fclose(fid_alerte);
            fprintf('%s\n',[ floatname ', cycle ' num2str(a_cycles(id)) ', ? PB DOUBLE CYCLE ?']);    %%parfois dû à erreur en fin de batterie, float qui reste en surface
            o_alertCyc_e4 = [o_alertCyc_e4 a_cycles_sorted(id)];
            %T.juld_qc.data(idCyc)=4;
            %T.position_qc.data(idCyc)=4;
            
        end
        
        
        if isequal(a_cycles_sorted,a_cycles)==0
            
            o_alerte4(a_cycles(id)+1)=str2double(floatname);
            o_alertCyc_e4 = [o_alertCyc_e4 a_cycles_sorted(id)];
            fid_alerte=fopen(file_alerte,'a');
            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(a_cycles_sorted(id)) ', CYCLES NON CROISSANTS']);
            fclose(fid_alerte);
            fprintf('%s\n',[ floatname ', cycle ' num2str(a_cycles_sorted(id)) ', CYCLES NON CROISSANTS']);
        end
        
        
        idMis = find(double(a_mission)==numMis);
        % on ne va considerer que les cycles superieurs au 1 (possibilite de prelude qui fausse)
        
        if(length(a_mission)<=length(M.CycleTime) & a_cycles(id)~=0 & ~isnan(numMis))    %% cycles(id)~=0 sinon numMis pouvait etre = NaN et erreur sur conditions suivantes. conditions sur numMis car peut avoir numMis=Nan et cycles(1)=1 quand cycle 0 dans T.cycle_number mais pas de measurement_code 703.
            %%date de la derniere loc du cycle donne - date de la
            %%derniere loc du cycle precedent
            %if id > 2 & ((T.juld.data(idCyc(end))-T.juld.data(idCyc_prec(end)))<=dureeMedianCycle(idMis)*(cycles_sorted(id)-cycles_sorted(id-1))-1 ...
            %     |(T.juld.data(idCyc(end))-T.juld.data(idCyc_prec(end)))>=dureeMedianCycle(idMis)*(cycles_sorted(id)-cycles_sorted(id-1))+1)
            Cyct= T.juld.data(idCyc(end))-T.juld.data(idCycprec{id}(end));
            
            %            if id > 2 & (Cyct<=a_dureeMedianCycle(idMis)*(a_cycles_sorted(id)-a_cycles_sorted(id-1))-1 ...
            %              | Cyct>=a_dureeMedianCycle(idMis)*(a_cycles_sorted(id)-a_cycles_sorted(id-1))+1) ...
            %              | Cyct<M.CycleTime(numMis)-PARAM.TIME_DUREE_CYCLE_M | Cyct>M.CycleTime(numMis)+PARAM.TIME_DUREE_CYCLE_M
            %keyboard
            % correction  cc 15/09/2020 : PARAM.TIME_DUREE_CYCLE_M  : faire un pourcentage de la duree de cycle?
            if id > 2 & (Cyct<=a_dureeMedianCycle(idMis)*(a_cycles_sorted(id)-a_cycles_sorted(id-1))-1 ...
                    | Cyct>=a_dureeMedianCycle(idMis)*(a_cycles_sorted(id)-a_cycles_sorted(id-1))+1) ...
                    | abs(Cyct-M.CycleTime(numMis))./M.CycleTime(numMis)*100>PARAM.TIME_DUREE_CYCLE_M
                
                o_alerte5(a_cycles(id)+1)=str2double(floatname);    %%%ne suffit pas quand il n'y a qu'un seul cycle.
                
                fid_alerte=fopen(file_alerte,'a');
                fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(a_cycles_sorted(id)) ', PB DUREE CYCLE, Date derniere loc: ' datestr(T.juld.data(idCyc(end))+datenum('01011950','ddmmyyyy')) ', Date derniere loc précédente:' datestr(T.juld.data(idCycprec{id}(end))+datenum('01011950','ddmmyyyy')) '. Avec M.CycleTime=' num2str(M.CycleTime(numMis)) 'j']);
                fclose(fid_alerte);
                fprintf('%s\n',[ floatname ', cycle ' num2str(a_cycles_sorted(id)) ', PB DUREE CYCLE. Date derniere loc: ' datestr(T.juld.data(idCyc(end))+datenum('01011950','ddmmyyyy')) ', Date derniere loc précédente:' datestr(T.juld.data(idCycprec{id}(end))+datenum('01011950','ddmmyyyy')) '. Avec M.CycleTime=' num2str(M.CycleTime(numMis)) 'j']);
                o_alertCyc_e4 =  [o_alertCyc_e4 a_cycles_sorted(id)];
                
            end
            
            clear Cyct;
            
            dureeCycle=M.CycleTime(isMis);
            
            if abs(a_dureeMedianCycle(isMis)-M.CycleTime(isMis))<0.1 & ((unique(floor(a_duree_cycle)))>2 ...
                    & isempty(a_duree_cycle(a_duree_cycle<0))==0)
                fid_alerte=fopen(file_alerte,'a');
                fprintf(fid_alerte,'%s\n',[ floatname,',, PB DUREE CYCLE PROVIENT PROBABLEMENT D''UN MAUVAIS ARRANGEMENT DES TRANSMISSIONS ARGOS']);
                fclose(fid_alerte);
                fprintf('%s\n',[ floatname, ' PB DUREE CYCLE PROVIENT PROBABLEMENT D''UN MAUVAIS ARRANGEMENT DES TRANSMISSIONS ARGOS']);
            end
            
        elseif(length(M.CycleTime)<isMis)                      %%% Cas où M.CycleTime mal renseigné.
            fid_alerte=fopen(file_alerte,'a');
            fprintf(fid_alerte,'%s\n',[ floatname, ',  cycle ' num2str(a_cycles_sorted(id)) ', NE PEUT PAS COMPARER LA VALEUR MEDIANE DE LA DUREE DU CYCLE AVEC LA DUREE FOURNIE PAR LE FICHIER META CAR INFO MANQUANTE POUR N° DE MISSION ' num2str(m)]);
            fclose(fid_alerte);
            fprintf('%s\n',[ floatname, ' , cycle ' num2str(a_cycles_sorted(id)) ', NE PEUT PAS COMPARER LA VALEUR MEDIANE DE LA DUREE DU CYCLE AVEC LA DUREE FOURNIE PAR LE FICHIER META CAR INFO MANQUANTE POUR N° DE MISSION ' num2str(m)]);
        end
        
    end
end           %%%fin de la boucle sur les cycles





% VERIFICATION DU DECOUPAGE DES CYCLES

if isequal(a_cycles_1,a_cycles_sorted)==0
    cycles_missing=setdiff(a_cycles_1,a_cycles_sorted);
    ipre = 1;
    ipost = 1;
    cycles_pre_missing = [];
    cycles_post_missing=[];
    for ic=1:length(cycles_missing)% on veut retrouver les cycles precedents et suivants le cycle manquant
        %sans qu'ils soient eux-memes
        %manquants pour estimer la duree ecoulee
        if isempty(find(cycles_missing==cycles_missing(ic)-1, 1))
            cycles_pre_missing(ipre)=a_cycles_1(a_cycles_1==cycles_missing(ic)-1);
            ipre=ipre+1;
        end
        if isempty(find(cycles_missing==cycles_missing(ic)+1, 1))
            cycles_post_missing(ipost)=a_cycles_1(a_cycles_1==cycles_missing(ic)+1);
            ipost=ipost+1;
        end
    end
    fid_alerte=fopen(file_alerte,'a');
    fprintf(fid_alerte, '%s\n',[ floatname ', ' num2str(length(cycles_missing)) ', cycles manquent,(', num2str(cycles_missing), ') ? PB DECOUPAGE CYCLE ?']);
    fclose(fid_alerte);
    fprintf('%s\n',[ floatname ', ' num2str(length(cycles_missing)) ' cycles manquent, ? PB DECOUPAGE CYCLE ?']);
    
    % Calcul de la duree entre les deux cycles encadrants le ou les
    % cycles manquants pour determiner si coherent ou probleme de
    % decoupage
    for im=1:length(cycles_pre_missing)-1
        isCycmpre = (T.cycle_number.data(a_idLoc) == cycles_pre_missing(im));    %indices correspondant aux cycles pre missing
        idCycmpre = a_idLoc(isCycmpre);
        isCycmpost = (T.cycle_number.data(a_idLoc) == cycles_post_missing(im));   %%indices correspondant aux cycles post missing
        idCycmpost = a_idLoc(isCycmpost);
        numMismpost = T.config_mission_number.data(a_cycles==cycles_post_missing(im));
        isMismpost = find(a_mission == numMismpost);
        dureeCyclem = a_dureeMedianCycle(isMismpost);
        if ~isnan(dureeCyclem) & (T.juld.data(idCycmpost(end))-T.juld.data(idCycmpre(end))>(cycles_post_missing(im)-cycles_pre_missing(im))*(dureeCyclem+PARAM.TIME_DUREE_CYCLE_JUMP) ...
                |T.juld.data(idCycmpost(end))-T.juld.data(idCycmpre(end))<(cycles_post_missing(im)-cycles_pre_missing(im))*(dureeCyclem-PARAM.TIME_DUREE_CYCLE_JUMP))
            o_alerte6=str2double(floatname);
            fid_alerte=fopen(file_alerte,'a');
            fprintf(fid_alerte,'%s\n',[ floatname ', cycles ' num2str(cycles_pre_missing(im)) ' - ' num2str(cycles_post_missing(im)) ', SAUTES CYCLES']);
            fclose(fid_alerte);
            fprintf('%s\n',[ floatname ', cycles ' num2str(cycles_pre_missing(im)) ' - ' num2str(cycles_post_missing(im)) ', SAUTES CYCLES']);
            o_alertCyc_e4 = [o_alertCyc_e4 a_cycles_sorted(id)];
        end
    end
end


end
