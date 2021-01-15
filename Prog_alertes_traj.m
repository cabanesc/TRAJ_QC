
%%%% Programme visant � qualifier les donn�es trajectoire des flotteurs
%%%% Argo via diff�rents tests, principalement:
%%%  - contr�le des cycles
%%%  - contr�le des localisations argos (dates et positions)
%%%  - contr�le des pressions/profondeur de parking, grounded ou non
%%%  - contr�lle des d�rives en surface et en profondeur
%%%
%%  Date:  04/2020
%%% Auteur: Gaelle Herbert: gaelle.herbert@ecomail.fr
%%  version Matlab: 2016b
%%----------------------------------------------------------------------
clear all
close all

global PARAM;
global P;
global floatname;
global TEMP_LONG;
global TEMP_LAT;
global TEMP_DEPTH;
global TEMP_DEPTH_STD;
global TEMP_TEMP;
global TEMP_TEMP_STD;
global PSAL_DEPTH_STD;
global PSAL_PSAL;
global PSAL_PSAL_STD;
global I_psal_std;
global I_psal;
global I_temp_std;
global I_temp;
global DIR_BATHY % add cc 15/09/2020
%%% fichier config - recuperation des path
P= config;
DIR_FTP = P.DIR_FTP;
DAC = P.DAC;
DIR_ISAS = P.DIR_ISAS;
DIR_HOME = P.DIR_HOME;
DIR_BATHY = P.DIR_BATHY;
DIR_TOOL = P.DIR_TOOL;
DIR_SOFT = P.DIR_SOFT;
DIR_VISU = P.DIR_VISU;
DIR_STAT = P.DIR_STAT;
Wr_Err  = P.Wr_Err;
Bathy = P.Bathy;
map = P.map;
Stat = P.Stat;
Liste_Float = P.Liste_Float


addpath(DIR_TOOL);
addpath([DIR_TOOL '/Lib_Argo/']);
addpath([DIR_TOOL '/Lib_Argo/Plots']);
addpath([DIR_TOOL 'Lib_Argo/Tests_TR']);
addpath([DIR_TOOL '/Lib_Argo/RWnetcdf/R2008b']);

addpath([DIR_SOFT '/decArgo_soft/soft/sub']);
addpath([DIR_SOFT '/decArgo_soft/soft/util']);
addpath([DIR_SOFT '/decArgo_soft/soft/util/sub']);
addpath([DIR_SOFT '/decArgo_soft/soft/util2']);
addpath([DIR_SOFT '/decArgo_soft/soft/misc']);
addpath([DIR_SOFT '/decArgo_soft/soft/sub_foreign']);
addpath([DIR_SOFT '/decArgo_soft/soft/m_map1.4e']);

addpath(DIR_VISU);

% recuperation des parametres
PARAM = param


%for ilist=1:length(Liste_Float)   % add boucle cc 02/11/2020
for ilist=1:1 % cc pour le moment on ne prend que la premiere liste: pas encore teste que ca marche bien d'enchainer les listes!!
    % Liste des flotteurs a tester
    
    fid1 = fopen(Liste_Float{ilist});
    %allfloats=textscan(fid1,'%s','Delimiter',',','CommentStyle','#');
    allfloats=textscan(fid1,'%s%s%s','Delimiter',',','CommentStyle','#');
    fclose(fid1);
    
    %% Creation du répertoire /STAT s'il n'existe pas
    if ~exist([DIR_STAT P.liste_rep_alerte{ilist}],'dir')
        mkdir([DIR_STAT P.liste_rep_alerte{ilist}]);
    end
    %% Creation du répertoire /HOME/alerts/... s'il n'existe pas  % add cc 02/11/2020
    if ~exist([DIR_HOME 'alerts/' P.liste_rep_alerte{ilist}],'dir')
        mkdir([DIR_HOME 'alerts/' P.liste_rep_alerte{ilist}]);
    end
    
    
    %% Fichier log des alertes
    global file_alerte;
    
    file_alerte = [DIR_HOME 'alerts/' P.liste_rep_alerte{ilist} '/RT_alertes.txt'];  % modif chemin alertes cc 02/11/2020
    fid_alerte =fopen(file_alerte,'wt');
    fclose(fid_alerte);
    %% Fichiers log des alertes par erreurs
    if(Bathy==1)
        file3=[DIR_HOME 'alerts/'  P.liste_rep_alerte{ilist} '/Alertes_grounded_etopo.txt'];
        fid3=fopen(file3,'w+');
    else
        file3=[DIR_HOME 'alerts/'  P.liste_rep_alerte{ilist} '/Alertes_grounded_srtm_b.txt'];
        fid3=fopen(file3,'w+');
    end
    file4=[DIR_HOME 'alerts/'  P.liste_rep_alerte{ilist} '/Alertes_cycle_b.txt'];
    file5=[DIR_HOME 'alerts/'  P.liste_rep_alerte{ilist} '/Alertes_locdate_b.txt'];
    file6=[DIR_HOME 'alerts/'  P.liste_rep_alerte{ilist} '/Alertes_launchdate_b.txt'];
    file7=[DIR_HOME 'alerts/'  P.liste_rep_alerte{ilist} '/Alertes_pressure_b.txt'];
    file8=[DIR_HOME 'alerts/'  P.liste_rep_alerte{ilist} '/Alertes_locterre_b.txt'];
    file9=[DIR_HOME 'alerts/'  P.liste_rep_alerte{ilist} '/Alertes_locpos_b.txt'];
    file10=[DIR_HOME 'alerts/'  P.liste_rep_alerte{ilist} '/Alertes_metatrajcyc_b.txt'];
    
    
    fid4 = fopen(file4,'w+');
    fid5= fopen(file5,'w+');
    fid6 = fopen(file6,'w+');
    fid7 = fopen(file7,'w+');
    fid8 = fopen(file8,'w+');
    fid9 = fopen(file9,'w+');
    fid10 = fopen(file10,'w+');
    
    fclose(fid3);
    fclose(fid4);fclose(fid5); fclose(fid6);fclose(fid7); fclose(fid8);fclose(fid9);
    fclose(fid10);
    %%
    istat=0;
    Ncvar.cycle_number.name='CYCLE_NUMBER';
    Ncvar.longitude.name='LONGITUDE';
    Ncvar.latitude.name='LATITUDE';
    Ncvar.juld.name='JULD';
    Ncvar.measurement_code.name='MEASUREMENT_CODE';
    Ncvar.pres='PRES';
    Ncvar.temp='TEMP';
    Ncvar.psal='PSAL';
    icounterfloat = 0;
    
    if(map==1)
        % Affichage d'une carte du monde pour localiser les aberrations
        figure(1)
        zone_visu=[-80 80 -179 180];
        reso='LR';proj='mercator';
        m_proj(proj,'long',[zone_visu(3) zone_visu(4)],'lat',[zone_visu(1) zone_visu(2)]);
        m_grid('box','fancy','tickdir','in');
        hold on
        
        % Trace des contours de bathy à 0 200 1000 et 2000m
        [bathy,lon_bathy,lat_bathy]=m_tbase([zone_visu(3) zone_visu(4) zone_visu(1) zone_visu(2)]);
        v=[-2000 -2000];
        [n,p]=m_contour(lon_bathy,lat_bathy,bathy,v,'color',[.6 .6 .6]);
        v=[-1000 -1000];
        [q,r]=m_contour(lon_bathy,lat_bathy,bathy,v,'color',[.7 .7 .7]);
        v=[-200 -200];
        [c,h]=m_contour(lon_bathy,lat_bathy,bathy,v,'color',[.8 .8 .8]);
        v=[0 0];
        [c,h]=m_contour(lon_bathy,lat_bathy,bathy,v,'k','linewidth',2);
        xlabel('Longitude');
        ylabel('Latitude');
        title('Localisation des alertes de positionnement et de mesures');
    end
    
    %% Lecture des parametres temperature et salinite de la climatologie ISAS13
    
    disp('lecture clim isas')
    isas_temp_std_file=[DIR_ISAS 'ISAS13_CLIM_ann_STD_TEMP.nc'];
    isas_temp_file=[DIR_ISAS 'ISAS13FD_ann_TEMP.nc'];
    I_temp_std=read_netcdf_allthefile(isas_temp_std_file);
    I_temp=read_netcdf_allthefile(isas_temp_file);
    TEMP_LONG=I_temp.longitude.data;
    TEMP_LAT=I_temp.latitude.data;
    TEMP_DEPTH=I_temp.depth.data;
    TEMP_TEMP=I_temp.temp.data;
    TEMP_DEPTH_STD=I_temp_std.depth.data;
    TEMP_TEMP_STD=I_temp_std.temp_std.data;
    isas_psal_std_file=[DIR_ISAS 'ISAS13_CLIM_ann_STD_PSAL.nc'];
    isas_psal_file=[DIR_ISAS 'ISAS13FD_ann_PSAL.nc'];
    I_psal_std=read_netcdf_allthefile(isas_psal_std_file);
    I_psal=read_netcdf_allthefile(isas_psal_file);
    PSAL_PSAL=I_psal.psal.data;
    PSAL_DEPTH_STD=I_psal_std.depth.data;
    PSAL_PSAL_STD=I_psal_std.psal_std.data;
    disp('fin lecture clim isas')
    
    % Initialisation des listes de flotteurs impactes par les differentes alertes
    alerte1=[];alerte2=[];alerte3={};alerte4={};alerte5={};alerte6={};alerte7={};alerte8={};
    alerte9={};alerte10=[];alerte11={};alerte12={};alerte13=[];alerte14=[];alerte15=[];
    alerte16=[];alerte17=[];alerte18=[];alerte19=[];alerte20=[];alerte21=[];alerte22=[];alerte23=[];
    alerte24=[];alerte25=[];alerte26=[];alerte27=[];alerte28=[];alerte29=[];alerte30=[];
    
    % Initialisation des tableaux/cellules/variables
    Num_Cyc_fl = {};
    ifloat = 0;
    ind=0;
    alertes3 = [];
    tab = [1:1:100];
    ite=0;
    pp=0;
    ppp=0;
    pp2=0;
    ecartMaxLocCycle_fl=[];
    ecartMeanLocCycle_fl = [];
    
    global T;
    global M;
    
    %% DEBUT DES TESTS
    %%************************************************************************
    
    %% BOUCLE SUR LES FLOTTEURS
    
    for k=1:length(allfloats{1})
        
        
        alertCyc_e1 = []; alertCyc_e2 = []; alertCyc_e3 = []; alertCyc_e4 = []; alertCyc_e5 = [];
        alertCyc_e6 = [];alertCyc_e7 = [];alertCyc_e8 = [];alertCyc_e9 = []; alertCyc_e10=[];
        alertCyc_e3_temp = []; %add cc 06/11/2020
        alertCyc_e7_temp = [];
        cycle_proche_surface=[];
        % RECUPERATION donees TRAJ (T) et META (M)
        % ----------------------------------------
        floatname = strtrim(allfloats{1}{k});
        doneT=0;
        doneM=0;
        ii=0;
        elev_all = [];
        %while ~doneT&&~doneM&& ii<=length(DAC)
        while ~doneT&&~doneM&& ii<=length(DAC)-1
            % lecture des traj
            ii=ii+1;
            traj_fileName_R = [DIR_FTP DAC{ii} '/' floatname '/' floatname '_Rtraj.nc'];
            traj_fileName_D = [DIR_FTP DAC{ii} '/' floatname '/' floatname '_Dtraj.nc'];
            if exist(traj_fileName_D,'file');
                [T,DimT,GlobT]=read_netcdf_allthefile(traj_fileName_D);
                traj_fileName_final=[floatname '_Dtraj.nc'];
                doneT=1;
            elseif exist(traj_fileName_R,'file');
                [T,DimT,GlobT]=read_netcdf_allthefile(traj_fileName_R);
                doneT=1;
                traj_fileName_final=[floatname '_Rtraj.nc'];
            end
            meta_fileName=[DIR_FTP DAC{ii} '/' floatname '/' floatname '_meta.nc'];
            prof_fileName=[DIR_FTP DAC{ii} '/' floatname '/' floatname '_prof.nc'];
            if exist(meta_fileName,'file');doneM=1;end
        end
        
        if  doneM
            % lecture des meta
            M=read_file_meta(meta_fileName);
        else
            disp(['NOT FOUND: ' meta_fileName]);
        end
        
        if  ~doneT
            disp(['NOT FOUND: ' traj_fileName_R]);
        else
            
            if isfield(M,'ProfilePressure')
               PARAM.PRESS_PARK_DUMB=max(max(M.ProfilePressure)+400, PARAM.PRESS_PARK_DUMB);
            end
            % RECUPERATION des indices de localisation (idLoc) et de derive en profondeur (idDrift)
            % -------------------------------------------------------------------------------------
            icounterfloat=icounterfloat+1;
            T_sav=T;
            T=format_flags_char2num(T);%change flag char strings to numerical vectors
            %ex '11111441111 ' -> [1 1 1 1 1 4 4 1 1 1 1 999]
            
            T=replace_fill_bynan(T);   % remplace les fillValues (9999.. par des NaN)
            
			%Variables temporaires utiles pour le calcul des vitesses. 
			% rajout de T.position_qc_koba pour stoker l'info sur les resultats des tests de koba
			 T.position_qc_koba=T.position_qc;
			 T.representative_park_temperature.data=NaN*T.representative_park_pressure.data;

		
            % VALEURS DES LOCALISATIONS ARGOS DE SURFACE
            if isfield(T,'measurement_code')==0
                idLoc =   find(~isnan(T.longitude.data)& ~isnan(T.latitude.data));
                fid_alerte=fopen(file_alerte,'a');
                fprintf(fid_alerte,'%s\n',[ floatname ', , Attention pas de measurement_code dans le fichier Traj.'])
                fclose(fid_alerte);
                fprintf('%s\n',[ floatname ' Attention pas de measurement_code dans le fichier Traj.'])
            else
                idLoc =   find(T.measurement_code.data==703);     %%%ne prend pas la première loca en surface car launch
            end
            
            % VALEURS DES MESURES DE PARKING
            % correction cc 21/09/2020 : idDrift  indice des mesures de
            % pression en cours de derive servant a calculer la RPP : ce sont
            % ces mesures qu'on doit flaguer si mauvaises: on récupère tout ce
            % qu'on peut
            
            idDrift=find(T.measurement_code.data==290|T.measurement_code.data==296|T.measurement_code.data==299|T.measurement_code.data==300|T.measurement_code.data==297|T.measurement_code.data==298);
            
            %         if isempty(find(T.measurement_code.data==290, 1)) & isempty(find(T.measurement_code.data==296))==0
            %             idDrift=find(T.measurement_code.data==296);% prend la moyenne des mesures de parking
            %         elseif isempty(find(T.measurement_code.data==290))==0
            %             idDrift=find(T.measurement_code.data==290);% prend l'ensemble des mesures de parking si existe
            %         else
            %             idDrift=find(T.measurement_code.data==299);
            %             if(isempty(idDrift)==0)
            %                  fid_alerte=fopen(file_alerte,'a');
            %                 fprintf(fid_alerte,'%s\n',[floatname ' PAS DE MEASUREMENT_CODE 290 ou 296 (MESURES DE PARKING), RECUPERE ID POUR MEASUREMENT_CODE 299.']);
            %                 fprintf('%s\n',[floatname 'PAS DE MEASUREMENT_CODE 290 ou 296 (MESURES DE PARKING), RECUPERE ID POUR MEASUREMENT_CODE 299.']);
            %                 fclose(fid_alerte);
            %             end
            %         end
            
            if(isempty(idDrift)==1)
                fid_alerte=fopen(file_alerte,'a');
                % correction cc 21/09/2020
                %fprintf(fid_alerte,'%s\n',[floatname ' !!!PAS DE MEASUREMENT_CODE 290 ou 296 (MESURES DE PARKING) OU 299.']);
                %fprintf('%s\n',[floatname '!!!PAS DE MEASUREMENT_CODE 290 ou 296 (MESURES DE PARKING), OU 299.']);
                fprintf(fid_alerte,'%s\n',[floatname ',, !!!PAS DE MEASUREMENT_CODE ASSOCIEE AUX MESURES DE PARKING (290/296/300/297/298)']);
                fprintf('%s\n',[floatname ' !!!PAS DE MEASUREMENT_CODE ASSOCIEE AUX MESURES DE PARKING (290/296/300/297/298)']);
                fclose(fid_alerte);
            end
            if isempty(idLoc)==0
                % on cree la liste des numéros de cycles théoriques (cycles_1) et reel (cycles)
                %-----------------------------------------------------------------------------
                CycLoc=T.cycle_number.data(idLoc);
                numCycleMaxTraj = max(CycLoc);
                numCycleMinTraj = min(CycLoc(CycLoc>=0));% cycle -1 dans les traj: date et position de mise à l'eau
                
                cycles_1 = [numCycleMinTraj:numCycleMaxTraj];% creation de la liste theorique des cycles (incrementation de 1)
                cycles = unique(CycLoc(CycLoc>=0))';% creation de la liste relle avec sauts de cycles
                
                clear CycLoc ;
                
                %On calcule la mediane des durees de cycles (dureeMedianCycle) pour chaque numero de mission.
                % ------------------------------------------------------------------------------------------
                dureeMedianCycle=[];
                end_cycles=[];
                pres_ok = [];
                cycles_m={};
                
                % on selectionne les missions renseignees
                allmiss=T.config_mission_number.data(~isnan(T.config_mission_number.data)); % correction cc T.config_mission_number.data=NaN pour iLaunch
                missions=unique(allmiss,'stable')';
                %missions=unique(T.config_mission_number.data,'stable')';
                %isnot = find(ismember(M.config_mission_number,unique(T.config_mission_number.data))==0);
                isnot = find(ismember(M.config_mission_number,allmiss)==0);
                if(isempty(isnot)==0) M.CycleTime(isnot)=NaN;end
                %%pour les cas ou valeur NaN pour un numero de mission donne dans T.config_mission_number et non dans M.CycleTime
                
                for cm=1:length(cycles)     %%%pour tous les cycles
                    
                    idCyc_m_d=idLoc(T.cycle_number.data(idLoc)==cycles(cm));
                    idCyc_m_p=idDrift(T.cycle_number.data(idDrift)==cycles(cm));
                    isok_d=~isnan(T.juld.data(idCyc_m_d));
                    isok_p=~isnan(T.pres.data(idCyc_m_p));
                    date_ok=T.juld.data(idCyc_m_d(isok_d));
                    end_cycles(cm)=date_ok(end);
                    
                end
                
                idCyc_m_d=idLoc(T.cycle_number.data(idLoc)==cycles(1));
                isok_d=~isnan(T.juld.data(idCyc_m_d));
                date_ok_1=T.juld.data(idCyc_m_d(isok_d));
                
                clear idCyc_m_d idCyc_m_p isok_d isok_p;
                
                duree_cycle(2:length(end_cycles)) = diff(end_cycles)./diff(cycles);    %%duree de chaque cycle  (sauf le cycle 0 ou le cycle 1 pour flotteur APEX)
                if cycles(1)==0|cycles(1)==1
                    duree_cycle(1)=date_ok_1(end)-date_ok_1(1);
                else
                    duree_cycle(1)=NaN;
                end
                
                
                for m = 1:length(missions)
                    
                    if(isnan(missions(m))==0)
                        
                        idcMis = find(T.config_mission_number.data==missions(m));
                        Tcycun  = unique(T.cycle_number.data);Tcycun = Tcycun(Tcycun>=0);  %%%ne prend pas les cycles -1
                        idcMisok  = find(ismember(cycles,Tcycun(idcMis))==1);                %% prend en compte les cas où le numéro de cycle est dans T.config_number mais pas dans cycles car pas de MC=703
                        
                        if(~isempty(idcMisok))
                            cycles_m{m}= cycles(idcMisok);
                            
                        else
                            cycles_m{m} = NaN;     %%% si pas de cycle attribué au n° de mission (peut arriver quand cyle 0 mais pas de Measurement Code 703)
                            
                        end
                        cyclesm = cycles_m{m};
                        index = find(cyclesm==0);
                        
                        % if(~isempty(index)) cyclesm(index)=NaN; cycles_m{m} = cyclesm; end    %%ne regarde qu'a� partir du cycle 1.
                        
                        if(length(find(~isnan(cycles_m{m})))>=1)    %pour les cas ou� cycles_m ~= juste le cycle 0
                            Begin = cycles(1);
                            %if Begin == 0
                            %    duree_cycle_m = duree_cycle(find(ismember(cycles,cycles_m{m})==1)-1);
                            %else                                   %%%pour le cas des fotteurs où cycles(1) = 1 comme apex ou cycles(1)=2,3...etc
                            %    cycles_m_Withoutfirst = cycles_m{m};
                            %    cycles_m_Withoutfirst = cycles_m_Withoutfirst(cycles_m_Withoutfirst>Begin);
                            %   duree_cycle_m = duree_cycle(find(ismember(cycles,cycles_m_Withoutfirst)==1)-1);
                            %    if cycles_m{m}==1; duree_cycle_m = NaN;end
                            %end
                            duree_cycle_m = duree_cycle(find(ismember(cycles,cycles_m{m})==1));
                            end_cycles_m = end_cycles(find(ismember(cycles,cycles_m{m})==1));
                            
                            %% calcul de la mediane des durees de cycles
                            
                            if length(duree_cycle_m)==1
                                dureeMedianCycle(m)=duree_cycle_m;
                                %elseif length(duree_cycle_m)~=1 & mynanstd(duree_cycle_m)<100
                            elseif length(duree_cycle_m)~=1 & mynanstd(duree_cycle_m)<PARAM.TIME_STD_DUREE_CYCLE_M
                                dureeMedianCycle(m)=medianoutnan(duree_cycle_m);
                            elseif isnan(duree_cycle_m)==1
                                dureeMedianCycle(m) = NaN;
                            elseif isempty(duree_cycle_m)==1
                                dureeMedianCycle(m) = NaN;
                            end
                        else
                            duree_cycle_m = NaN;
                            
                            dureeMedianCycle(m) = NaN;
                            
                        end
                        
                        
                        if length(missions)<=length(M.CycleTime)
                            
                            %if abs(dureeMedianCycle(m)-M.CycleTime(missions(m)))>0.1
                            if abs(dureeMedianCycle(m)-M.CycleTime(missions(m)))>PARAM.TIME_DIFF_CYCLE
                                fid_alerte=fopen(file_alerte,'a');
                                fprintf(fid_alerte,'%s\n',[ floatname ', , INCOHERENCE DUREE CYCLE ENTRE META, (', num2str(M.CycleTime(m)),' j) ET TRAJ (',num2str(dureeMedianCycle(m)), ' j)']);
                                fclose(fid_alerte);
                                fprintf('%s\n',[ floatname ', INCOHERENCE DUREE CYCLE ENTRE META (',num2str(M.CycleTime(m)),' j) ET TRAJ (',num2str(dureeMedianCycle(m)), ' j)']);
                                if(Stat==1)
                                    alertCyc_e10 = [alertCyc_e10 cycles_m{m}];
                                    alerte1(k)=str2double(floatname);
                                    
                                end
                            end
                            
                            %%%pour STATS
                            ite = ite+1;
                            Diff_Medcyc_Mcyc(k,m) = dureeMedianCycle(m)-M.CycleTime(missions(m));
                            
                            idLoc(T.cycle_number.data(idLoc)==cycles(cm));
                        else
                            Diff_Medcyc_Mcyc(k,m)=NaN;
                            
                        end
                        
                    else
                        dureeMedianCycle(m)=NaN;
                        Diff_Medcyc_Mcyc(k,m)=NaN;
                    end      %%fin de la boucle sur la condition isnan
                    
                end   %%% fin de la boucle sur les missions
                
                if length(missions)>length(M.CycleTime)
                    
                    fid_alerte=fopen(file_alerte,'a');
                    fprintf(fid_alerte,'%s\n',[ floatname ',, INCOHERENCE META & TRAJ,: Nbre de  MISSIONS DANS META (', num2str(length(M.CycleTime)),' ) < a celui dans TRAJ (', num2str(length(missions)), ' ). Ne peut pas correctement verifier la coherence ni pour la duree des cycles ni pour P_park.']);
                    fclose(fid_alerte);
                    fprintf('%s\n',[ floatname ', INCOHERENCE META & TRAJ: Nbre de MISSIONS DANS META (',num2str(length(M.CycleTime)),' ) < a celui dans TRAJ (', num2str(length(missions)), ' ). Ne peut pas correctement verifier la coherence ni pour la duree des cycles ni pour P_park.']);
                    if (Stat==1)
                        alertCyc_e10 = [alertCyc_e10 999];
                    end
                end
                
                
                
                %% -------------------------------------------------------------
                % alerte 3, 4, 5, 6
                % BOUCLE SUR CHAQUE CYCLE POUR VERIFIER:
                % - double de cycles
                % - pb croissance cycle
                % - mauvaise duree du cycle
                
                [cycles_sorted,idCyc_sorted]= sort(cycles);
                
                
                [o_alerte3, o_alerte4, o_alerte5,...
                    o_alerte6,o_alertCyc_e4, o_alertCyc_e5,idCycprec]= Test_cycles_new(idLoc, cycles, cycles_1,cycles_sorted,duree_cycle, ...
                    missions, dureeMedianCycle);
                if(isempty(o_alerte3)==0)  alerte3{k} = o_alerte3;end;
                if(isempty(o_alerte4)==0)  alerte4{k}  = o_alerte4; end;
                if(isempty(o_alerte5)==0)  alerte5{k} = o_alerte5;end;
                if(isempty(o_alerte6)==0) alerte6{k} = o_alerte6; end;
                
                if(isempty(o_alertCyc_e4)==0)  alertCyc_e4 = [alertCyc_e4 o_alertCyc_e4]; end;
                if(isempty(o_alertCyc_e5)==0)  alertCyc_e5 = [alertCyc_e5 o_alertCyc_e5]; end;
                
                
                %% ---------------------------------------------------------------
                % alerte 7, 8,9
                % VERIFICATION DES DATES DE LOC ARGOS: DATES REALISTES: later than 1 st January 1997
                % VERIFICATION DES LOC ARGOS:
                % VERIFICATION DES DATES ET LOCS DE LANCEMENT (coherence entre les fichiers de trajectoires et de metadonnees)
                
                iLaunch=find(T.cycle_number.data==-1);
                LaunchDate=T.juld.data(iLaunch);
                LaunchLat=T.latitude.data(iLaunch);
                LaunchLon=T.longitude.data(iLaunch);
                
                [o_alerte7, o_alerte8, o_alerte9, o_alertCyc_e5,o_alertCyc_e6, o_alertCyc_e9, o_Diff_LaunchD,...
                    o_Diff_LaunchLat, o_Diff_LaunchLon] =Test_locs(idLoc,...
                    cycles_sorted, iLaunch, LaunchDate, LaunchLat, LaunchLon);
                
                if(isempty(o_alerte7)==0)  alerte7{k} = o_alerte7;end;
                if(isempty(o_alerte8)==0)  alerte8{k}  = o_alerte8; end;
                if(isempty(o_alerte9)==0)  alerte9{k} = o_alerte9;end;
                
                
                if(isempty(o_alertCyc_e5)==0)  alertCyc_e5 = [alertCyc_e5 o_alertCyc_e5]; end;
                if(isempty(o_alertCyc_e6)==0)  alertCyc_e6 = [alertCyc_e6 o_alertCyc_e6]; end;
                if(isempty(o_alertCyc_e9)==0)  alertCyc_e9 = [alertCyc_e9 o_alertCyc_e9]; end;
                
                
                if(isempty(o_Diff_LaunchD)==0) Diff_LaunchD(icounterfloat) = o_Diff_LaunchD;end
                if(isempty(o_Diff_LaunchLat)==0) Diff_LaunchLat(icounterfloat) = o_Diff_LaunchLat;end
                if(isempty(o_Diff_LaunchLon)==0) Diff_LaunchLon(icounterfloat) = o_Diff_LaunchLon;end
                
                %% ---------------------------------------------------------------
                
                % BOUCLE SUR CHAQUE CYCLE POUR VERIFIER:
                % - coherence de la climatologie (temperature, salinite)
                % - calcul pression mediane de derive
                % - coherence climatologie pression de derive.
                % pres_ok: moyenne de la pression a la profondeur de parking par cycle
                % elev_all: bathy la premi�re position de chaque cycle
                % elev_end_all: bathy pour la derni�re position du cycle pr�c�dent
                
                
                idCyc=(1);
                pres_drift_mes=[];
                pres_drift_th=[];
                temp_drift_mes=[];
                temp_drift_th=[];
                psal_drift_mes=[];
                psal_drift_th=[];
                i_rpp=0;
                cyc_wrong_pres=[];
                pres_wrong=[];
                itemp=0;
                cyc_wrong_temp=[];
                temp_wrong=[];
                ipsal=0;
                cyc_wrong_psal=[];
                psal_wrong=[];
                iground=0;
                cyc_wrong_ground=[];
                ground_wrong=[];
                speed_Cyc = {};
                speed_Cyc_isokd= {};
                speed_Cyc_isok={};
                Tab_Num_Cyc = {};
                CyclesGrounded = [];
                it=0;
                pres_alert=0;
                isas_alert=0;
                isas_alert2=0;
                isas_non_ref=0;
                proche_surface=0; % add cc 09/11/2020
                elev_end_all = [];
                pres_drift_mes=[];
                pres_drift_th=[];
                elev_end_all = [];
                elev_all = [];
                
                
                % on convertit les classes Argos en numérique pour Grounded
                GroundedNum = ones(length(T.grounded.data), 1)*-1;
                GroundedNum(find(T.grounded.data == 'Y')) = 1;
                GroundedNum(find(T.grounded.data == 'N')) = 2;
                GroundedNum(find(T.grounded.data == 'U')) = 3;
                
                i_rpp=0;
                itemp=0;
                ipsal=0;
                
                for id=1:length(cycles_sorted)
                    istat=istat+1;
                    
                    idCyc=(1);
                    
                    %temp_drift_mes=[]; % cc comment 01/10/2020
                    temp_drift_th=[];
                    %psal_drift_mes=[]; % cc comment 01/10/2020
                    psal_drift_th=[];
                    cyc_wrong_temp=[];
                    temp_wrong=[];
                    cyc_wrong_psal=[];
                    psal_wrong=[];
                    speed_Cyc = {};
                    speed_Cyc_isokd= {};
                    speed_Cyc_isok={};
                    Tab_Num_Cyc = {};
                    it=0;
                    
                    
                    
                    
                    %istat=istat+1; % cc dejà incrementé ligne 570??
                    % numero du cycle
                    numCycle = cycles_sorted(id);
                    numMis=T.config_mission_number.data(id);
                    % id des mesures associees  a ce cycle
                    isCyc = (T.cycle_number.data(idLoc) == numCycle);
                    isCyc_drift = (T.cycle_number.data(idDrift) == numCycle);
                    idCyc = idLoc(isCyc);
                    idCyc_drift=idDrift(isCyc_drift);
                    
                    locDate = T.juld.data(idCyc);
                    locLon = T.longitude.data(idCyc);
                    locLat = T.latitude.data(idCyc);
                    locQc = T.position_accuracy.data(idCyc);
                    locTemp = T.temp.data(idCyc);
                    locDate_qc=T.juld_qc.data(idCyc);
                    locPosition_qc=T.position_qc.data(idCyc);
                    
                    pres_drift_mes2(id)=abs(mynanmean(T.pres.data(idCyc_drift)));   %%a supprimer
                    temp_drift_mes2(id)=mynanmean(T.temp.data(idCyc_drift));
                    
                    
                    %if ~isempty(idCyc_drift)&~isnan(T.pres.data(idCyc_drift))  % cc voir comment ca se passe si length(idCyc_drift)>1 attention ca ne passe pas dans boucle si un de pres est NaN
                    %if ~isempty(idCyc_drift)&sum(~isnan(T.pres.data(idCyc_drift)))>0       % correction cc 18/09/2020
                    if ~isempty(idCyc)  % correction cc 03/11/2020 : on ne passe  pas si pas de LOC
                        i_rpp=i_rpp+1;
                        %                 long_first_curr = locLon(1);
                        %                 lat_first_curr = locLat(1);
                        % correction cc 21/09/2020 attention loc argos pas forcement
                        % ranges par date croissantes!
                        [locDateSort,isort]=sort(locDate);
                        long_first_curr = locLon(isort(1));
                        lat_first_curr = locLat(isort(1));
                        
                        LocDateprec = T.juld.data(idCycprec{id});
                        [locDateSort,isortprec]=sort(LocDateprec);
                        
                        lonprec = T.longitude.data(idCycprec{id});
                        lonprec = lonprec(isortprec);
                        long_last_prec= lonprec(max(find(~isnan(lonprec))));
                        latprec = T.latitude.data(idCycprec{id});
                        latprec = latprec(isortprec);
                        lat_last_prec = latprec(max(find(~isnan(latprec)))) ;
                        
                        % cc: correction :21/09/2020 remplace les noms "pres_long"
                        % peu  explicites par les "long_first_curr",
                        % "lat_first_curr" etc
                        
                        
                        %pres_drift_mes(i_rpp)=abs(mynanmean(T.pres.data(idCyc_drift)));
                        [tabFinalParkPres,tabFinalParkTemp,tabFinalParkEtat]= compute_rpp(T,M,idCyc_drift,idCyc,prof_fileName); %add cc 29/09/2020 add a function to compute rpp
                        pres_drift_mes(i_rpp)=tabFinalParkPres;
                        
                        %%   DETERMINE l'elevation a la derniere loc du cycle precedent et premiere loc du cycle courant
                        %   elev_all et elev_end_all
                        
                        if(P.Bathy==1)  % ETOPO
                            
                            ilong = round(mean(find(LONG(:,1)<=long_first_curr+max(diff(LONG)/2) & ...
                                LONG(:,1)>=long_first_curr-max(diff(LONG)/2))));
                            if(long_first_curr>max(LONG(:,1)))
                                ilong = find(LONG(:,1) == max(LONG(:,1)));   %%pour eviter ilong=NaN;
                            end
                            
                            ilat = round(mean(find(LAT(:,1)<=lat_first_curr+max(diff(LAT)/2) & ...
                                LAT(:,1)>=lat_first_curr-max(diff(LAT)/2)))) ;
                            
                            if(lat_first_curr>max(LAT(:,1)))
                                ilong = find(LAT(:,1) == max(LAT(:,1)));       %%%pour eviter ilat=NaN;
                            end
                            
                            elev_first_curr = ELEV(ilat,ilong);     %%%altitude de la premiere pos du cycle n
                            
                            
                            if(id>1)
                                ilong_end = round(mean(find(LONG(:,1)<=long_last_prec+max(diff(LONG)/2) ...
                                    & LONG(:,1)>=long_last_prec-max(diff(LONG)/2))));
                                if(long_last_prec>max(LONG(:,1)))
                                    ilong_end = find(LONG(:,1) == max(LONG(:,1)));   %%pour eviter ilong=NaN;
                                end
                                ilat_end = round(mean(find(LAT(:,1)<=lat_last_prec+max(diff(LAT)/2) ...
                                    & LAT(:,1)>=lat_last_prec-max(diff(LAT)/2))));
                                if(lat_last_prec>max(LAT(:,1)))
                                    ilat_end = find(LAT(:,1) == max(LAT(:,1)));       %%%pour eviter ilat=NaN;
                                end
                                elev_last_prec = ELEV(ilat_end,ilong_end);    %%%altitude de la derniere pos du cycle precedent.
                            end
                            
                        else
                            %%% recupere la bathy de la premiere pos du cycle
                            
                            [o_elev, o_lon, o_lat] = get_srtm_elev(long_first_curr, long_first_curr,...
                                lat_first_curr, lat_first_curr);
                            
                            if(length(o_lon)>1 & length(o_lat)>1 & abs(long_first_curr)<max(abs(o_lon)) ...
                                    & abs(long_first_curr)>min(abs(o_lon)))
                                elev_first_curr = interp2(o_lon,o_lat,o_elev,long_first_curr,lat_first_curr);
                                %                     elseif(length(o_lon)<=1 || length(o_lat)<=1) %remove cc : 15/09/2020
                                %                         elev_first_curr = mean(o_elev);
                            else
                                elev_first_curr = mean(mean(o_elev));
                            end
                            
                            
                            if(id>1 & isempty(idCycprec{id})==0)
                                %%%recupere la bathy de la derniere pos du cycle
                                %%%precedent s'il existe
                                [o_elev_end, o_lon_end, o_lat_end] = get_srtm_elev(long_last_prec, long_last_prec,...
                                    lat_last_prec, lat_last_prec);
                                if(length(o_lon_end)>1 & length(o_lat_end)>1 & abs(long_last_prec)<max(abs(o_lon_end)) ...
                                        & abs(long_last_prec)>min(abs(o_lon_end)))
                                    elev_last_prec = interp2(o_lon_end,o_lat_end,o_elev_end,long_last_prec,lat_last_prec);
                                    %                         elseif(length(o_lon_end)<=1 ||length(o_lat_end)<=1) %remove cc : 15/09/2020
                                    %                             elev_last_prec = mean(o_elev_end);
                                else
                                    elev_last_prec = mean(mean(o_elev_end));
                                end
                            else
                                elev_last_prec = NaN;
                                
                            end
                        end
                        
                        %%% recuperation de l'altitude par i_rpp pour pouvoir la
                        %%% reutiliser dans une autre boucle sur les cycles.
                        elev_all(i_rpp) = elev_first_curr;
                        if exist ('elev_last_prec','var') == 1
                            elev_end_all(i_rpp) = elev_last_prec;
                        else
                            elev_end_all(i_rpp) = NaN;  % add cc 15/09/2020
                        end
                        
                        
                        
                        %% On determine les indices ISAS en longitude, latitude, pression : ilon_drift, ilat_drift, idepth_drift, idepth_std_drift
                        
                        % if max(TEMP_LONG)>=long_first_curr% probleme a la limite du 180°
                        if long_first_curr <= max(TEMP_LONG) & long_first_curr >= min(TEMP_LONG)   % correction cc 15/09/2020: meme probleme a -180°
                            %                     ilong_drift=round(mean(find(TEMP_LONG(:,1)<=long_first_curr+max(diff(TEMP_LONG)/2) ...
                            %                         & TEMP_LONG(:,1)>=long_first_curr-max(diff(TEMP_LONG)/2))));
                            [themin,ilong_drift]=min(abs((TEMP_LONG(:,1)-long_first_curr))); % correction cc 15/09/2020
                        else
                            ilong_drift=1;
                        end
                        %                 ilat_drift = round(mean(find(TEMP_LAT(:,1)<=lat_first_curr+max(diff(TEMP_LAT)/2) ...
                        %                     & TEMP_LAT(:,1)>=lat_first_curr-max(diff(TEMP_LAT)/2)))); % cc: attention l'indice ne correspond pas forcement a la latitude la plus proche
                        [themin,ilat_drift]=min(abs((TEMP_LAT(:,1)-lat_first_curr))); % correction cc 15/09/2020
                        
                        
                        % On fait un premier check des pressions en cours de derive : Test_PTS_isas   add cc 29/09/2020
                        % et on recalcule pres_drift_mes en tenant compte des flags
                        [o_alerte11, o_alerte12,one_isas_alert,one_isas_non_ref,one_proche_surface] = Test_PTS_isas(idCyc_drift, id, cycles_sorted,ilong_drift, ilat_drift,i_rpp);
                        isas_alert(i_rpp)=one_isas_alert;
                        isas_non_ref(i_rpp)=one_isas_non_ref;
                        proche_surface(i_rpp)=one_proche_surface;
                        % if unique(T.cycle_number.data(idCyc_drift))==111
                        % keyboard
                        % end
                        % on recalcule rpp en tenant compte des flags
                        [tabFinalParkPres,tabFinalParkTemp,tabFinalParkPsal,tabFinalParkEtat]= compute_rpp(T,M,idCyc_drift,idCyc,prof_fileName); %add cc 29/09/2020 add a function to compute rpp
                        pres_drift_mes(i_rpp)=tabFinalParkPres;
                        temp_drift_mes(i_rpp)=tabFinalParkTemp;
                        psal_drift_mes(i_rpp)=tabFinalParkPsal;
						T.representative_park_pressure.data(idCyc_sorted(id))=tabFinalParkPres;
                        T.representative_park_pressure_status.data(idCyc_sorted(id))=tabFinalParkEtat;
						T.representative_park_temperature.data(idCyc_sorted(id))=tabFinalParkTemp; % stokage pour fichier yo
                        %if ( isas_alert(i_rpp)==1)
                        %keyboard
                        % fid_alerte=fopen(file_alerte,'a');
                        % fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', SOME DRIFT PRESSURE ARE FLAG 6  (from Test_P_ISAS)']);
                        % fclose(fid_alerte);
                        % fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', SOME DRIFT PRESSURE ARE FLAG 6  (from Test_P_ISAS)']);
                        % T.pres.data(idCyc_drift)
                        % T.pres_qc.data(idCyc_drift)
                        % T.temp.data(idCyc_drift)
                        % T.temp_qc.data(idCyc_drift)
                        
                        % tabFinalParkPres
                        % tabFinalParkTemp
                        % lat_first_curr
                        % long_first_curr
                        % end
                       
                        
                        
                        
                        
                        %% alerte 11, 12
                        % VERIFICATION DE L'ACORD AVEC LA CLIMATOLOGIE POUR T et S
                        % ---------------------------------------------------------
                        % on veut verifier que la rpp qu'on calcule a partir des
                        % pression de derive est coherente avec les données T/S mesuree,
                        % Pour cela, on estime T et S a partir de la clim ISAS à la
                        % profondeur =rpp
                        
                        if max(TEMP_DEPTH)>=pres_drift_mes(i_rpp)
                            %                     idepth_drift=round(mean(find(TEMP_DEPTH(:,1)<=abs(pres_drift_mes(i_rpp))+max(diff(TEMP_DEPTH)/2) ...       % estimation de l'indice : a recalculer
                            %                         & TEMP_DEPTH(:,1)>=abs(pres_drift_mes(i_rpp))-max(diff(TEMP_DEPTH)/2))));
                            [themin,idepth_drift]=min(abs((TEMP_DEPTH(:,1)-pres_drift_mes(i_rpp)))); % correction cc 15/09/2020
                        else
                            idepth_drift=length(TEMP_DEPTH(:,1));
                        end
                        
                        
                        if max(TEMP_DEPTH_STD)>=pres_drift_mes(i_rpp)% certaines profondeurs extremes ne sont pas referencees
                            %                     idepth_std_drift=round(mean(find(TEMP_DEPTH_STD(:,1)<=abs(pres_drift_mes(i_rpp))+max(diff(TEMP_DEPTH_STD)/2) ...
                            %                         & TEMP_DEPTH_STD(:,1)>=abs(pres_drift_mes(i_rpp))-max(diff(TEMP_DEPTH_STD)/2))));
                            [themin,idepth_std_drift]=min(abs((TEMP_DEPTH_STD(:,1)-pres_drift_mes(i_rpp)))); % correction cc 15/09/2020
                        else
                            idepth_std_drift=length(TEMP_DEPTH_STD(:,1));
                        end
                        
                        
                        
                        %                                 [o_alerte11_bis, o_alerte12_bis, isas_alert2, isas_non_ref2] = Test_TS(idCyc_drift, id, cycles_sorted,numMis, missions,i_rpp, idepth_drift,...
                        %                                     ilong_drift, ilat_drift, idepth_std_drift,temp_drift_mes,psal_drift_mes);
                        
                        
                        if(isempty(o_alerte11)==0)  alerte11{k} = o_alerte11;end;
                        if(isempty(o_alerte12)==0)  alerte12{k}  = o_alerte12; end;
                        
                        
                        
                        %% Met le qc à 4 pour les pressions < -2800 et calcule  pres_ok qui correspond a la pression de parking pour le cylcle (id)
                        %  remarque : on calcule pres_ok seulement pour les cycles
                        %  pour lesquel la bathy est assez profonde > 1500m, sinon
                        %  il y a des chances que le flotteur n'ait pas pu
                        %  atteindre sa profondeur de derive
                        
                        pres_drift_th(i_rpp)=-TEMP_DEPTH(idepth_drift);
                        pres_drift_mes(i_rpp)=-pres_drift_mes(i_rpp);
                        
                        
                        %                 if(pres_drift_mes(i_rpp)<=elev_first_curr-PARAM.PRESS_PARK_DIFF_BATHY_QC || pres_drift_mes(i_rpp) <-PARAM.PRESS_PARK_DUMB)
                        %                     T.pres_qc.data(idCyc_drift)=4;
                        %                 end
                        %comentaire cc 21/09/2020 : ce test me semble trop strict:
                        % pres_drift_mes peut etre mesuré tout le long de la derive, ca
                        % peut etre inferieur a la bathy au premier point de loc
                        % (ex d'un flotteur qui arrive sur un talus)
                        % on a deja fait des test sur P dans Test_P_isas
                        % pour calculer pres_ok on se sert de pres_drift_mes
                        % dans pres_ok on recupere la pression de derive si la bathy est
                        % suffisament profonde
                        
                        % remove cc 01/10/2020 : ca ne sert plus
                        %                 qcCyc = T.pres_qc.data(idCyc_drift);
                        %                 isok_p = (qcCyc~=4 & qcCyc~=6 & ~isnan(qcCyc));  % prend en compte le nouveau qc
                        %                 selec = idCyc_drift(isok_p); select = T.pres.data(selec);
                        
                        %recuperation des pressions de derive et moyenne par cycle
                        if(elev_first_curr<-PARAM.PRESS_FOND) & -pres_drift_mes(i_rpp)>PARAM.PRESS_SURF  %ne recupere les pressions que quand la bathy est <-1500
                            pres_ok(id) = -pres_drift_mes(i_rpp);
                            %pres_ok(id) = abs(mynanmean(select(find(select>PARAM.PRESS_SURF))));     %ne prend pas les pressions <50m
                        else
                            pres_ok(id) = NaN;
                        end
                        
                    else
                        pres_ok(id)=NaN;   %%% si le cycle n'a pas de LOC
                    end
                    
                    fclose all;
                    
                end   %%% fin de la boucle sur les cycles (id)
                
                
                
                %% ------------------------------------------------------------------
                
                % on convertit les classes Argos en numerique pour Grounded
                GroundedNum = ones(length(T.grounded.data), 1)*-1;
                GroundedNum(find(T.grounded.data == 'Y')) = 1;
                GroundedNum(find(T.grounded.data == 'N')) = 2;
                GroundedNum(find(T.grounded.data == 'U')) = 3;
                
                
                
                % - Test grounded
                % - presence de doublons de date
                % - croissance des dates
                % - ecart max entre deux loc consecutives
                % - ecart max entre la premiere et derniere loc
                % -  date de lancement anterieure a� permiere loc
                % - loca lisation sur l'ocean
                % - derive en surface et en profondeur
                % - ecart entre les locs de surface
                
                
                
                %%%%% CALCUL DE LA VALEUR MEDIANE DE LA PRESSION PAR numero DE MISSION EN PRENANT EN
                %%%%%  COMPTE LES NOUVEAUX QC  ET TEST SUR LES INCOHERENCES
                %%%%%  TRAJ/META POUR LA PRESSION.
                
                presMedianDrift=[];
                
                %if ~isempty(idDrift)& isempty(find(~isnan(T.pres.data(idDrift)), 1))==0  %% si dispose de mesure en profondeur et si au moins une valeur nonNan parmi les pressions
                if ~isempty(idLoc) % correction cc 03/11/2020 on ne passe pas si pas de loc
                    % par nunero de mission
                    for m = 1:length(missions)
                        
                        idcMis=find(T.config_mission_number.data==missions(m));
                        Tcycun  = unique(T.cycle_number.data);Tcycun = Tcycun(Tcycun>=0);  %%% ne prend pas les cycles -1
                        idcMisok  = find(ismember(cycles,Tcycun(idcMis))==1);                %% prend en compte les cas où le numéro de cycle est dans T.config_number mais pas dans cycles car pas de MC=703
                        
                        % remarque cc: attention ici  idcMisok n'est pas un indice
                        % sur les cycles croissants, contrairement a pres_ok: voir
                        % si ca ne pose pas pb
                        
                        if(~isempty(idcMisok))      %%% si pas de cycle attribué au n� de mission (peut arriver quand cyle 0 mais pas de MC 703)
                            pres_ok_m{m} = pres_ok(idcMisok);
                        else
                            pres_ok_m{m} = NaN;
                        end
                        
                        pres = pres_ok_m{m};
                        
                        %% Calcul de la pression mediane par numero de mission
                        if(~isempty(idcMisok))
                            if sum(~isnan(cycles_m{m}))>=1    %% pour les cas ou cycles_m{m} ~= juste le cycle 0
                                %%prise en compte du qc mis suite a� la verif avec la bathy
                                %  cc 21/09/2020 : tester avec PARAM.PRESS_STD plus
                                %  grand
                                if mynanstd(pres(pres>PARAM.PRESS_STD_MIN&pres<PARAM.PRESS_STD_MAX))<PARAM.PRESS_STD
                                    %presMedianDrift(m)=-abs(median(pres_ok_m(pres_ok_m>500&pres_ok_m<1500),2));    %%qu'en est-il quand Ppark ~ 300m ?
                                    presMedianDrift(m)=-abs(median(pres(pres>PARAM.PRESS_STD_MIN&pres<PARAM.PRESS_STD_MAX),2));
                                else
                                    presMedianDrift(m)= -M.ParkPressure(missions(m));
                                end
                                
                                
                            else
                                presMedianDrift(m) = -M.ParkPressure(missions(m));
                            end
                        else
                            presMedianDrift(m) = NaN;
                        end
                        
                        %%pour stats
                        presall_traj(icounterfloat,m) = presMedianDrift(m) ;
                        if(isnan(missions(m))==0)
                            if(length(missions)<=length(M.ParkPressure))
                                ParkPress_all(icounterfloat,m) = -M.ParkPressure(missions(m));
                            end
                        else
                            ParkPress_all(icounterfloat,m) = NaN;
                        end
                        
                        %%%Comparaison entre pression mediane et pression indiquee dans fichier meta
                        
                        if length(missions)<=length(M.ParkPressure) & isnan(presMedianDrift(m))==0
                            
                            if (abs(presMedianDrift(m)- (-M.ParkPressure(missions(m))))>PARAM.PRESS_PARK_DIFF_M) % NC laisser une marge d'erreur
                                fid_alerte=fopen(file_alerte,'a');
                                fprintf(fid_alerte,'%s\n',[ floatname ', ,INCOHERENCE ENTRE META ET TRAJ POUR P_park : |TPmed-MP| >30 db']);
                                fclose(fid_alerte);
                                fprintf('%s\n',[ floatname ', ,INCOHERENCE ENTRE META ET TRAJ POUR P_park : |TPmed-MP| > 30 db']);
                                if Stat==1
                                    alerte2(k)=str2double(floatname);
                                    
                                end
                                
                            end
                            %%%pour STATS
                            if Stat==1
                                presall_met(icounterfloat,m) = M.ParkPressure(missions(m));
                                Diff_Ppark(icounterfloat,m) = abs(presMedianDrift(m) - (-M.ParkPressure(missions(m))));
                            end
                        else
                            if Stat==1
                                Diff_Ppark(icounterfloat,m)= NaN;
                                presall_met(icounterfloat,m) = NaN;
                            end
                        end
                        
                        %%%recuperation de la profondeur de parking a partir de
                        %%%l'info dans les meta si elle existe. sinon prend la
                        %%%valeur mediane calculee pour chaque mission. On prend
                        %%%aussi cette derniere dans le cas ou la valeur de la
                        %%%pression mediane est bien differente de la valeur indiquee dans
                        %%%les meta.
                        
                        if m<=length(M.ParkPressure) & abs(presMedianDrift(m) - -M.ParkPressure(m)) <  PARAM.PRESS_PARK_DIFF_M    %%si nbre de mission dans meta et traj est <= et si la mediane de la pression est proche de celle indiquee dans les meta
                            park_prof(m) = -M.ParkPressure(m);
                        elseif m>length(M.ParkPressure) || presMedianDrift(m)<-M.ParkPressure(m)-PARAM.PRESS_PARK_DIFF_M ...  %%%si nbre de mission dans meta manquant ou si la mediane de la pression trop eloignee de celle indiquee dans les meta
                                || presMedianDrift(m)>-M.ParkPressure(m)+PARAM.PRESS_PARK_DIFF_M
                            
                            park_prof(m) = presMedianDrift(m);
                        end
                        
                    end  %%fin de la boucle sur les missions
                    
                end   %%fin de la condition sur idLoc
                
                %%% alerte 10, 13, 14
                %%%%BOUCLE SUR LES CYCLES
                i_rpp=0;
                istat=0;
                ipres=0;
                %keyboard
                for id=1:length(cycles_sorted)
                    
                    istat=istat+1;
                    % numero du cycle
                    numCycle = cycles_sorted(id);
                    numMis=T.config_mission_number.data(id);
                    idMis = find(double(missions)==numMis);
                    % id des mesures associees  à ce cycle (et au cycle suivant
                    % pour utilisation ulterieure
                    isCyc = (T.cycle_number.data(idLoc) == numCycle);
                    isCyc_post = (T.cycle_number.data(idLoc) == numCycle+1);
                    isCyc_drift = (T.cycle_number.data(idDrift) == numCycle);
                    idCyc = idLoc(isCyc);
                    idCyc_post = idLoc(isCyc_post);
                    idCyc_drift=idDrift(isCyc_drift);
                    
                    locDate = T.juld.data(idCyc);
                    locDate_post = T.juld.data(idCyc_post);
                    locLon = T.longitude.data(idCyc);
                    locLat = T.latitude.data(idCyc);
                    locQc = T.position_accuracy.data(idCyc);
                    locDate_qc=T.juld_qc.data(idCyc);
                    locPosition_qc=T.position_qc.data(idCyc);
                    
                    iground=0;
                    pres_alert=0;
                    alerte_grounded = 0;
                    CyclesGrounded = [];
                    cyc_wrong_ground=[];
                    ground_wrong=[];
                    cyc_wrong_pres=[];
                    pres_wrong=[];
                    
                    %             if ~isempty(idCyc_drift)&~isnan(T.pres.data(idCyc_drift))
                    %if ~isempty(idCyc_drift)&sum(~isnan(T.pres.data(idCyc_drift)))>0       % correction cc 21/09/2020
                    if ~isempty(idCyc)% correction cc 03/11/2020  on ne passe pas si pas de loc
                        i_rpp=i_rpp+1;
                        
                        %pres_drift_mes(i_rpp)=abs(mynanmean(T.pres.data(idCyc_drift)));
                        
                        % Verification de la pression de derive  (par rapport à la
                        % pression mediane)     (remplacer par park_prof ?)
                        
                        if (pres_drift_mes(i_rpp)<presMedianDrift(idMis)-PARAM.PRESS_PARK_DIFF_ISAS ...
                                | pres_drift_mes(i_rpp)>presMedianDrift(idMis)+PARAM.PRESS_PARK_DIFF_ISAS)        %| presMedianDrift(idMis)<=elev_all(i_rpp)-50   %%%< elev_first_curr-100 plutôt ? ou enlever cette condition car vérif grounded plus loin
                            
                            % commentaire cc 21/09/2020 attention ici
                            % pres_drift_mes n'est pas forcement aberrantes: le
                            % flotteur peut simplement passer sur de la bathy, il
                            % ne descend pas a sa prof de parking
                            
                            ipres=ipres+1;
                            % recuperation des donnees aberrantes
                            cyc_wrong_pres(ipres)=i_rpp;
                            pres_wrong(ipres)=pres_drift_mes(i_rpp);
                            fid_alerte=fopen(file_alerte,'a');
                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DIFF PRESSION / PROFONDEUR PARKING , ( ' num2str(pres_drift_mes(i_rpp)) 'm) ~= Mediane: ' num2str(presMedianDrift(idMis)) 'm.']);
                            fclose(fid_alerte);
                            fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DIFF PRESSION / PROFONDEUR PARKING , (' num2str(pres_drift_mes(i_rpp)) 'm) ~= Mediane: ' num2str(presMedianDrift(idMis)) 'm.']);% NC test très large, teste l'ecart à la pression de derive mediane
                            %keyboard
                            pres_alert(i_rpp)=1;
                            if Stat==1
                                % alertCyc_e7 = [alertCyc_e7 cycles_sorted(id)];
                                alerte10(k,cycles_sorted(id)+1)=str2double(floatname);
                            end
                        else
                            pres_alert(i_rpp)=0;
                            it = it+1;
                            if Stat==1 & isempty(idMis)==0
                                Diff_Pm_Pmed(it) =  pres_drift_mes(i_rpp)-presMedianDrift(idMis);
                            end
                        end
                        %keyboard
                        if  isas_alert(i_rpp)==1  % modif cc 26/10/2020
                            %if pres_alert(i_rpp)==1 & isas_alert(i_rpp)==1
                            %                     fid_alerte=fopen(file_alerte,'a');
                            %                     fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB CAPTEUR PRESSION ? pression non compatible avec pression de parking et temperature']);
                            %                     fclose(fid_alerte);
                            %                     fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB CAPTEUR PRESSION ? pression non compatible avec pression de parking et temperature']); % NC verifier si on doit garder cette alerte sous cette forme
                            if Stat==1
                                alertCyc_e7_temp = [alertCyc_e7_temp cycles_sorted(id)];
                                alerte13(k,id)=str2double(floatname);
                            end
                            
                        end
                        if proche_surface(i_rpp)==1
                            cycle_proche_surface = [cycle_proche_surface cycles_sorted(id)];  % cc 09/11/2020 recuperation des cycles pour lesquels le flotteur reste en surface pour les ajouter aux alertes pression lors de la comparaison avec Altran.
                        end
                        
                        
                        if  (elev_end_all(i_rpp)>=park_prof(idMis)+PARAM.PRESS_PARK_DIFF_BATH | elev_all(i_rpp)>=park_prof(idMis)+PARAM.PRESS_PARK_DIFF_BATH) ...
                                | (isas_non_ref(i_rpp)==1 & pres_alert(i_rpp)==1)
                            %if(GroundedNum(id)~=2)    %%% 2 = 'N'	% correction cc 03/11/2020 Si le flotteur dit qu'il n'est pas grounded alors, pas d'alerte grounded
                            alerte_grounded = 1;
                            iground=iground+1;
                            ifloat = ifloat+1;
                            % recuperation des donnees aberrantes
                            cyc_wrong_ground(iground)=i_rpp;
                            ground_wrong(iground)=pres_drift_mes(i_rpp);
                            %float_groundedgrounded(ifloat) = str2double(floatname);
                            %CyclesGrounded(iground) = cycles_sorted(id);
                            fid_alerte=fopen(file_alerte,'a');
                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', GROUNDED ?']);
                            fclose(fid_alerte);
                            fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', GROUNDED ?']);
                            %%% appliquer flags ?
                            if Stat==1
                                alertCyc_e3_temp = [alertCyc_e3_temp cycles_sorted(id)];   % on elimine les locterrre pour comparaison avec altran
                                alerte14(k,cycles_sorted(id)+1)=str2double(floatname);
                            end
                            %fprintf(fid3,'%s\n',[floatname ', [' num2str(cycles_sorted(id)) ']']);
                            
                            %T.grounded.data(id) = 'Y'; % add cc 05/10/2020
                            T.grounded.data(idCyc_sorted(id)) = 'B'; % correction cc 08/01/2021
                            %else
                            if(GroundedNum(idCyc_sorted(id))==2)
                                fid_alerte=fopen(file_alerte,'a');
                                fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', grounded selon bathy mais incohérence avec la variable T.grounded =''N''.']);
                                fclose(fid_alerte);
                                fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', grounded selon bathy mais incohérence avec la variable T.grounded =''N''.']);                      %%% appliquer flags ?
                                %%%mettre la variable grounded à 'Y' ??
                                %T.grounded.data(id) = 'Y'; % remove cc 05/10/2020
                            end
                            
                            if(map==1)
                                m_plot(TEMP_LONG(ilong_drift,1),TEMP_LAT(ilat_drift,1),'color',[0.4660 0.6740 0.1880],'marker','p','markersize',8)
                            end
                        elseif(GroundedNum(idCyc_sorted(id))==1)    %%%1 = 'Y'
                            % alertCyc_e3 = [alertCyc_e3 cycles_sorted(id)]; %
                            % remove cc 05/10/2020
                            % fid_alerte=fopen(file_alerte,'a');
                            % fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', pas d''alerte grounded alors que oui d''apres la variable T.grounded.']);
                            % fclose(fid_alerte);
                            % fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', pas d''alerte grounded alors que oui d''apres la variable T.grounded.']);                      %%% appliquer flags ?
                            T.grounded.data(idCyc_sorted(id))='N';   %  cc 15/01/2021
                        elseif(GroundedNum(idCyc_sorted(id))==2)    %%%2 = 'N'  % add cc 05/10/2020
                            T.grounded.data(idCyc_sorted(id))='N';
                        %elseif(GroundedNum(idCyc_sorted(id))==3)       %%%3 = 'U'
						else       %%%'U' and all other
                            %                     fid_alerte=fopen(file_alerte,'a');
                            %                     fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', pas d''alerte grounded et ''U'' d''apres la variable T.grounded.']);
                            %                     fclose(fid_alerte);
                            % fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', pas d''alerte grounded et ''U'' d''apres la variable T.grounded.']);                      %%% appliquer flags ?
                            T.grounded.data(idCyc_sorted(id))='N';
                        end
                        
                        
                        
                    end   %% fin de la condition sur l'existence des LOC
                    
                    if(isnan(T.config_mission_number.data(idCyc_sorted(id)))==0)
                        %%id mission du flotteur pour ce cycle
                        numMis=T.config_mission_number.data(idCyc_sorted(id));   %%%attention pour id=1 peut etre = NaN;
                        idMis = find(missions==numMis);
                        % recuperation des durees de cycle
                        if idMis<=length(dureeMedianCycle)
                            dureeCycle=dureeMedianCycle(idMis);
                        else
                            dureeCycle = duree_cycle(id);              %%%recupere duree_cycle quand M.CycleTime mal renseigné donc dureeMedianCycle aussi.
                        end
                    else
                        dureeCycle = NaN;
                    end
                    
                    % choix des bornes des alertes en fonction du type de
                    % positionnement
                    if strfind('IRIDIUM',M.trans_system)
                        bornesurf=20000;
                        ecartcycle=PARAM.TIME_DIFF_FLLOCS_IR;
                    elseif strfind('ARGOS',M.trans_system)
                        bornesurf=100000;
                        ecartcycle=PARAM.TIME_DIFF_FLLOCS_AR;
                    end
                    
                    
                    if (sum(idCyc) > 0)
                        
                        
                        % VERIF DOUBLE DE DATE DE LOC
                        %-----------------------------
                        % trouve les doubles de date (qd <=30s) et choisi la date
                        % A� garder selon la precision du positionnement ARGOS
                        % voir ce que ca donne avec iridium!  (idDoublon_date
                        % renseigne uniquement si idDoublon_latlon associe)
                        
                        %[isdouble_toremove,idDoublon_date,idDoublon_latlon] = ...
                        %    find_doubles_datelatlon(locDate, locLon, locLat, locQc, locTemp);   %%%recherche les doublons. Attention, il peut y avoir un double de date pas forcement associe à double de lon ou de lat
                        
                        % cc 05/10/2020 : si critere de double de date <30s on
                        % trouve enormement de doubles. Peut être normal
                        % Je pense qu'on peut se contenter de checher les
                        % doubles exacts seulement
                        % nouveau parametre :PARAM.TIME_DIFF_DOUBLE_DATE_LOC=0
                        % pour pouvoir regler cette difference
                        
                        [isdouble_toremove,idDoublon_date,idDoublon_latlon] = ...
                            find_doubles_datelatlon_30s(locDate, locLon, locLat, locQc, locTemp);   %%%recherche les doublons. Attention, il peut y avoir un double de date pas forcement associe à double de lon ou de lat
                        
                        
                        
                        if isempty(isdouble_toremove)==0
                            
                            %keyboard
                            %
                            
                            fid_alerte=fopen(file_alerte,'a');
                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DOUBLE DATES DE LOC']);
                            fclose(fid_alerte);
                            fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DOUBLE DATES DE LOC']);
                            
                            %T.position_qc.data(idCyc(idDoublon_date(find(T.position_qc.data(idCyc(idDoublon_date)))>1)))=6;   %%flag à 6 les dates doublées dont le qc est >1
                            % correction cc 24/09/2020
                            isdb = zeros(length(idCyc),1);
                            isdb(isdouble_toremove)=1;
                            isflag4=T.juld_qc.data(idCyc)==4;
                            T.juld_qc.data(idCyc(isdb&~isflag4))=6;
                            isflag4=T.position_qc.data(idCyc)==4;
                            T.position_qc.data(idCyc(isdb&~isflag4))=6;
                            %T.position_qc.data(idCyc(idDoublon_date(find(T.position_qc.data(idCyc(idDoublon_date))>1))))=6;
                            
                            if Stat==1
                                %alertCyc_e9 = [alertCyc_e9 cycles_sorted(id)];  cc 06/11/2020 ce n'est pas dans les alertes locpos altran
                                alerte15(k,cycles_sorted(id)+1)=str2double(floatname);
                            end
                            
                            % comment cette partie cc 05/10/2020
                            %                     if isequal(idDoublon_date,idDoublon_latlon)
                            %                         fid_alerte=fopen(file_alerte,'a');
                            %                         fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ' id LOC ' num2str(idDoublon_latlon(1)) '&' num2str(idDoublon_latlon(2))  ', DOUBLE DATES DE LOC ASSOCIE A DOUBLES LON-LAT']);
                            %                         fclose(fid_alerte);
                            %                         fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ' id LOC ' num2str(idDoublon_latlon(1)) '&' num2str(idDoublon_latlon(2)) ', DOUBLE DATES DE LOC ASSOCIE A DOUBLES LON-LAT']);
                            %                         % T.position_qc.data(idCyc(idDoublon_latlon(find(T.position_qc.data(idCyc(idDoublon_latlon)))>1)))=6;   %%flag à 6 les lon-lat doublés dont le qc est >1
                            %                         % correctio, cc 24/09/2020
                            %                         T.position_qc.data(idCyc(idDoublon_latlon(find(T.position_qc.data(idCyc(idDoublon_latlon))>1))))=6;
                            %                         if Stat==1
                            %                             alertCyc_e9 = [alertCyc_e9 cycles_sorted(id)];
                            %                         end
                            %                     end
                            
                            % question cc 05/10/2020? si on a un doublon de position sans doublons de
                            % date pour un meme cycle?  que fait-on? => les
                            % positions sont mauvaise?
                            
                            
                        end
                        
                        %                 % QUE FAIRE avec des doubles (FLAG ?)
                        %                 %T.juld_qc.data(isdouble_toremove)=6;
                        %                 T.juld_qc.data(idCyc(isdouble_toremove))=6;
                        
                        locPosition_qc = T.position_qc.data(idCyc);     %%%de nouveau car les flags ont change
                        locDate_qc = T.juld_qc.data(idCyc);             %%%idem
                        
                        %
                        % VERIF CROISSANCE DES DATES DE LOC (qui ont un flag bon)
                        %--------------------------------------------------------
                        isbad = T.juld_qc.data(idCyc)==6|T.juld_qc.data(idCyc)==4;
                        isbad_post = T.juld_qc.data(idCyc_post)==6|T.juld_qc.data(idCyc_post)==4;% add cc 15/09/2020
                        locDateok=locDate(~isbad);
                        locDate_postok=locDate_post(~isbad_post);   % add cc 15/09/2020
                        [LocDateok_sorted, idSorted] = sort(locDateok);
                        [locDate_postok_sorted, idSorted_post] = sort(locDate_postok);% add cc 15/09/2020
                        if isequal(LocDateok_sorted,locDateok)==0
                            fid_alerte=fopen(file_alerte,'a');
                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DATES DE LOC NON CROISSANTES']);
                            fclose(fid_alerte);
                            fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DATES DE LOC NON CROISSANTES']);
                            if Stat==1
                                %alertCyc_e5 = [alertCyc_e5 cycles_sorted(id)];
                                alerte16(k,cycles_sorted(id)+1)=str2double(floatname);
                            end
                        end
                        [LocDate_sorted, idSorted] = sort(locDate);
                        locLon_sorted=locLon(idSorted);
                        locLat_sorted=locLat(idSorted);
                        locQc_sorted=locQc(idSorted);
                        locDate_qc_sorted=locDate_qc(idSorted);
                        locPosition_qc_sorted=locPosition_qc(idSorted);
                        %idCyc_sorted=idCyc(idSorted);
                        
                        
                        %
                        % CALCUL DE L'ECART MAX (en temps) ENTRE DEUX LOCS CONSECUTIVES ARGOS POUR UN MEME
                        % CYCLE (detecter des cycles qui ont ete accoles par erreur, comme
                        %  le cycle correspondant au prelude et le cycle DPF)
                        %----------------------------------------------------------------
                        %if length(LocDate_sorted)>1&(locDate_qc_sorted~=6|locDate_qc_sorted~=4);     % A considere dans les tests sur les numeros de cycles.
                        if length(LocDate_sorted)>1   &    length(LocDate_sorted((locDate_qc_sorted<6)))>1
                            %ecartMeanLocCycle(istat)=median(diff(LocDate_sorted))*24;
                            ecartMeanLocCycle(istat)=median(diff(LocDate_sorted(locDate_qc_sorted<6)))*24;
                            %  ecartMaxLocCycle(istat)= max(diff(LocDate_sorted))*24;
                            ecartMaxLocCycle(istat)= max(diff(LocDate_sorted((locDate_qc_sorted<6))))*24; % en heure
                        else
                            ecartMeanLocCycle(istat)=NaN;
                            ecartMaxLocCycle(istat)=NaN;
                        end
                        
                        if ecartMaxLocCycle(istat)>PARAM.TIME_DIFF_2LOCS %en heure
                            if numCycle==0|numCycle==1   % add cc 25/09/2020
                                fid_alerte=fopen(file_alerte,'a');
                                fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', cycle DPF accolé, ECART MAX 2 LOC > 10h :' num2str(ecartMaxLocCycle(istat)) ' heures']);
                                fclose(fid_alerte);
                                fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', cycle DPF accolé, ECART MAX 2 LOC > 10h :' num2str(ecartMaxLocCycle(istat)) ' heures']);
								
                                if Stat==1
                                    alertCyc_e4 = [alertCyc_e4 cycles_sorted(id)];
                                end
								
								% % cc ajout 10/01/2021 pour ne pas utiliser dans calcul des vitesses (test)
								% igoodd=find(locDate_qc_sorted<6);
								% ijump=find(diff(LocDate_sorted(locDate_qc_sorted<6))==ecartMaxLocCycle(istat)/24);
								% if numCycle==0
								   % T.cycle_number_qc.data(idCyc(idSorted(igoodd(ijump+1:end))))=6;
								% end
								% if numCycle==1
								   % T.cycle_number_qc.data(idCyc(idSorted(igoodd(1:ijump))))=6;
								% end
                            
                                
                            else
                                fid_alerte=fopen(file_alerte,'a');
                                fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ECART MAX 2 LOC > 10h, :' num2str(ecartMaxLocCycle(istat)) ' heures']);
                                fclose(fid_alerte);
                                fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ECART MAX 2 LOC > 10h, :' num2str(ecartMaxLocCycle(istat)) ' heures']);
                            end
                            %keyboard
                            if Stat==1
                                %alertCyc_e5 = [alertCyc_e5 cycles_sorted(id)];
                                % cc remove 26/10/2020 c'est une alerte pour DMOP
                                alerte17(k,id)=str2double(floatname);
                            end
                            
                            
                        end
                        
                        if Stat==1
                            ecartMaxLocCycle_fl(icounterfloat,istat) = ecartMaxLocCycle(istat);
                            ecartMeanLocCycle_fl(icounterfloat,istat) = ecartMeanLocCycle(istat);
                        end
                        
                        % CALCUL DE L'ECART ENTRE FIRST ET LAST LOCS ARGOS POUR UN MEME
                        % CYCLE
                        %----------------------------------------------------------------
                        
                        %keyboard
                        %if (locDate_qc_sorted~=6|locDate_qc_sorted~=4);
                        if length(LocDate_sorted)>1 & ~(locDate_qc_sorted(1)>=4 &&  locDate_qc_sorted(end)>=4) % correction cc 24/09/2020 : si les premieres dates et dernieres dates sont bonnes
                            ecartFirstLastLoc(istat)= LocDate_sorted(end)-LocDate_sorted(1);
                            %if ecartFirstLastLoc(istat)*24 > ecartcycle & id<length(cycles_sorted) % on est pas sur le dernier cycle qui peut etre EOL
                            if ecartFirstLastLoc(istat)*24 > ecartcycle & ecartFirstLastLoc(istat)*24 < PARAM.TIME_DIFF_FLLOCS_EOL & id<length(cycles_sorted) % correction cc 24/09/2020
                                fid_alerte=fopen(file_alerte,'a');
                                fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ECART FIRST LAST LOC >,' num2str(ecartcycle) 'h :' num2str(ecartFirstLastLoc(istat)) ' jours']);
                                fclose(fid_alerte);
                                fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ECART FIRST LAST LOC >,' num2str(ecartcycle) 'h :' num2str(ecartFirstLastLoc(istat)) ' jours']);
                                if Stat==1
                                    %alertCyc_e5 = [alertCyc_e5 cycles_sorted(id)];
                                    % cc remove 26/10/2020
                                    alerte18(k,id)=str2double(floatname);
                                end
                                %elseif ecartFirstLastLoc(istat)*24 > PARAM.TIME_DIFF_FLLOCS_EOL & id<length(cycles_sorted) % en heure
                            elseif (id>1&ecartFirstLastLoc(istat)*24) >= PARAM.TIME_DIFF_FLLOCS_EOL | (ecartFirstLastLoc(istat)*24 > ecartcycle & id==length(cycles_sorted))% en heure   % correction cc 24/09/2020
                                fid_alerte=fopen(file_alerte,'a');
                                fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',? EOL ?']);
                                fclose(fid_alerte);
                     
                                fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',? EOL ?']); % peut-on considerer EOL à partir de 5 j à la surface ?
                                fleet_EOF(k,id)=str2double(floatname);
                                if Stat==1
                                    alertCyc_e4 = [alertCyc_e4 cycles_sorted(id)];
                                end
								
                            end
                            ecartFirstLastLoc(icounterfloat,istat) = ecartFirstLastLoc(istat);
                        end
                        
                        % faire des stats et definir un seuil, ARGOS VS IRIDIUM?,
                        
                        %
                        % CONTROLE DES PREMIERES LOCALISATIONS (MISE EN ROUTE SUR BATEAU)
                        %---------------------------------------------------------------
                        % Il arrive que le flotteur soit mis en route sur le
                        % bateau et localisee par ARGOS avant d'etre mis à l'eau:
                        % on veut eliminer ces localisations qui refletent le
                        % deplacement du bateau
                        
                        
                        % if length(LocDate_sorted)>1&(locDate_qc_sorted~=6|locDate_qc_sorted~=4);
                        %if (numCycleMinTraj==0|numCycleMinTraj==1)&id==1&str2double(M.launch_qc)==1 % on regarde si cycle 0 ou 1
                        if (numCycle==0|numCycle==1)&id==1  % on regarde si cycle 0 ou 1
                            if length(LocDateok_sorted)>=1   % correction cc 15/09/2020 : verifie qu'on a au moins une LOC avec qc valide
                                LaunchDate=T.juld.data(iLaunch);
                                
                                %%%reperer si LaunchDate a une valeur realiste
                                launch_date = LaunchDate+datenum('01011950','ddmmyyyy');
                                is_ok = (launch_date >= datenum('01011997','ddmmyyyy')) & (launch_date < datenum(date));
                                str2double(strtrim(allfloats{1}));
                                
                                if is_ok == 0
                                    fid_alerte=fopen(file_alerte,'a');
                                    fprintf(fid_alerte,'%s\n',[ floatname,', , PB LAUNCH DATE NON REALISTE, ', datestr(launch_date), '. FLAG LAUNCH MIS A 6.']);
                                    fclose(fid_alerte);
                                    fprintf('%s\n',[ floatname, ', ,PB LAUNCH DATE NON REALISTE, ', datestr(launch_date), '. FLAG LAUNCH MIS A 6.']);
                                    %M.launch_qc = 6;
                                    M.launch_qc = '6';          %correction cc 15/09/2020
                                    if T.juld_qc.data(iLaunch)~=4;
                                        T.juld_qc.data(iLaunch)= 6; %add cc 15/09/2020
                                    end
                                    if Stat==1
                                        alertCyc_e6 = [alertCyc_e6 1];
                                        alerte19(k)=str2double(floatname);
                                    end
                                    ecartLaunchDateFirstLoc(k)= NaN;
                                    
                                else
                                    ecartLaunchDateFirstLoc(k)=LocDate_sorted(1)-LaunchDate;
                                    
                                    %% On verifie la Startuo Date (allumage du flotteur qui doit etre anterieur
                                    %% au lancement) dans le cas ou elle est présente dans le fichier meta
                                    if isfield(M,'startup_date')==1
                                        StartDate = datenum(M.startup_date,'yyyymmddHHMMSS');
                                        Diff_LS = LaunchDate+datenum('01011950','ddmmyyyy') - StartDate;
                                        
                                        if M.startup_date_qc=='1' & Diff_LS<0
                                            %M.launch_qc = 6;           %%% Flag à 6 pour la launch date.
                                            M.launch_qc = '6';          %correction cc 15/09/2020
                                            if T.juld_qc.data(iLaunch)~=4;
                                                T.juld_qc.data(iLaunch)= 6; %add cc 15/09/2020
                                            end
                                            fid_alerte=fopen(file_alerte,'a');
                                            fprintf(fid_alerte,'%s\n',[ floatname,', , PB LAUNCH DATE ANTERIEURE A START DATE, DE ', num2str(Diff_LS), ' j. FLAG LAUNCH MIS A 6.']);
                                            fclose(fid_alerte);
                                            fprintf('%s\n',[ floatname, ', , PB LAUNCH DATE ANTERIEURE A START DATE DE ', num2str(Diff_LS), ' j. FLAG LAUNCH MIS A 6.']);
                                            if Stat == 1
                                                alertCyc_e6 = [alertCyc_e6 1];
                                                alerte20(k)=str2double(floatname);
                                            end
                                            
                                        elseif M.startup_date_qc=='1' && Diff_LS>1
                                            fid_alerte=fopen(file_alerte,'a');
                                            fprintf(fid_alerte,'%s\n',[floatname, ', , ECART ENTRE LAUNCH DATE et START DATE > 1jour, PB LAUNCHDATE ?']);
                                            fclose(fid_alerte);
                                            M.launch_qc = '6';          %correction cc 9/10/2020
                                            if T.juld_qc.data(iLaunch)~=4;
                                                T.juld_qc.data(iLaunch)= 6; %add cc 9/10/2020
                                            end
                                            fprintf('%s\n',[ floatname,' ECART ENTRE LAUNCH DATE et START DATE > 1jour, PB LAUNCHDATE ?']);
                                            if Stat == 1
                                                alertCyc_e6 = [alertCyc_e6 1];
                                            end
                                        end
                                    end
                                    % il se peut que la LaunchDate (metadonnee rentree à la main) puisse etre erronee:
                                    % on definit des seuils (en jours) au dela desquels il est fort probable qu'elle soit erronee, auquel cas il
                                    % ne faut pas flagguer a mauvais les dates de loc.
                                    % Si LaunchDate postérieure aux LOCS de plus de 3 jours (ou 3h ???), on se dit que c'est la launchdate qui est
                                    % douteuse (les mises en route sur le bateau se font qq minutes ou heures avant)
                                    %isokLaunchDate(k)=ecartLaunchDateFirstLoc(k)>-PARAM.TIME_DIFF_FLLOCS_AR; %  NC faire des stats sur ces ecarts entre launch date et premiere position: definir un seuil
                                    
                                    isokLaunchDate(k)=ecartLaunchDateFirstLoc(k)*24 > -PARAM.TIME_LAUNCH_FIRST_LOC; % correction cc 15/09/2020: mettre tous les param en heure
                                    
                                    %if isokLaunchDate(k)  & M.launch_qc==1
                                    if isokLaunchDate(k)  & str2double(M.launch_qc)~=6 % correction cc 15/09/2020
                                        
                                        %on verifie que les dates de loc du premier cycle sont sup àla Launch date sinon, on les flaggue
                                        % isbad = T.juld.data(idSorted) <= LaunchDate;
                                        isbad = T.juld.data(idCyc(idSorted)) <= LaunchDate;
                                        isbad_post = T.juld.data(idCyc_post(idSorted_post)) <= LaunchDate;
                                        %T.juld.data(idCyc(idSorted(isbad)))=6;
                                        % T.juld_adjusted_qc(idCyc(idSorted(isbad)))=6;
                                        
                                        if sum(isbad)>0  & sum(isbad_post)== 0  % SI L'ERREUR NE VIENT PAS D'UNE LAUNCH DATE ERRONEE
                                            fid_alerte=fopen(file_alerte,'a');
                                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DATES DE LOC DU 1er CYCLE ANTERIEURES A LAUNCH DATE > 3h']);
                                            fclose(fid_alerte);
                                            fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DATES DE LOC DU 1er CYCLE ANTERIEURES A LAUNCH DATE > 3h']);
                                            %T.juld.data(idSorted(isbad))=2;                 %%%%FLAG � 2 les dates Loc du 1er cycle (douteuses)
											
                                            notQC4=T.juld_qc.data(idCyc(idSorted))~=4;
                                            T.juld_qc.data(idCyc(idSorted(isbad&notQC4)))=6;   % correction cc 15/09/2020
                                            
                                            if Stat == 1
                                                alertCyc_e5 = [alertCyc_e5 cycles_sorted(id)];
                                                alerte21(k,cycles_sorted(id)+1)=str2double(floatname);
                                                
                                            end
                                        end
                                        %elseif isokLaunchDate(k) < -PARAM.TIME_DIFF_FLLOCS_AR       %% SI L'ERREUR VIENT D'UNE LAUNCH DATE ERRONEE
                                    elseif isokLaunchDate(k) ==0 & str2double(M.launch_qc)~=6    %% L'ERREUR PEUT VENIR D'UNE LAUNCH DATE ERRONEE, NON DETECTEE PRECEDEMMENT    %  correction cc 15/09/2020
                                        % fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ? PB LAUNCH DATE POSTERIEURE à DATE DE LOC DE PLUS DE 3 h?']);
                                        fid_alerte=fopen(file_alerte,'a');
                                        fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB LAUNCH DATE ? : POSTERIEURE à DATE DE LOC, DE PLUS DE ' num2str(PARAM.TIME_LAUNCH_FIRST_LOC) ' h']); % correction cc 15/09/2020
                                        fclose(fid_alerte);
                                        %fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ? PB LAUNCH DATE POSTERIEURE à DATE DE LOC DE PLUS DE 3 h?']);
                                        fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB LAUNCH DATE ? : POSTERIEURE à DATE DE LOC DE PLUS DE ' num2str(PARAM.TIME_LAUNCH_FIRST_LOC) ' h']);% correction cc 15/09/2020
                                        
                                        if Stat==1
                                            %alertCyc_e5 = [alertCyc_e5
                                            %cycles_sorted(id)];  % cc 26/10/2020
                                            %ce n'est pas un pb de loc date
                                            alertCyc_e6 = [alertCyc_e6 1]; % add cc
                                            alerte22(k,id)=str2double(floatname);
                                        end
                                    end
                                    
                                    
                                    %if ecartLaunchDateFirstLoc(k)>PARAM.TIME_LAUNCH_FIRST_LOC_DPF % Cas � part ou ecart entre Launch date et dates de LOC tres grand
                                    if ecartLaunchDateFirstLoc(k)*24>PARAM.TIME_LAUNCH_FIRST_LOC_DPF & ecartLaunchDateFirstLoc(k)*24< 30    % correction cc 15/09/2020
                                        %fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ? PB ECART LAUNCH DATE ET DATE DE LOC DE PLUS DE 10.5 j. DPF sans PRELUDE ?']);
                                        %fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ? PB ECART LAUNCH DATE ET DATE DE LOC DE PLUS DE 10.5 j. DPF sans PRELUDE ?']);
                                        fid_alerte=fopen(file_alerte,'a');
                                        fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DPF sans PRELUDE ? : DATE FIRST LOC > LAUNCH DATE, DE PLUS DE '  num2str(PARAM.TIME_LAUNCH_FIRST_LOC_DPF)  'h.']);% correction cc 15/09/2020
                                        fclose(fid_alerte);
                                        fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DPF sans PRELUDE ? : DATE FIRST LOC > LAUNCH DATE DE PLUS DE '  num2str(PARAM.TIME_LAUNCH_FIRST_LOC_DPF)  'h.']);% correction cc 15/09/2020
                                        
                                        if Stat==1
                                            %alertCyc_e4 = [alertCyc_e4 cycles_sorted(id)];
                                            alerte23(k,cycles_sorted(id)+1)=str2double(floatname);
                                            
                                        end
                                        
                                    end
                                    
                                    % add cc 15/09/2020 : j'ajoute cette boucle
                                    if ecartLaunchDateFirstLoc(k)*24>PARAM.TIME_LAUNCH_FIRST_LOC_DUMB
                                        fid_alerte=fopen(file_alerte,'a');
                                        fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB LAUNCH DATE ? : DATE FIRST LOC > LAUNCH DATE, DE PLUS DE '  num2str(PARAM.TIME_LAUNCH_FIRST_LOC_DUMB)  'h.']);
                                        fclose(fid_alerte);
                                        fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB LAUNCH DATE ? : DATE FIRST LOC > LAUNCH DATE DE PLUS DE '  num2str(PARAM.TIME_LAUNCH_FIRST_LOC_DUMB)  'h.']);
                                        M.launch_qc = '6';          %correction cc 9/10/2020
                                        if T.juld_qc.data(iLaunch)~=4;
                                            T.juld_qc.data(iLaunch)= 6; %add cc 9/10/2020
                                        end
                                        %                                 if Stat==1   % cc Alerte A definir
                                        %                                     %alertCyc_e4 = [alertCyc_e4 cycles_sorted(id)];
                                        %                                     alerte23(k,cycles_sorted(id)+1)=str2double(floatname);
                                        %
                                        %                                 end
                                        
                                    end
                                end   % en isok
                                
                                
                                %%%VERIFICATION PHASE DE PRELUDE (si pas agglomere avec cycle 1
                                %%%dans le cas des Deep Profile First (DPF).
                                
                                %%%vérification de la présence du cycle 0 (prélude)
                                
                                % cc 15/09/2020 : je commente cette partie
                                
                                %                         if(isempty(find(cycles==0)))
                                %                             if(T.measurement_code.data(2)==100) %%%DPF (Deep Profile first)
                                %                                 %if ecartLaunchDateFirstLoc(k)*24>24
                                %                                 fprintf(fid_alerte,'%s\n',[ floatname ', cycle '...
                                %                                     num2str(cycles(id)) ', PAS DE CYCLE 0, DPF sans prelude']);
                                %                                 fprintf('%s\n',[ floatname ', cycle '...
                                %                                     num2str(cycles(id)) ', PAS DE CYCLE 0, DPF sans prelude']); % ', LOC plus d''un jour après => DPF sans prelude
                                %                                 if Stat==1
                                %                                     alertCyc_e4=[alertCyc_e4 cycles_sorted(id)];
                                %                                     alerte24(k,cycles_sorted(id)+1)=str2double(floatname);
                                %
                                %                                 end
                                %                             end
                                %                         end
                            else % cc 15/09/2020 :si aucune date de loc valide pour cycle 0 ou 1
                                ecartLaunchDateFirstLoc(k)=NaN;
                            end
                        elseif numCycle>1 & id==1
                            fid_alerte=fopen(file_alerte,'a');
                            fprintf(fid_alerte, '%s\n',[ floatname ', Commence au cycle,  ', num2str(cycles(1)), ' ecartLaunchDateFirstLoc mis a NaN']);
                            fclose(fid_alerte);
                            fprintf('%s\n',[ floatname ', Commence au cycle  ', num2str(cycles(1)), ' ecartLaunchDateFirstLoc mis a NaN']);
                            ecartLaunchDateFirstLoc(k)= NaN;
                        end
                        
                        %                 else  % cc 15/09/2020 :si aucune date de loc valide pour cycle 0 ou 1
                        %
                        %                     %ecartLaunchDateFirstLoc(k)=LocDate_sorted(1)-LaunchDate;
                        %                     ecartLaunchDateFirstLoc(k)=NaN;
                        %
                        %                     % remove cc 15/09/2020 : je commente cette boucle
                        %                     %                     if(ecartLaunchDateFirstLoc(k)>PARAM.TIME_LAUNCH_FIRST_LOC_DUMB)
                        %                     %
                        %                     %                         fprintf(fid_alerte,'%s\n',[ floatname ', cycle '...
                        %                     %                             num2str(cycles(id)) ', UNE SEULE LOC et DATE INCOHERENTE.']);
                        %                     %                         fprintf('%s\n',[ floatname ', cycle '...
                        %                     %                             num2str(cycles(id)) ', UNE SEULE LOC et DATE INCOHERENTE.']);
                        %                     %                         if Stat == 1
                        %                     %                             alertCyc_e4 = [alertCyc_e4 cycles_sorted(id)];
                        %                     %
                        %                     %                         end
                        %                     %                     else
                        %                     %                         fprintf(fid_alerte,'%s\n',[ floatname ', cycle '...
                        %                     %                             num2str(cycles(id)) ', UNE SEULE LOC.']);
                        %                     %                         fprintf('%s\n',[ floatname ', cycle '...
                        %                     %                             num2str(cycles(id)) ', UNE SEULE LOC.']);
                        %                     %                         if Stat == 1
                        %                     %                             alertCyc_e4 = [alertCyc_e4 cycles_sorted(id)];
                        %                     %                             alerte25(k,cycles_sorted(id)+1)=str2double(floatname);
                        %                     %
                        %                     %                         end
                        %                     %                      end
                        %                     %
                        %                     %                     ecartLaunchDateFirstLoc(k) = NaN;
                        %
                        %                 end   %  fin de la condition sur locDateok_sorted
                        
                        
                        % CONTROLE DES POSITIONS
                        %-----------------------------------
                        % si est bien dans l'eau et non sur la terre
                        % on flagge 4 si ce n'est pas le cas
                        
                        iscylocterre=0;
                        isprocheterre=0;
                        for iiloc=1:length(locLat_sorted)  % add boucle cc 23/10/2020
                            %dlat=locLat_sorted(i);    %commentaire cc 25/09/2020 pourquoi ne controller que la 1ere loc du cycle???
                            %dlon=locLon_sorted(1);
                            dlat=locLat_sorted(iiloc);    %commentaire cc 25/09/2020 pourquoi ne controller que la 1ere loc du cycle???
                            dlon=locLon_sorted(iiloc);
                            
                            
                            %%pour recuperer   (temporaire)
                            loclon_first(id) = dlon;
                            loclat_first(id) = dlat;
                            
                            
                            dlatend = locLat_sorted(end);
                            dlonend = locLon_sorted(end);
                            %dlocPosition_qc_sorted = locPosition_qc_sorted(1);
                            dlocPosition_qc_sorted = locPosition_qc_sorted(iiloc);   % cc 23/10/2020
                            dlocPosition_qc_sortedend = locPosition_qc_sorted(end);
                            % on ne verifie que les mesures qui ont passe les tests
                            % temporels
                            % if (locPosition_qc_sorted~=6|locPosition_qc_sorted~=4)&(length(locLon_sorted)>1&length(locLat_sorted)>1); %on ne regarde que les flags corrects ou non traités
                            
                            
                            if ~isnan(dlat)&~isnan(dlon)
                                % cc 23/10/2020 remove condition on qc to get alert
                                % for DM opertor
                                %if (dlocPosition_qc_sorted<4)&(dlocPosition_qc_sortedend<4)% &(length(locLon_sorted)>1&length(locLat_sorted)>1); %on ne regarde que les flags corrects ou non traités (de la 1ère et dernière loc)
                                if(Bathy==1)
                                    ilong = round(mean(find(LONG(:,1)<=dlon+max(diff(LONG)/2) & LONG(:,1)>=dlon-max(diff(LONG)/2))));
                                    ilat = round(mean(find(LAT(:,1)<=dlat+max(diff(LAT)/2) & LAT(:,1)>=dlat-max(diff(LAT)/2))));
                                    elev_float = ELEV(ilat,ilong);
                                else
                                    
                                    
                                    %%%Recupere indices de la premiere loc
                                    [o_elev, o_lon, o_lat] = get_srtm_elev(dlon-0.01, dlon+0.01,...
                                        dlat-0.01, dlat+0.01);
                                    [mini,ilong] = min(abs(o_lon-dlon));
                                    [mini,ilat] = min(abs(o_lat-dlat));
                                    LONG = o_lon'; LAT = o_lat';
                                    elev_float = o_elev(ilat,ilong);
                                    
                                end
                                
                                if ~isempty(ilong)&~isempty(ilat)% condition necessaire ?
                                    
                                    
                                    if ilat>1 & ilat<length(LAT(:,1)) & ilong>1 & ilong<length(LONG(:,1))
                                        
                                        if(Bathy==1)
                                            voisin1=ELEV(ilat(1)-1,ilong(1));
                                            voisin2=ELEV(ilat(1)+1,ilong(1));
                                            voisin3=ELEV(ilat(1),ilong(1)-1);
                                            voisin4=ELEV(ilat(1),ilong(1)+1);
                                        else
                                            voisin1=o_elev(ilat-1,ilong);
                                            voisin2=o_elev(ilat+1,ilong);
                                            voisin3=o_elev(ilat,ilong-1);
                                            voisin4=o_elev(ilat,ilong+1);
                                        end
                                        voisins=[voisin1,voisin2,voisin3,voisin4];
                                        % on regarde si les voisins sont aussi positifs
                                        % au moins 4 voisins doivent avoir une
                                        % elevation positive pour que l'on
                                        % puisse considerer le point comme sur
                                        % le continent    (attention pas le cas
                                        % quand se trouve coinc� dans un atoll
                                        % par exemple)
                                    end
                                    if elev_float>0 % limite par la resolution de la carte à disposition (environ 2°)
                                        
                                        if numel(find(voisins>0))>=4
                                            if (dlocPosition_qc_sorted<4)  % add cc 25/09/2020 fait une alerte que si flag pas deja a 4
                                                fid_alerte=fopen(file_alerte,'a');
                                                fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB DE LOCALISATION SUR LE CONTINENT. FLAG POSITION MIS A 6.']);
                                                fclose(fid_alerte);
                                                fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB DE LOCALISATION SUR LE CONTINENT.FLAG POSITION MIS A 6.']);
                                                if Stat == 1
                                                    alertCyc_e8 = [alertCyc_e8 cycles_sorted(id)];
                                                    alerte26(k,cycles_sorted(id)+1)=str2double(floatname);
                                                end
                                                if(map==1)
                                                    m_plot(LONG(ilong,1),LAT(ilat,1),'color','g','marker','o','markersize',5);
                                                end
												%keyboard
                                                T.position_qc.data(idCyc(idSorted(iiloc)))=6; % correction cc 09/10/2020
                                                %T.position_qc.data(idCyc)=4;
                                            end
                                        else
                                            if iscylocterre==0  % add cc 23/10/2020 1 alerte par cycle pour DM operator
                                                fid_alerte=fopen(file_alerte,'a');
                                                fprintf(fid_alerte,'%s\n',[floatname ', cycle ' num2str(cycles_sorted(id)) ' , LOC A TERRE: altitude du flotteur > 0, (' (num2str(elev_float)) 'm) , mais au moins 1 point voisin < 0. Peut être coincé dans un atoll ou récupéré et proche de la terre']);
                                                fclose(fid_alerte);
                                                fprintf('%s\n',[floatname ', cycle ' num2str(cycles_sorted(id)) ' , LOC A TERRE: altitude du flotteur > 0 (' (num2str(elev_float)) 'm) , mais au moins 1 point voisin < 0. Peut être coincé dans un atoll ou récupéré et proche de la terre']);
                                                if Stat == 1
                                                    alertCyc_e8 = [alertCyc_e8 cycles_sorted(id)];
                                                    alerte27(k,cycles_sorted(id)+1)=str2double(floatname);
                                                    
                                                end
                                                iscylocterre=1;
                                            end
                                            
                                        end
                                        
                                    elseif elev_float<0&elev_float>PARAM.ALT_GROUND_NEAR_COAST &isempty(find(voisins>0))==0    %%%si l'atitude du point est >-200 et l'altitude d'au moins un voisin >0 est probablement proche de la terre.
                                        if isprocheterre==0
                                            fid_alerte=fopen(file_alerte,'a');
                                            fprintf(fid_alerte,'%s\n',[floatname ', cycle ' num2str(cycles_sorted(id)) ' ,  altitude du flotteur proche de 0 ,(' (num2str(elev_float)) 'm) , et au moins 1 point voisin > 0. Probablement proche terre.']);
                                            fclose(fid_alerte);
                                            fprintf('%s\n',[floatname ', cycle ' num2str(cycles_sorted(id)) ' ,  -20 < altitude du flotteur < 0 (' (num2str(elev_float)) 'm) , et au moins 1 point voisin > 0. Probablement proche terre.']);
                                            if Stat == 1
                                                alertCyc_e8 = [alertCyc_e8 cycles_sorted(id)];   
                                            end
                                            isprocheterre=1;
                                        end
                                    end
                                    
                                    
                                    %       if(elev_float>-150 & length(alerte28(k,cycles_sorted(id-15):end))>5 & mean(voisins)>-150)  %%%si elev_first_curr proche surface et plusieurs loc terre aux cycles précédents
                                    %         fprintf(fid_alerte,'%s\n',[floatname ', cycle ' num2str(cycles_sorted(id)) ' ,  altitude du flotteur > -150 (' (num2str(elev_float)) 'm) , et au moins 5 cycles précédents en loc terre. Probablement proche terre.']);
                                    %         fprintf('%s\n',[floatname ', cycle ' num2str(cycles_sorted(id)) ' ,  altitude du flotteur > -150 (' (num2str(elev_float)) 'm) , et au moins 5/15 cycles précédents en loc terre. Probablement proche terre.']);
                                    %         if Stat == 1
                                    %             alertCyc_e8 = [alertCyc_e8 cycles_sorted(id)];
                                    %         end
                                    %      end
                                    
                                    %end
                                    
                                end
                                
                                %                         end    %%fin de la condition sur le qc    %
                                %                         remove condition cc 23/10/2020
                                
                            end       %%% fin de la condition sur dlat et dlon
                            
                        end    %% fin boucle sur toutes les locs
                        
                        
                        
                        % alerte 29,30
                        % CONTROLE DES DERIVES (REALISTES)
                        % --------------------------------
                        % Etablit des statistiques sur la derive du flotteur en
                        % profondeur (entre dernière localisation du cycle N-1 et
                        % la premiere du cycle N) et en surface (entre la premiere
                        % et la derniere localisation du cycle N)
                        % Si superieure aux bornes donnees, on se dit que la
                        % localisation est douteuse
                        %if id>1 & (T.position_qc.data(idCyc)~=6|T.position_qc.data(idCyc)~=4) & ~isnan(dureeCycle)
                        a=6378137;%repere de reference est le WGS84, a� convertir en cartesiennes pour calculer les distances)
                        % dlatlas=dlat;
                        % dlonlas=dlon;
                        
                        % keyboard
                        dlatlas=NaN; %add cc 09/10/2020
                        dlonlas=NaN; %add cc 09/10/2020
                        ddates=NaN;
                        distancederiveprof(1)=NaN; %add cc 12/11/2020
                        distancederivesurf(1)=NaN; %add cc 12/11/2020
                        if  id>1 &  (dlocPosition_qc_sorted<4)&(dlocPosition_qc_sortedend<4) & ~isnan(dureeCycle) & sum(idCycprec{id})>0
                            dlatlas=T.latitude.data(idCycprec{id}(end));
                            dlonlas=T.longitude.data(idCycprec{id}(end));
                            ddates=T.juld.data(idCycprec{id}(end));
                            distanceprof=a*acos(sin(dlat*pi/180)*sin(dlatlas*pi/180)+cos(dlat*pi/180)*cos(dlatlas*pi/180)*cos((dlonlas-dlon)*pi/180));
                            distancederiveprof(istat)=distanceprof; %regroupe les valeurs de dérive en profondeur
                            % la borne varie en fonction de la duree de cycle
                            % calculee et il ne faut pas que la vitesse du
                            % profileur depasse 3m/s
                            
                            
                            %%%%EN PROFONDEUR
                            if distanceprof> dureeCycle*(cycles_sorted(id)-cycles_sorted(id-1))*(24*60*60)*PARAM.SPEED_MAX % borne a adapter en fonction de la duree du cycle (vitesse maximale du flotteur : 3 m/s)
                                fid_alerte=fopen(file_alerte,'a');
                                fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',DERIVE EN PROFONDEUR, > ' num2str(dureeCycle*(cycles(id)-cycles(id-1))*(24*60*60)*PARAM.SPEED_MAX/1000) 'km :' num2str(distanceprof/1000) ' km. Loc Argos probablement erronée.']);
                                fclose(fid_alerte);
                                fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',DERIVE EN PROFONDEUR > ' num2str(dureeCycle*(cycles(id)-cycles(id-1))*(24*60*60)*PARAM.SPEED_MAX/1000) 'km :' num2str(distanceprof/1000) ' km. Loc Argos probablement erronée.']);
                                if Stat == 1
                                    %alertCyc_e9 = [alertCyc_e9 cycles_sorted(id)];  % cc remove 06/11/2020
                                    alerte28(k,cycles_sorted(id)+1)=str2double(floatname);
                                    
                                end
                                if(map==1)
                                    m_plot(dlon,dlat,'color','r','marker','^','markersize',5)
                                    m_plot(dlonlas,dlatlas,'color','r','marker','x','markersize',5)
                                    m_line([dlonlas,dlon],[dlatlas,dlat],'linewi',1,'color','b','linest','--');
                                end
                            end
                            
                            
                            %%%EN SURFACE
                            
                            %%Entre la 1a 1ere et la derniere loc
                            distancesurf=a*acos(sin(dlatend*pi/180)*sin(dlat*pi/180)+cos(dlatend*pi/180)*cos(dlat*pi/180)*cos((dlon-dlonend)*pi/180));
                            distancederivesurf(istat)=distancesurf;%regroupe les valeurs de dérive en surface
                            % la borne varie selon le systeme de positionnement
                            % (1h en surface pour Iridium contre 6h en surface
                            % pour ARGOS)
                            
                            if distancesurf>bornesurf % borne dependant du système de positionnement
                                fid_alerte=fopen(file_alerte,'a');
                                fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DERIVE EN SURFACE ENTRE 1er et DERNIERE LOC D'' UN CYCLE, >' num2str(bornesurf/1000) 'km :' num2str(distancesurf/1000) ' km. Loc Argos probablement erronée.']);
                                fclose(fid_alerte);
                                fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DERIVE EN SURFACE ENTRE 1er et DERNIERE LOC D''UN CYCLE>' num2str(bornesurf/1000) 'km :' num2str(distancesurf/1000) ' km. Loc Argos probablement erronée.']);
                                if Stat == 1
                                    alertCyc_e9 = [alertCyc_e9 cycles_sorted(id)];
                                    alerte29(k,cycles_sorted(id)+1)=str2double(floatname);
                                    
                                end
                                if(map==1)
                                    m_plot(dlonend,dlatend,'color','r','marker','^','markersize',5)
                                    m_plot(dlon,dlat,'color','r','marker','x','markersize',5)
                                    m_line([dlon,dlonend],[dlat,dlatend],'linewi',1,'color','c');
                                end
                            end
                            
                            
                            
                        end   %%fin de la condition sur les qc
                        %end       %%fin de la condition sur idCyc et idCycprec
                        
                        %%% Test de Koba et al. pour localiser les positions Argos
                        %%% erronees a� partir du calcul de vitesse entre les loc de
                        %%% surface. (ne considere pas dans le resultat final les
                        %%% duplicata de dates/lon/lat).
                        %if(id>1 & sum(idCycprec{id})>0)
                        %keyboard
                        fid_koba = fopen([DIR_HOME,'log_koba.txt'],'wt');
                        %locDate_prec = T.juld.data(idCycprec{id}); isok = find(locDate_prec(T.position_qc.data(idCycprec{id})<4));
                        
                        %if(isempty(isok)==0)
                        
                        %%% exclu le cas ou tous les qc seraient a� 4 (qd flotteur sur le continent par ex)
                        %locDate_prec_ok = locDate_prec(isok);
                        %locLat_sorted(find(locPosition_qc_sorted==6))=NaN; locLon_sorted(find(locPosition_qc_sorted==6))=NaN;
                        %locLat_sorted(find(locPosition_qc_sorted==3))=NaN; locLon_sorted(find(locPosition_qc_sorted==3))=NaN; % remarque cc 06/11/2020: des tests koba sont realisés en TR avec flag des positions a 3
                        %locLat_sorted(find(locPosition_qc_sorted==4))=NaN; locLon_sorted(find(locPosition_qc_sorted==4))=NaN;
                        %%test de Koba sur les vitesses
                        [o_date, o_longitude, o_latitude, o_posAcc, o_posQcIn, o_posQcOut ...
                            o_idBadPos] = check_argos_positions_koba_all(LocDate_sorted,...
                            locLon_sorted, locLat_sorted, locQc_sorted, locPosition_qc_sorted, ...
                            ddates, dlonlas, dlatlas, fid_koba);
                        % cc 06/11/2020 add a test that verifie the distance between consecutive Argo LOC (same as the test use by altran to detec bad pos)
                        idGoodPos = find((locQc_sorted == '1') | (locQc_sorted == '2') | (locQc_sorted == '3') | (locQc_sorted == 'G'));
                        locLon = locLon_sorted(idGoodPos);
                        locLat = locLat_sorted(idGoodPos);
                        ranges = distance_lpo(locLat, locLon);
                        ECART_MAX_LOC_ARGOS=30; % remarque pour altran, le seuil est 30km. Mais visualisation apres. On met donc un seuil plus grand pour flag auto
                        idKo = find(ranges > ECART_MAX_LOC_ARGOS*1000);
						% sauvegarde des resultats koba pour calcul des vitesses
                        T.position_qc_koba.data(idCyc(idSorted(o_idBadPos)))=6;
						
                        if  (~isempty(idKo))&~isempty(o_idBadPos) % if that test fail, use koba to flag bad positions => erreur locpos d'altran
                            
							fid_alerte=fopen(file_alerte,'a');
                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', LOC DISTANTES DE PLUS DE 30 km']);
                            fclose(fid_alerte);
                            fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', LOC DISTANTE DE PLUS DE 30 km']);
                            if Stat == 1
                                alertCyc_e9 = [alertCyc_e9 cycles_sorted(id)];
                            end
                            
                            % if T.position_qc.data(idCyc(o_idBadPos))~=4;
                                % T.position_qc.data(idCyc(o_idBadPos))=6;
                            % end
							if T.position_qc.data(idCyc(idSorted(o_idBadPos)))~=4;  % correction cc 07/01/2021
                                T.position_qc.data(idCyc(idSorted(o_idBadPos)))=6;
                            end
                        end
                        % [o_date, o_longitude, o_latitude, o_posAcc, o_posQcIn, o_posQcOut, ...
                        % o_speed,o_idBadPos, Little_Time, Little_Time_Sp, deltaTime_Sp, normal_Time_Sp] = check_argos_positions_koba_all_gh(LocDate_sorted,...
                        % locLon_sorted, locLat_sorted, locQc, locPosition_qc, ...
                        % locDate_prec_ok(end), dlonlas, dlatlas, fid_koba);
                        %if(o_posQcOut(find(o_speed==max(o_speed)))==1)
                        %if(isempty(find(o_speed>76.88 & o_speed<76.90))==0)
                        %keyboard
                        %end
                        
                        
                        % %%Pour l'affichage des vitesse brutes et flaguees
                        
                        % if(isempty(o_idBadPos)==0 & o_idBadPos~=length(o_date))
                        % o_idBad_speed = [o_idBadPos, o_idBadPos+1];
                        % else
                        % o_idBad_speed = o_idBadPos;
                        % end
                        % o_speed_isok = o_speed;
                        % if(isempty(o_idBad_speed)==0)
                        % o_speed_isok(o_idBad_speed)=NaN;
                        % end
                        % speed_Cyc =  [speed_Cyc, o_speed_isok];
                        % if(isempty(idDoublon_date)==0)
                        % ppp = ppp+length(idDoublon_date);
                        % Speed_doublon(ppp-(length(idDoublon_date)-1):ppp,:) = o_speed_isok(idDoublon_date);
                        % end
                        % o_speed_isok(idDoublon_date)=NaN;  %%car koba ne prend pas en compte les doublons de date si pas associés à doublons de lon/lat
                        % o_speed_isokd = o_speed_isok;
                        
                        if(isempty(o_idBadPos)==0)&isequal(LocDateok_sorted,locDateok)==0
						%keyboard
						end
                        if(isempty(o_idBadPos)==0)
						    
                            fid_alerte=fopen(file_alerte,'a');
                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', LOCA ARGOS ERRONEE (d''apres test de Koba), Position ' num2str(o_idBadPos') ' ']);
                            fclose(fid_alerte);
                            %fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', LOCA ARGOS ERRONEE (d''apres test de Koba) Position '  num2str(o_idBadPos')  ' ']);
                            if Stat == 1
                                alerte30(k,cycles_sorted(id)+1)=str2double(floatname);
                                
                            end
                            locPosition_qc(o_idBadPos)=6;   %%%flag à 6 les positions reperees comme erronees.
                            
                            if(Stat == 1 & isempty(find(o_idBadPos==1))==0)
                                % alertCyc_e5 = [alertCyc_e5 cycles_sorted(id)] ; % cc remove 26/10/2020
                            end
                            %T.position_qc.data(idCyc(o_idBadPos))=6; % %
                            %cc 26/10/2020: a voir ... plutot ne pas les
                            %utiliser dans calcul vitesse mais laisser les
                            %flags tels quels
                        end
                        
                        %end
                        %end   %%fin de la condition sur id>1
                    end  %%fin de la condition (sum(idCyc) > 0 & id==1) || (sum(idCyc)>0  & sum(idCycprec{id}) > 0)
                    % idCycprec=idCyc_sorted;
                    
                end   %%%Fin de la boucle sur les cycles
                %
                
                alertCyc_e3=alertCyc_e3_temp;
                % add cc 06/11/2020Elimination des alertes locterre des alertes grounded (pour comparaison avec altran seulement qui ne prend pas toujours en compte le grounded si le flotteur est a terre).
                % islocterre_grounded=ismember(alertCyc_e3_temp,alertCyc_e8);
                % alertCyc_e3=alertCyc_e3_temp(~islocterre_grounded); %
                
                alertCyc_e7=alertCyc_e7_temp;
                % add cc 07/11/2020 Rajout des alertes proche surface aux alertes  pression  (pour comparaison avec altran seulement)
                % alertCyc_e7=union(alertCyc_e7_temp,cycle_proche_surface);  %
                % if size(alertCyc_e7,1)>1
                % alertCyc_e7=alertCyc_e7';
                % end
                
                %%%Ecriture des alertes dans un fichier txt pour STATS
                
                if(Wr_Err==1)
                    if(isempty(alertCyc_e3)==0)
                        fid3=fopen(file3,'a');
                        %fprintf(fid3,'%s\n',[floatname ', [' num2str((alertCyc_e3)) ']']);   %%% Grounded  cc 18/12/2020 : elimination des doubles dans les alertes
						fprintf(fid3,'%s\n',[floatname ', [' num2str(unique(alertCyc_e3)) ']']);   %%% Grounded
                        fclose(fid3);
                    end
                    if(isempty(alertCyc_e4)==0)
                        fid4=fopen(file4,'a');
                        fprintf(fid4,'%s\n',[floatname ', [' num2str(unique(alertCyc_e4)) ']']);   %%% alerte cycle
                        fclose(fid4);
                    end
                    if(isempty(alertCyc_e5)==0)
                        fid5=fopen(file5,'a');
                        fprintf(fid5,'%s\n',[floatname ', [' num2str(unique(alertCyc_e5)) ']']);    %%% alerte locdate
                        fclose(fid5);
                    end
                    if(isempty(alertCyc_e6)==0)
                        fid6=fopen(file6,'a');
                        fprintf(fid6,'%s\n',[floatname ', [' num2str(unique(alertCyc_e6)) ']']);   %%% alerte launchdate
                        fclose(fid6);
                    end
                    if(isempty(alertCyc_e7)==0)
                        fid7=fopen(file7,'a');
                        fprintf(fid7,'%s\n',[floatname ', [' num2str(unique(alertCyc_e7)) ']']);   %%% alerte pression
                        fclose(fid7);
                    end
                    if(isempty(alertCyc_e8)==0)
                        fid8=fopen(file8,'a');
                        fprintf(fid8,'%s\n',[floatname ', [' num2str(unique(alertCyc_e8)) ']']);   %%% alerte loc terre
                        fclose(fid8);
                    end
                    if(isempty(alertCyc_e9)==0)
                        fid9=fopen(file9,'a');
                        fprintf(fid9,'%s\n',[floatname ', [' num2str(unique(alertCyc_e9)) ']']);   %%% alerte loc pos
                        fclose(fid9);
                    end
                    if(isempty(alertCyc_e10)==0)
                        fid10=fopen(file10,'a');
                        fprintf(fid10,'%s\n',[floatname ', [' num2str(unique(alertCyc_e10)) ']']);   %%% Incoherence meta/traj pour durée de cycles
                        fclose(fid10);
                    end
                end
                
                
                
                
                
                %   figure
                %subplot(3,1,1)
                %plot(1:1:i_rpp,pres_drift_mes,'r*',1:1:i_rpp,pres_drift_th,'.',cyc_wrong_pres,pres_wrong,'go',cyc_wrong_ground,ground_wrong,'yo');
                %title('Pression moyenne de parking')
                %xlabel('cycles')
                %ylabel('dbar ~ m')
                %legend('Pression mesurée','Profondeur estimée par la climatologie','Mesure aberrante','Grounded')
                %subplot(3,1,2)
                %plot(1:1:i_rpp,temp_drift_mes,'r*',1:1:i_rpp,temp_drift_th,'.',cyc_wrong_temp,temp_wrong,'go')
                %title('Température moyenne mesurée')
                %xlabelfopen('Alertes_cycle.txt','w+');('cycles')
                %ylabel('°C')
                %legend('Température mesurée','Température estimée par la climatologie','Mesure aberrante')
                %subplot(3,1,3)
                %plot(1:1:i_rpp,psal_drift_mes,'r*',1:1:i_rpp,psal_drift_th,'.',cyc_wrong_psal,psal_wrong,'go')
                %title('Salinité moyenne mesurée')
                %xlabel('cycles')
                %ylabel('-')
                %legend('Salinité mesurée','Salinité estimée par la climatologie','Mesure aberrante')
                
                % figure
                % plot(1:1:id,elev_all,'k',1:1:id,-pres_drift_mes2,'*r');
                % title('Pression mesuree (r) et bathy (k)');
                % ylabel('Pression');
                
                
                %figure
                %plot(distancederivesurf,'k-*');
                
                
                %figure
                % plot(1:1:id,temp_drift_th2,'k',1:1:id,temp_drift_mes2,'*r');
                %plot(temp_drift_mes,'k',temp_drift_th,'*r');
                %title('Temp mesuree(k) et theorique (r)');
                %ylabel('T');

                
                % Sauvegarde des nouveaux fichiers traj contenant les flags et les vitesses (deep, surface , errors,...)
                if P.save_traj_file==1
                    filename_new =[P.DIR_NEW_TRAJ_FILE traj_fileName_final];
                    if exist(P.DIR_NEW_TRAJ_FILE)==0
                        mkdir(P.DIR_NEW_TRAJ_FILE)
                    end
					
                    T =  create_new_variables(T,DimT,'deep_velocity'); % cree les champs vides pour les vitesses (deep, surface), les erreurs 
					T =  create_new_variables(T,DimT,'first_surface_velocity'); 
					T =  create_new_variables(T,DimT,'last_surface_velocity'); 
					
					% Compute velocities and write in netcdf file
					T= compute_velocities(T,unique([alertCyc_e4 alertCyc_e8] ));
					
					% Compute velocities and write in andro file
					compute_velocities_to_yo(T,P,unique([alertCyc_e4 alertCyc_e8]));
					
					% remove temporary variable
                    T=rmfield(T,'position_qc_koba');
					T=rmfield(T,'representative_park_temperature');
					% save new netcdf files 
                    create_netcdf_allthefile(T,DimT,filename_new,GlobT)
                    
                end
            end  % condition existence de loc
        end % condition existence fichier traj
        
    end  %% fin de la boucle sur les flotteurs
    
    
    disp('*****FIN DES TESTS******')
    
    % fclose(fid3)
    % fclose(fid4);fclose(fid5); fclose(fid6);fclose(fid7); fclose(fid8);fclose(fid9);
    % fclose(fid10);
    
    if Stat==1&isempty(idLoc)==0
        Diff_Cyc = Diff_Medcyc_Mcyc(Diff_Medcyc_Mcyc>0);
        Name_fileStat = [DIR_STAT P.liste_rep_alerte{ilist} '/variables_stat.mat'];   % modif chemin sauvegarde cc 02/11/2020
        save(Name_fileStat,'Diff_Medcyc_Mcyc','Diff_Ppark','Diff_Cyc',...
            'Diff_LaunchD','Diff_LaunchLat','Diff_LaunchLon','ecartMaxLocCycle_fl',...
            'ecartMeanLocCycle_fl','ecartFirstLastLoc','ecartLaunchDateFirstLoc',...
            'ecartLaunchDateFirstLoc','distancederiveprof',...
            'distancederivesurf','pres_drift_mes','temp_drift_mes','psal_drift_mes');
        
        Name_fileStat2 = [DIR_STAT P.liste_rep_alerte{ilist} '/alertes_stat.mat'];     % modif chemin sauvegarde cc 02/11/2020
        save(Name_fileStat2,'alerte1', 'alerte2', 'alerte3', 'alerte4', 'alerte5', 'alerte6', ...
            'alerte7', 'alerte8', 'alerte9', 'alerte10', 'alerte11', 'alerte12', 'alerte13', 'alerte14', ...
            'alerte15', 'alerte16', 'alerte17', 'alerte18','alerte19', 'alerte20', 'alerte21', ...
            'alerte21', 'alerte22', 'alerte23', 'alerte24', 'alerte25', 'alerte26', 'alerte27', ...
            'alerte28', 'alerte29', 'alerte30');
    end
    
end  % boucle sur ilist



