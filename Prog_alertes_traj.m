
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

% Liste des flotteurs a� tester
fid1 = fopen(Liste_Float);
allfloats=textscan(fid1,'%s\n','Delimiter',',','CommentStyle','#');
fclose(fid1);

%% Creation du r�pertoire /STAT s'il n'existe pas
if ~exist(DIR_STAT,'dir')
  mkdir(DIR_STAT);
end

%% Fichier log des alertes
global fid_alerte;
fid_alerte = fopen([DIR_HOME,'RT_alertes.txt'],'wt');    

%% Fichiers log des alertes par erreurs
if(Bathy==1)
   fid3=fopen('Alertes_grounded_etopo.txt','w+');
else
   fid3=fopen('Alertes_grounded_srtm_b.txt','w+');
end
fid4 = fopen('Alertes_cycle_b.txt','w+');
fid5= fopen('Alertes_locdate_b.txt','w+');
fid6 = fopen('Alertes_launchdate_b.txt','w+');
fid7 = fopen('Alertes_pressure_b.txt','w+');
fid8 = fopen('Alertes_locterre_b.txt','w+');
fid9 = fopen('Alertes_locpos_b.txt','w+');
fid10 = fopen('Alertes_metatrajcyc_b.txt','w+');


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

% Lecture des parametres temperature et salinite de la climatologie ISAS13
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

%%boucle sur les flotteurs

for k=1:length(allfloats{1})
    

    alertCyc_e1 = []; alertCyc_e2 = []; alertCyc_e3 = []; alertCyc_e4 = []; alertCyc_e5 = [];
    alertCyc_e6 = [];alertCyc_e7 = [];alertCyc_e8 = [];alertCyc_e9 = []; alertCyc_e10=[];


   
    floatname = strtrim(allfloats{1}{k});
    doneT=0;
    doneM=0;
    ii=0;
    elev_all = [];
    while ~doneT&&~doneM&& ii<=length(DAC)
        % lecture des traj
        ii=ii+1;
        traj_fileName_R = [DIR_FTP DAC{ii} '/' floatname '/' floatname '_Rtraj.nc'];
        traj_fileName_D = [DIR_FTP DAC{ii} '/' floatname '/' floatname '_Dtraj.nc'];
        if exist(traj_fileName_D,'file');
            T=read_netcdf_allthefile(traj_fileName_D);
            doneT=1;
        elseif exist(traj_fileName_R,'file');
            T=read_netcdf_allthefile(traj_fileName_R);
            doneT=1;
        end
        meta_fileName=[DIR_FTP DAC{ii} '/' floatname '/' floatname '_meta.nc'];
      
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
        icounterfloat=icounterfloat+1;
        T_sav=T;
        T=format_flags_char2num(T);%change flag char strings to numerical vectors
                                   %ex '11111441111 ' -> [1 1 1 1 1 4 4 1 1 1 1 999]

        T=replace_fill_bynan(T);   % remplace les fillValues (9999.. par des NaN)
        
        % VALEURS DES LOCALISATIONS ARGOS DE SURFACE
        if isfield(T,'measurement_code')==0
            idLoc =   find(~isnan(T.longitude.data)& ~isnan(T.latitude.data));
            fprintf(fid_alerte,'%s\n',[ floatname ' Attention pas de measurement_code dans le fichier Traj.'])
            fprintf('%s\n',[ floatname ' Attention pas de measurement_code dans le fichier Traj.'])
        else
            idLoc =   find(T.measurement_code.data==703);     %%%ne prend pas la première loca en surface car launch 
        end 
        
        % VALEURS DES MESURES DE PARKING
        if isempty(find(T.measurement_code.data==290, 1)) & isempty(find(T.measurement_code.data==296))==0
            idDrift=find(T.measurement_code.data==296);% prend la moyenne des mesures de parking
        elseif isempty(find(T.measurement_code.data==290))==0
            idDrift=find(T.measurement_code.data==290);% prend l'ensemble des mesures de parking si existe
        else
            idDrift=find(T.measurement_code.data==299);
            if(isempty(idDrift)==0)
              fprintf(fid_alerte,'%s\n',[floatname ' PAS DE MEASUREMENT_CODE 290 ou 296 (MESURES DE PARKING), RECUPERE ID POUR MEASUREMENT_CODE 299.']);
              fprintf('%s\n',[floatname 'PAS DE MEASUREMENT_CODE 290 ou 296 (MESURES DE PARKING), RECUPERE ID POUR MEASUREMENT_CODE 299.']);
            end
        end
        
         if(isempty(idDrift)==1)
            fprintf(fid_alerte,'%s\n',[floatname ' !!!PAS DE MEASUREMENT_CODE 290 ou 296 (MESURES DE PARKING) OU 299.']);
            fprintf('%s\n',[floatname '!!!PAS DE MEASUREMENT_CODE 290 ou 296 (MESURES DE PARKING), OU 299.']);
         end  
        
        % CREATION DE LA LISTE DES NUMEROS DE CYCLES EXISTANTS
        CycLoc=T.cycle_number.data(idLoc);
        numCycleMaxTraj = max(CycLoc);
        numCycleMinTraj = min(CycLoc(CycLoc>=0));% cycle -1 dans les traj: date et position de mise à l'eau
        
        cycles_1 = [numCycleMinTraj:numCycleMaxTraj];% creation de la liste theorique des cycles (incrementation de 1)
        cycles = unique(CycLoc(CycLoc>=0))';% creation de la liste relle avec sauts de cycles
        
        clear CycLoc ;
        
        % CALCUL DE LA MOYENNE DES MESURES
        dureeMedianCycle=[];
        end_cycles=[];
        pres_ok = [];
        cycles_m={};   
           
        % on selectionne les missions renseignees 
        missions=unique(T.config_mission_number.data,'stable')';  
        isnot = find(ismember(M.config_mission_number,unique(T.config_mission_number.data))==0); 
        if(isempty(isnot)==0) M.CycleTime(isnot)=NaN;end
%%pour les cas o� valeur NaN pour un numero de mission donne dans T.config_mission_number et non dans M.CycleTime
                            
        for cm=1:length(cycles)     %%%pour tous les cycles
                   
            idCyc_m_d=idLoc(T.cycle_number.data(idLoc)==cycles(cm));
            idCyc_m_p=idDrift(T.cycle_number.data(idDrift)==cycles(cm));
            isok_d=~isnan(T.juld.data(idCyc_m_d));
            isok_p=~isnan(T.pres.data(idCyc_m_p));
            date_ok=T.juld.data(idCyc_m_d(isok_d));
            end_cycles(cm)=date_ok(end);       
        
        end
        clear idCyc_m_d idCyc_m_p isok_d isok_p;
        
        duree_cycle = diff(end_cycles)./diff(cycles);    %%duree de chaque cycle  (sauf le cycle 0 ou le cycle 1 pour flotteur APEX)
       
        
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
            
                if(~isempty(index)) cyclesm(index)=NaN; cycles_m{m} = cyclesm; end    %%ne regarde qu'a� partir du cycle 1.
            
                if(length(find(~isnan(cycles_m{m})))>=1)    %% pour les cas ou� cycles_m ~= juste le cycle 0
                     Begin = cycles(1);
                     if Begin == 0
                          duree_cycle_m = duree_cycle(find(ismember(cycles,cycles_m{m})==1)-1);
                     else                                   %%%pour le cas des fotteurs où cycles(1) = 1 comme apex ou cycles(1)=2,3...etc
                         cycles_m_Withoutfirst = cycles_m{m};
                         cycles_m_Withoutfirst = cycles_m_Withoutfirst(cycles_m_Withoutfirst>Begin);
                         duree_cycle_m = duree_cycle(find(ismember(cycles,cycles_m_Withoutfirst)==1)-1);
                         if cycles_m{m}==1; duree_cycle_m = NaN;end
                     end
        
                     end_cycles_m = end_cycles(find(ismember(cycles,cycles_m{m})==1));
            
                    %% calcul de la mediane des durees de cycles 
                   
                     if length(duree_cycle_m)==1
                          dureeMedianCycle(m)=duree_cycle_m;    
                     elseif length(duree_cycle_m)~=1 & mynanstd(duree_cycle_m)<100 
                          dureeMedianCycle(m)=median(duree_cycle_m); 
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
              
                     if abs(dureeMedianCycle(m)-M.CycleTime(missions(m)))>0.1        
                                                
                        fprintf(fid_alerte,'%s\n',[ floatname ', INCOHERENCE DUREE CYCLE ENTRE META (', num2str(M.CycleTime(m)),' j) ET TRAJ (',num2str(dureeMedianCycle(m)), ' j)']); 
                        fprintf('%s\n',[ floatname ', INCOHERENCE DUREE CYCLE ENTRE META (',num2str(M.CycleTime(m)),' j) ET TRAJ (',num2str(dureeMedianCycle(m)), ' j)']);
                        if(Stat==1)
                            alertCyc_e10 = [alertCyc_e10 cycles_m{m}];
                            alerte1(k)=str2double(floatname);
                        end
                     end
                     
                     %%%pour STATS
                     ite = ite+1;
                     Diff_Medcyc_Mcyc(k,m) = dureeMedianCycle(m)-M.CycleTime(missions(m));
                     
                     idLoc(T.cycle_number.data(idLoc)==cycles(cm))         
               else
                   Diff_Medcyc_Mcyc(k,m)=NaN;  
                      
               end
           
           else
              dureeMedianCycle(m)=NaN;
              Diff_Medcyc_Mcyc(k,m)=NaN;
           end      %%fin de la boucle sur la condition isnan
                                  
        end   %%% fin de la boucle sur les missions
        
         if length(missions)>length(M.CycleTime)
               
                fprintf(fid_alerte,'%s\n',[ floatname ', INCOHERENCE META & TRAJ: Nbre de N� MISSION DANS META (', num2str(length(M.CycleTime)),' ) < � celui dans TRAJ (', num2str(length(missions)), ' ). Ne peut pas correctement verifier la coherence ni pour la duree des cycles ni pour P_park.']); 
                fprintf('%s\n',[ floatname ', INCOHERENCE META & TRAJ: Nbre de N� MISSION DANS META (',num2str(length(M.CycleTime)),' ) < � celui dans TRAJ (', num2str(length(missions)), ' ). Ne peut pas correctement verifier la coherence ni pour la duree des cycles ni pour P_park.']);
                if (Stat==1)
                    alertCyc_e10 = [alertCyc_e10 999];
                end
         end  
        
            
        
        % -------------------------------------------------------------
        % alerte 3, 4, 5, 6
        % BOUCLE SUR CHAQUE CYCLE POUR VERIFIER:
        % - double de cycles
        % - pb croissance cycle
        % - mauvaise duree du cycle
        
        [cycles_sorted,idCyc_sorted]= sort(cycles);
        
       [o_alerte3, o_alerte4, o_alerte5,...
           o_alerte6,o_alertCyc_e4, idCycprec]= Test_cycles(idLoc, cycles, cycles_1,cycles_sorted,duree_cycle, ...
                          missions, dureeMedianCycle);
        if(isempty(o_alerte3)==0)  alerte3{k} = o_alerte3;end;
        if(isempty(o_alerte4)==0)  alerte4{k}  = o_alerte4; end;
        if(isempty(o_alerte5)==0)  alerte5{k} = o_alerte5;end;
        if(isempty(o_alerte6)==0) alerte6{k} = o_alerte6; end;
        
        if(isempty(o_alertCyc_e4)==0)  alertCyc_e4 = [alertCyc_e4 o_alertCyc_e4]; end;
       
       
             
       %%---------------------------------------------------------------
       % alerte 7, 8,9
       % VERIFICATION DES DATES DE LOC ARGOS: DATES REALISTES: later than 1 st January 1997
       % VERIFICATION DES LOC ARGOS: 
       % VERIFICATION DES DATES ET LOCS DE LANCEMENT (coherence entre les fichiers de trajectoires et de metadonnees)
       
       iLaunch=find(T.cycle_number.data==-1)
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
        
        % ---------------------------------------------------------------     
        
        % BOUCLE SUR CHAQUE CYCLE POUR VERIFIER:
        % - coherence de la climatologie (temperature, salinite)
        % - calcul pression mediane de derive
        % - coherence climatologie pression de derive.
        %% pres_ok: moyenne de la pression a la profondeur de parking par cycle
        %% elev_all: bathy la premi�re position de chaque cycle
        %% elev_end_all: bathy pour la derni�re position du cycle pr�c�dent
        
        
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
        temp_alert=0;
        temp_non_ref=0;
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

           temp_drift_mes=[];
           temp_drift_th=[];
           psal_drift_mes=[];
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
          
           
        

           istat=istat+1;
           % numero du cycle
           numCycle = cycles_sorted(id);
           numMis=T.config_mission_number.data(id);
           % id des mesures associees  a� ce cycle
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
             
           if ~isempty(idCyc_drift)&~isnan(T.pres.data(idCyc_drift))
            
             i_rpp=i_rpp+1;
             pres_long=locLon(1);
             pres_lat=locLat(1);
             lonprec = T.longitude.data(idCycprec{id});
             pres_lonend= lonprec(max(find(~isnan(lonprec))));
             latprec = T.latitude.data(idCycprec{id});
             pres_latend = latprec(max(find(~isnan(latprec)))) ;   
             pres_drift_mes(i_rpp)=abs(mynanmean(T.pres.data(idCyc_drift)));       
             
             if(P.Bathy==1)
                   
                ilong = round(mean(find(LONG(:,1)<=pres_long+max(diff(LONG)/2) & ...
                LONG(:,1)>=pres_long-max(diff(LONG)/2))));               
                if(pres_long>max(LONG(:,1)))
                    ilong = find(LONG(:,1) == max(LONG(:,1)));   %%pour eviter ilong=NaN;
                end
                   
                ilat = round(mean(find(LAT(:,1)<=pres_lat+max(diff(LAT)/2) & ...
                LAT(:,1)>=pres_lat-max(diff(LAT)/2)))) ;                 
                    
                if(pres_lat>max(LAT(:,1)))
                    ilong = find(LAT(:,1) == max(LAT(:,1)));       %%%pour eviter ilat=NaN;
                end
                   
                elev = ELEV(ilat,ilong);     %%%altitude de la premiere pos du cycle n
                   
                   
                if(id>1)
                    ilong_end = round(mean(find(LONG(:,1)<=pres_lonend+max(diff(LONG)/2) ...
                     & LONG(:,1)>=pres_lonend-max(diff(LONG)/2))));
                     if(pres_lonend>max(LONG(:,1)))
                        ilong_end = find(LONG(:,1) == max(LONG(:,1)));   %%pour eviter ilong=NaN;
                     end
                     ilat_end = round(mean(find(LAT(:,1)<=pres_latend+max(diff(LAT)/2) ...
                      & LAT(:,1)>=pres_latend-max(diff(LAT)/2))));
                     if(pres_latend>max(LAT(:,1)))
                         ilat_end = find(LAT(:,1) == max(LAT(:,1)));       %%%pour eviter ilat=NaN;
                     end
                     elev_end = ELEV(ilat_end,ilong_end);    %%%altitude de la derniere pos du cycle precedent.
                end                   
                   
             else
                   %%%recupere la bathy de la premiere pos du cycle
                  
                [o_elev, o_lon, o_lat] = get_srtm_elev(pres_long, pres_long,...
                 pres_lat, pres_lat);
                   
                 if(length(o_lon)>1 & length(o_lat)>1 & abs(pres_long)<max(abs(o_lon)) ...
                  & abs(pres_long)>min(abs(o_lon)))
                      elev = interp2(o_lon,o_lat,o_elev,pres_long,pres_lat);
                 elseif(length(o_lon)<=1 || length(o_lat)<=1)
                     elev = mean(o_elev);
                 else
                     elev = mean(mean(o_elev));
                 end   
                  
                   
                 if(id>1 & isempty(idCycprec{id})==0)                          
                  %%%recupere la bathy de la derniere pos du cycle
                  %%%precedent s'il existe
                      [o_elev_end, o_lon_end, o_lat_end] = get_srtm_elev(pres_lonend, pres_lonend,...
                      pres_latend, pres_latend);    
                      if(length(o_lon_end)>1 & length(o_lat_end)>1 & abs(pres_lonend)<max(abs(o_lon_end)) ...
                       & abs(pres_lonend)>min(abs(o_lon_end)))
                           elev_end = interp2(o_lon_end,o_lat_end,o_elev_end,pres_lonend,pres_latend);
                      elseif(length(o_lon_end)<=1 || length(o_lat_end)<=1)
                           elev_end = mean(o_elev_end);
                      else
                           elev_end = mean(mean(o_elev_end));
                      end   
                 else
                      elev_end = NaN;
                                  
                 end
             end
              
             %%%recuperation de l'altitude par i_rpp pour pouvoir la
             %%%reutiliser dans une autre boucle sur les cycles.
             elev_all(i_rpp) = elev;
             if exist ('elev_end','var') == 1
                elev_end_all(i_rpp) = elev_end;    
             end
                 

  
             %On determine les indices ISAS en longitude, latitude, pression
             if max(TEMP_LONG)>=pres_long% probleme a la limite du 180°
                ilong_drift=round(mean(find(TEMP_LONG(:,1)<=pres_long+max(diff(TEMP_LONG)/2) ...
                & TEMP_LONG(:,1)>=pres_long-max(diff(TEMP_LONG)/2))));
             else
                ilong_drift=1;
             end
             ilat_drift = round(mean(find(TEMP_LAT(:,1)<=pres_lat+max(diff(TEMP_LAT)/2) ...
               & TEMP_LAT(:,1)>=pres_lat-max(diff(TEMP_LAT)/2))));

             if max(TEMP_DEPTH)>=pres_drift_mes(i_rpp)
                 idepth_drift=round(mean(find(TEMP_DEPTH(:,1)<=abs(pres_drift_mes(i_rpp))+max(diff(TEMP_DEPTH)/2) ...       % estimation de l'indice : a recalculer 
               & TEMP_DEPTH(:,1)>=abs(pres_drift_mes(i_rpp))-max(diff(TEMP_DEPTH)/2))));
             else
                 idepth_drift=length(TEMP_DEPTH(:,1));
             end

             if max(TEMP_DEPTH_STD)>=pres_drift_mes(i_rpp)% certaines profondeurs extremes ne sont pas referencees
                 idepth_std_drift=round(mean(find(TEMP_DEPTH_STD(:,1)<=abs(pres_drift_mes(i_rpp))+max(diff(TEMP_DEPTH_STD)/2) ... 
                 & TEMP_DEPTH_STD(:,1)>=abs(pres_drift_mes(i_rpp))-max(diff(TEMP_DEPTH_STD)/2))));
             else
                 idepth_std_drift=length(TEMP_DEPTH_STD(:,1));
             end



             %% alerte 11, 12
             % VERIFICATION DE L'ACORD AVEC LA CLIMATOLOGIE POUR T et S
             % ---------------------------------------------------------
             % on veut verifier que les mesures des parametres physiques
             % sont coherentes en comparant avec le modele donne (ici ISAS)
              
             [o_alerte11, o_alerte12, temp_alert, temp_non_ref] = Test_TS(idCyc_drift, idCycprec, id, cycles_sorted,pres_long,pres_lat,numMis, missions,i_rpp, idepth_drift,...
             ilong_drift, ilat_drift, idepth_std_drift);
                  
             
             if(isempty(o_alerte11)==0)  alerte11{k} = o_alerte11;end;
             if(isempty(o_alerte12)==0)  alerte12{k}  = o_alerte12; end;
                
                                                 
                             
              %%% Met le qc à 4 pour les pressions < -2800 ou < bathy  
      
              pres_drift_th(i_rpp)=-TEMP_DEPTH(idepth_drift);
              pres_drift_mes(i_rpp)=-pres_drift_mes(i_rpp); 
             
              if(pres_drift_mes(i_rpp)<=elev-100 || pres_drift_mes(i_rpp) <-2200)
                 T.pres_qc.data(idCyc_drift)=4;    
              end
                
              qcCyc = T.pres_qc.data(idCyc_drift);
              isok_p = (qcCyc~=4 & qcCyc~=6 & ~isnan(qcCyc));  %%prend en compte le nouveau qc
              selec = idCyc_drift(isok_p); select = T.pres.data(selec);   
                
              %%recuperation des pressions de derive et moyenne par cycle
              if(elev<-1500)   %%%ne recupere les pressions que quand la bathy est <-1500
                 
                  pres_ok(id) = abs(mynanmean(select(find(select>50))));     %%%ne prend pas les pressions <50m 
              else
                  pres_ok(id) = NaN;
              end
               
           else
               pres_ok(id)=NaN;   %%%si le cycle n'a pas de idDrift
           end
            
       
       end   %%% fin de la boucle sur les cycles
        
                             
        %------------------------------------------------------------------
        
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
        
       
                      
        %%%%% CALCUL DE LA VALEUR MEDIANE DE LA PRESSION PAR n� DE MISSION EN PRENANT EN
        %%%%%  COMPTE LES NOUVEAUX QC  ET TEST SUR LES INCOHERENCES
        %%%%%  TRAJ/META POUR LA PRESSION.
        
        presMedianDrift=[];
        
        if ~isempty(idDrift)& isempty(find(~isnan(T.pres.data(idDrift))))==0  %% si dispose de mesure en profondeur et si au moins une valeur nonNan parmi les pressions
             
            %%par n� de mission
             for m = 1:length(missions)
          
                 idcMis=find(T.config_mission_number.data==missions(m));
                 Tcycun  = unique(T.cycle_number.data);Tcycun = Tcycun(Tcycun>=0);  %%% ne prend pas les cycles -1
                 idcMisok  = find(ismember(cycles,Tcycun(idcMis))==1);                %% prend en compte les cas où le numéro de cycle est dans T.config_number mais pas dans cycles car pas de MC=703
            
                 if(~isempty(idcMisok))      %%% si pas de cycle attribué au n� de mission (peut arriver quand cyle 0 mais pas de MC 703)
                    pres_ok_m{m} = pres_ok(idcMisok);
                 else
                    pres_ok_m{m} = NaN;
                 end
            
                 pres = pres_ok_m{m}; 
                  
                 %% Calcul de la pression mediane par n� de mission
                 if(~isempty(idcMisok))
                      if(length(find(~isnan(cycles_m{m})))>=1)    %% pour les cas ou cycles_m{m} ~= juste le cycle 0
                   %%prise en compte du qc mis suite a� la verif avec la bathy
                         if mynanstd(pres(pres>200&pres<2000))<200  
                     %presMedianDrift(m)=-abs(median(pres_ok_m(pres_ok_m>500&pres_ok_m<1500),2));    %%qu'en est-il quand Ppark ~ 300m ?
                            presMedianDrift(m)=-abs(median(pres(pres>200&pres<2000),2));
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
                               
                     if (abs(presMedianDrift(m)- (-M.ParkPressure(missions(m))))>20) % NC laisser une marge d'erreur
                        fprintf(fid_alerte,'%s\n',[ floatname ', INCOHERENCE ENTRE META ET TRAJ POUR P_park : |TPmed-MP| >20 db']);
                        fprintf('%s\n',[ floatname ', INCOHERENCE ENTRE META ET TRAJ POUR P_park : |TPmed-MP| > 20 db']);
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
                
                %%%recuperation de la profondeur de parking a� partir de
                %%%l'info dans les meta si elle existe. sinon prend la
                %%%valeur mediane calculee pour chaque mission. On prend
                %%%aussi cette derniere dans le cas ou la valeur de la
                %%%pression mediane est bien differente de la valeur indiquee dans
                %%%les meta.    
                
                if m<=length(M.ParkPressure) & abs(presMedianDrift(m) - -M.ParkPressure(m)) <  30    %%si nbre de mission dans meta et traj est <= et si la mediane de la pression est proche de celle indiquee dans les meta
                    park_prof(m) = -M.ParkPressure(m);
                elseif m>length(M.ParkPressure) || presMedianDrift(m)<-M.ParkPressure(m)-30 ...  %%%si nbre de mission dans meta manquant ou si la mediane de la pression trop eloignee de celle indiquee dans les meta
                    || presMedianDrift(m)>-M.ParkPressure(m)+30
                
                    park_prof(m) = presMedianDrift(m);
                end
                
             end  %%fin de la boucle sur les missions

        end   %%fin de la boucle sur idDrift
        
       %%% alerte 10, 13, 14
       %%%%BOUCLE SUR LES CYCLES
       
       i_rpp=0;
       istat=0;
       ipres=0;
      
       for id=1:length(cycles_sorted)
             
            istat=istat+1;
            % numero du cycle
            numCycle = cycles_sorted(id);
            numMis=T.config_mission_number.data(id);
            idMis = find(double(missions)==numMis);
            % id des mesures associees  à ce cycle (et au cycle suivant
            % pour utilisation ult�rieure
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
             
            if ~isempty(idCyc_drift)&~isnan(T.pres.data(idCyc_drift))
                i_rpp=i_rpp+1;  

                %pres_drift_mes(i_rpp)=abs(mynanmean(T.pres.data(idCyc_drift)));       

                % Verification de la pression de derive  (par rapport à la
                % pression mediane)     (remplacer par park_prof ?)
                %keyboard     
                if (pres_drift_mes(i_rpp)<presMedianDrift(idMis)-100 ... 
                   | pres_drift_mes(i_rpp)>presMedianDrift(idMis)+100)        %| presMedianDrift(idMis)<=elev_all(i_rpp)-50   %%%< elev-100 plutôt ? ou enlever cette condition car vérif grounded plus loin
                  
                    ipres=ipres+1;
                    % recuperation des donnees aberrantes
                    cyc_wrong_pres(ipres)=i_rpp;
                    pres_wrong(ipres)=pres_drift_mes(i_rpp);
                    fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB PRESSION / PROFONDEUR PARKING ( ' num2str(pres_drift_mes(i_rpp)) 'm) ~= Mediane: ' num2str(presMedianDrift(idMis)) 'm.']);
                    fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB PRESSION / PROFONDEUR PARKING (' num2str(pres_drift_mes(i_rpp)) 'm) ~= Mediane: ' num2str(presMedianDrift(idMis)) 'm.']);% NC test très large, teste l'ecart à la pression de derive mediane 
                    pres_alert(i_rpp)=1;
                    if Stat==1
                       alertCyc_e7 = [alertCyc_e7 cycles_sorted(id)];
                       alerte10(k,cycles_sorted(id)+1)=str2double(floatname);
                    end
                else
                    pres_alert(i_rpp)=0;
                    it = it+1;  
                    if Stat==1 & isempty(idMis)==0
                        Diff_Pm_Pmed(it) =  pres_drift_mes(i_rpp)-presMedianDrift(idMis);
                    end
                end
                     
                if pres_alert(i_rpp)==1 & temp_alert(i_rpp)==1
                    fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB CAPTEUR PRESSION ? CAR PB PRESSION ET TEMPERATURE']);
                    fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB CAPTEUR PRESSION ?CAR PB PRESSION ET TEMPERATURE']); % NC verifier si on doit garder cette alerte sous cette forme
                    if Stat==1
                        alertCyc_e7 = [alertCyc_e7 cycles_sorted(id)];
                        alerte13(k,id)=str2double(floatname);
                    end
                    
                 end 

                 

                 if  (elev_end_all(i_rpp)>=park_prof(idMis)+30 | elev_all(i_rpp)>=park_prof(idMis)+30) ...
                     | (temp_non_ref(i_rpp)==1 & pres_alert(i_rpp)==1) 
                     alerte_grounded = 1;
                     iground=iground+1;
                     ifloat = ifloat+1;
                     % recuperation des donnees aberrantes
                     cyc_wrong_ground(iground)=i_rpp;
                     ground_wrong(iground)=pres_drift_mes(i_rpp);
                     %float_groundedgrounded(ifloat) = str2double(floatname);
                     %CyclesGrounded(iground) = cycles_sorted(id);
                     fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', GROUNDED ?']);
                     fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', GROUNDED ?']);                      %%% appliquer flags ?
                     if Stat==1
                        alertCyc_e3 = [alertCyc_e3 cycles_sorted(id)];
                        alerte14(k,cycles_sorted(id)+1)=str2double(floatname);
                     end
                     %fprintf(fid3,'%s\n',[floatname ', [' num2str(cycles_sorted(id)) ']']);
                       
                     if(GroundedNum(id)==2)    %%% 2 = 'N'
                         fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', alerte grounded mais incohérence avec la variable T.grounded.']);
                         fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', alerte grounded mais incohérence avec la variable T.grounded.']);                      %%% appliquer flags ?
                        %%%mettre la variable grounded à 'Y' ??
                           %T.grounded.data(id) = 'Y';
                     end
                        
                     if(map==1)
                           m_plot(TEMP_LONG(ilong_drift,1),TEMP_LAT(ilat_drift,1),'color',[0.4660 0.6740 0.1880],'marker','p','markersize',8)
                     end
                elseif(GroundedNum(id)==1) & alerte_grounded==0    %%%1 = 'Y'
                     %alertCyc_e3 = [alertCyc_e3 cycles_sorted(id)];
                     fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', pas d''alerte grounded alors que oui d''apres la variable T.grounded.']);
                     fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', pas d''alerte grounded alors que oui d''apres la variable T.grounded.']);                      %%% appliquer flags ?
                    %T.grounded.data(id)='P';
                elseif(GroundedNum(id)==3) & alerte_grounded==0      %%%3 = 'U'
                     fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', pas d''alerte grounded et ''U'' d''apres la variable T.grounded.']);
                     fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', pas d''alerte grounded et ''U'' d''apres la variable T.grounded.']);                      %%% appliquer flags ?
                     %T.grounded.data(id)='N';          
                end
                        
                
             end   %% fin de la condition sur l'existence des idDrift
                
             if(isnan(T.config_mission_number.data(id))==0)           
                 %%id mission du flotteur pour ce cycle
                 numMis=T.config_mission_number.data(id);   %%%attention pour id=1 peut etre = NaN;
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
                 ecartcycle=10;
             elseif strfind('ARGOS',M.trans_system)
                 bornesurf=100000;
                 ecartcycle=24;
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
              
                    [isdouble_toremove,idDoublon_date,idDoublon_latlon] = ...
                    find_doubles_datelatlon_30s(locDate, locLon, locLat, locQc, locTemp);   %%%recherche les doublons. Attention, il peut y avoir un double de date pas forcement associe à double de lon ou de lat
              
               
                    if isempty(isdouble_toremove)==0
                       fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DOUBLE DATES DE LOC']);
                       fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DOUBLE DATES DE LOC']);
                       T.position_qc.data(idCyc(idDoublon_date(find(T.position_qc.data(idCyc(idDoublon_date)))>1)))=6;   %%flag à 6 les dates doublées dont le qc est >1
                       if Stat==1
                         alertCyc_e9 = [alertCyc_e9 cycles_sorted(id)];
                         alerte15(k,cycles_sorted(id)+1)=str2double(floatname);
                       end
                    
                       if isequal(idDoublon_date,idDoublon_latlon)
                          fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ' id LOC ' num2str(idDoublon_latlon(1)) '&' num2str(idDoublon_latlon(2))  ', DOUBLE DATES DE LOC ASSOCIE A DOUBLES LON-LAT']);
                          fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ' id LOC ' num2str(idDoublon_latlon(1)) '&' num2str(idDoublon_latlon(2)) ', DOUBLE DATES DE LOC ASSOCIE A DOUBLES LON-LAT']);
                          T.position_qc.data(idCyc(idDoublon_latlon(find(T.position_qc.data(idCyc(idDoublon_latlon)))>1)))=6;   %%flag à 6 les lon-lat doublés dont le qc est >1
                          if Stat==1
                              alertCyc_e9 = [alertCyc_e9 cycles_sorted(id)];
                          end
                       end
                                
                    end
             
                     % QUE FAIRE avec des doubles (FLAG ?)
                    %T.juld_qc.data(isdouble_toremove)=6;        
                     T.juld_qc.data(idCyc(isdouble_toremove))=6;
               
                     locPosition_qc = T.position_qc.data(idCyc);     %%%de nouveau car les flags ont change
                     locDate_qc = T.juld_qc.data(idCyc);             %%%idem
                
                    % 
                    % VERIF CROISSANCE DES DATES DE LOC (qui ont un flag bon)
                    %--------------------------------------------------------
                    isbad = T.juld_qc.data(idCyc)==6|T.juld_qc.data(idCyc)==4;    
                    locDateok=locDate(~isbad);
                   [LocDateok_sorted, idSorted] = sort(locDateok);
                   if isequal(LocDateok_sorted,locDateok)==0
                      fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DATES DE LOC NON CROISSANTES']);
                      fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DATES DE LOC NON CROISSANTES']);
                      if Stat==1
                         alertCyc_e5 = [alertCyc_e5 cycles_sorted(id)];
                         alerte16(k,cycles_sorted(id)+1)=str2double(floatname);
                      end
                   end
                   [LocDate_sorted, idSorted] = sort(locDate);
                   locLon_sorted=locLon(idSorted);
                   locLat_sorted=locLat(idSorted);
                   locQc_sorted=locQc(idSorted);
                   locDate_qc_sorted=locDate_qc(idSorted);
                   locPosition_qc_sorted=locPosition_qc(idSorted);
                   idCyc_sorted=idCyc(idSorted);
                
               
                   % 
                   % CALCUL DE L'ECART MAX (en temps) ENTRE DEUX LOCS CONSECUTIVES ARGOS POUR UN MEME
                   % CYCLE (detecter des cycles qui ont ete accoles par alerte, comme
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
                
                   if ecartMaxLocCycle(istat)>10
                       fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ECART MAX 2 LOC > 10h :' num2str(ecartMaxLocCycle(istat)) ' heures']);
                       fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ECART MAX 2 LOC > 10h :' num2str(ecartMaxLocCycle(istat)) ' heures']);
                       if Stat==1
                           alertCyc_e5 = [alertCyc_e5 cycles_sorted(id)];
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
                   if (locDate_qc_sorted~=6|locDate_qc_sorted~=4);
                       ecartFirstLastLoc(istat)= LocDate_sorted(end)-LocDate_sorted(1);
                       if ecartFirstLastLoc(istat)*24 > ecartcycle & id<length(cycles_sorted) % on est pas sur le dernier cycle qui peut etre EOF
                           fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ECART FIRST LAST LOC >' num2str(ecartcycle) 'h :' num2str(ecartFirstLastLoc(istat)) ' jours']);
                           fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ECART FIRST LAST LOC >' num2str(ecartcycle) 'h :' num2str(ecartFirstLastLoc(istat)) ' jours']);
                           if Stat==1
                             alertCyc_e5 = [alertCyc_e5 cycles_sorted(id)];
                             alerte18(k,id)=str2double(floatname);
                           end
                      elseif ecartFirstLastLoc(istat)*24 > 60 & id<length(cycles_sorted)
                           fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',? EOL ?']);
                           fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',? EOL ?']); % peut-on considerer EOL à partir de 5 j à la surface ?
                           fleet_EOF(k,id)=str2double(floatname);
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
                
               
                  if length(LocDate_sorted)>1&(locDate_qc_sorted~=6|locDate_qc_sorted~=4);
                   
                    if (numCycleMinTraj==0|numCycleMinTraj==1)&id==1&str2double(M.launch_qc)==1 % on regarde si cycle 0 ou 1
                        
                        LaunchDate=T.juld.data(iLaunch);
                        
                        %%%reperer si LaunchDate a une valeur realiste
                        launch_date = LaunchDate+datenum('01011950','ddmmyyyy');
                        is_ok = (launch_date >= datenum('01011997','ddmmyyyy')) & (launch_date < datenum(date));
                        str2double(strtrim(allfloats{1}));
                        
                        if is_ok == 0
                              
                              fprintf(fid_alerte,'%s\n',[ floatname,' PB LAUNCH DATE NON REALISTE ', datestr(launch_date), '. FLAG LAUNCH MIS A 6.']);
                              fprintf('%s\n',[ floatname, 'PB LAUNCH DATE NON REALISTE ', datestr(launch_date), '. FLAG LAUNCH MIS A 6.']);
                              M.launch_qc = 6;
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
                                      M.launch_qc = 6;                    %%% Flag à 6 pour la launch date.
                                      fprintf(fid_alerte,'%s\n',[ floatname,' PB LAUNCH DATE ANTERIEURE A START DATE DE ', num2str(Diff_LS), ' j. FLAG LAUNCH MIS A 6.']);
                                      fprintf('%s\n',[ floatname, 'PB LAUNCH DATE ANTERIEURE A START DATE DE ', num2str(Diff_LS), ' j. FLAG LAUNCH MIS A 6.']);
                                      if Stat == 1
                                          alertCyc_e6 = [alertCyc_e6 1];
                                          alerte20(k)=str2double(floatname);
                                      end
                                      
                                  elseif M.startup_date_qc=='1' && Diff_LS>1
                                      fprintf(fid_alerte,'%s\n',[floatname, ' ECART ENTRE LAUNCH DATE et START DATE > 1jour, PB LAUNCHDATE ?']);
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
                               isokLaunchDate(k)=ecartLaunchDateFirstLoc(k)>-0.125; %  NC faire des stats sur ces ecarts entre launch date et premiere position: definir un seuil
                        
                               
                              if isokLaunchDate(k)  & M.launch_qc==1
                                 %on verifie que les dates de loc du premier cycle sont sup àla Launch date sinon, on les flaggue
                                 % isbad = T.juld.data(idSorted) <= LaunchDate;
                                 isbad = T.juld.data(idCyc(idSorted)) <= LaunchDate;
                                 isbad_post = T.juld.data(idCyc_post(idSorted_post)) <= LaunchDate;
                                 %T.juld.data(idCyc(idSorted(isbad)))=6;  
                                 % T.juld_adjusted_qc(idCyc(idSorted(isbad)))=6;                           
                                 
                                 if sum(isbad)>0  & sum(isbad_post)== 0  % SI L'ERREUR NE VIENT PAS D'UNE LAUNCH DATE ERRONEE
                                     fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DATES DE LOC DU 1er CYCLE ANTERIEURES A LAUNCH DATE > 3h']);
                                     fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DATES DE LOC DU 1er CYCLE ANTERIEURES A LAUNCH DATE > 3h']);
                                     T.juld.data(idSorted(isbad))=2;                 %%%%FLAG � 2 les dates Loc du 1er cycle (douteuses)
                                     if Stat == 1
                                        alertCyc_e5 = [alertCyc_e5 cycles_sorted(id)];
                                        alerte21(k,cycles_sorted(id)+1)=str2double(floatname);
                                       
                                     end
                                 end
                              elseif isokLaunchDate(k) < -0.125       %% SI L'ERREUR VIENT D'UNE LAUNCH DATE ERRONEE
                                 
                                 fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ? PB LAUNCH DATE POSTERIEURE à DATE DE LOC DE PLUS DE 3 h?']);
                                 fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ? PB LAUNCH DATE POSTERIEURE à DATE DE LOC DE PLUS DE 3 h?']);
                                 if Stat==1
                                     alertCyc_e5 = [alertCyc_e5 cycles_sorted(id)];
                                     alerte22(k,id)=str2double(floatname);
                                 end
                             end
                           
                        
                           if ecartLaunchDateFirstLoc(k)>10.5 % Cas � part ou ecart entre Launch date et dates de LOC tres grand
                              
                              fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ? PB ECART LAUNCH DATE ET DATE DE LOC DE PLUS DE 10.5 j. DPF sans PRELUDE ?']);
                              fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', ? PB ECART LAUNCH DATE ET DATE DE LOC DE PLUS DE 10.5 j. DPF sans PRELUDE ?']);
                              
                              if Stat==1
                                  %alertCyc_e4 = [alertCyc_e4 cycles_sorted(id)];
                                  alerte23(k,cycles_sorted(id)+1)=str2double(floatname);
                                  
                              end
                                  
                           end
                        end
                            
                      
                        %%%VERIFICATION PHASE DE PRELUDE (si pas agglomere avec cycle 1
                         %%%dans le cas des Deep Profile First (DPF).
        
                         %%%vérification de la présence du cycle 0 (prélude)
                                             
                        if(isempty(find(cycles==0)))    
                           if(T.measurement_code.data(2)==100) %%%DPF (Deep Profile first) 
                        %if ecartLaunchDateFirstLoc(k)*24>24
                              fprintf(fid_alerte,'%s\n',[ floatname ', cycle '...
                              num2str(cycles(id)) ', PAS DE CYCLE 0, DPF sans prelude']);
                              fprintf('%s\n',[ floatname ', cycle '...
                              num2str(cycles(id)) ', PAS DE CYCLE 0, DPF sans prelude']); % ', LOC plus d''un jour après => DPF sans prelude 
                              if Stat==1
                                alertCyc_e4=[alertCyc_e4 cycles_sorted(id)];
                                alerte24(k,cycles_sorted(id)+1)=str2double(floatname);
                              
                              end
                           end
                        end
                        
                   elseif numCycleMinTraj>1 & id==1
                        fprintf(fid_alerte, '%s\n',[ floatname ', Commence au cycle n� ', num2str(cycles(1)), ' ecartLaunchDateFirstLoc mis a� NaN']);
                        fprintf('%s\n',[ floatname ', Commence au cycle n� ', num2str(cycles(1)), ' ecartLaunchDateFirstLoc mis a� NaN']);
                        ecartLaunchDateFirstLoc(k)= NaN;  
                    end
                   
                else
                    
                     ecartLaunchDateFirstLoc(k)=LocDate_sorted(1)-LaunchDate;
                     
                     if(ecartLaunchDateFirstLoc(k)>100)
                          
                          fprintf(fid_alerte,'%s\n',[ floatname ', cycle '...
                          num2str(cycles(id)) ', UNE SEULE LOC et DATE INCOHERENTE.']);   
                          fprintf('%s\n',[ floatname ', cycle '...
                          num2str(cycles(id)) ', UNE SEULE LOC et DATE INCOHERENTE.']);
                          if Stat == 1
                              alertCyc_e4 = [alertCyc_e4 cycles_sorted(id)];
                              
                          end
                     else
                         fprintf(fid_alerte,'%s\n',[ floatname ', cycle '...
                         num2str(cycles(id)) ', UNE SEULE LOC.']);   
                         fprintf('%s\n',[ floatname ', cycle '...
                         num2str(cycles(id)) ', UNE SEULE LOC.']); 
                         if Stat == 1
                              alertCyc_e4 = [alertCyc_e4 cycles_sorted(id)];
                              alerte25(k,cycles_sorted(id)+1)=str2double(floatname);
                             
                         end
                     end
                
                        ecartLaunchDateFirstLoc(k) = NaN;
                        
                  end   %  fin de la condition sur locDate_sorted
               
                
                % CONTROLE DES POSITIONS                       
                %-----------------------------------
                % si est bien dans l'eau et non sur la terre
                % on flagge 4 si ce n'est pas le cas
                dlat=locLat_sorted(1);
                dlon=locLon_sorted(1);
               
                %%pour recuperer   (temporaire)
                loclon_first(id) = dlon;
                loclat_first(id) = dlat;
                
                
                dlatend = locLat_sorted(end);
                dlonend = locLon_sorted(end);
                dlocPosition_qc_sorted = locPosition_qc_sorted(1);
                dlocPosition_qc_sortedend = locPosition_qc_sorted(end);
                % on ne verifie que les mesures qui ont passe les tests
                % temporels
               % if (locPosition_qc_sorted~=6|locPosition_qc_sorted~=4)&(length(locLon_sorted)>1&length(locLat_sorted)>1); %on ne regarde que les flags corrects ou non traités
              
                                           
                if ~isnan(dlat)&~isnan(dlon)
                    
                   if (dlocPosition_qc_sorted<4)&(dlocPosition_qc_sortedend<4)% &(length(locLon_sorted)>1&length(locLat_sorted)>1); %on ne regarde que les flags corrects ou non traités (de la 1ère et dernière loc)
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
                                        fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB DE LOCALISATION SUR LE CONTINENT. FLAG POSITION MIS A 4.']);
                                        fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB DE LOCALISATION SUR LE CONTINENT.FLAG POSITION MIS A 4.']);
                                        if Stat == 1
                                             alertCyc_e8 = [alertCyc_e8 cycles_sorted(id)];
                                             alerte26(k,cycles_sorted(id)+1)=str2double(floatname);
                                        end
                                        if(map==1)
                                           m_plot(LONG(ilong,1),LAT(ilat,1),'color','g','marker','o','markersize',5);
                                        end
                                        T.position_qc.data(idCyc)=4;
                                    else
                                        fprintf(fid_alerte,'%s\n',[floatname ', cycle ' num2str(cycles_sorted(id)) ' , LOC A TERRE: altitude du flotteur > 0 (' (num2str(elev_float)) 'm) , mais au moins 1 point voisin < 0. Peut être coincé dans un atoll ou récupéré et proche de la terre']); 
                                        fprintf('%s\n',[floatname ', cycle ' num2str(cycles_sorted(id)) ' , LOC A TERRE: altitude du flotteur > 0 (' (num2str(elev_float)) 'm) , mais au moins 1 point voisin < 0. Peut être coincé dans un atoll ou récupéré et proche de la terre']); 
                                        if Stat == 1
                                            alertCyc_e8 = [alertCyc_e8 cycles_sorted(id)];
                                            alerte27(k,cycles_sorted(id)+1)=str2double(floatname);    
                                          
                                        end
                                        
                                    end
                                
                            elseif elev_float<0&elev_float>-200&isempty(find(voisins>0))==0    %%%si l'atitude du point est >-200 et l'altitude d'au moins un voisin >0 est probablement proche de la terre.
                                   fprintf(fid_alerte,'%s\n',[floatname ', cycle ' num2str(cycles_sorted(id)) ' ,  altitude du flotteur proche de 0 (' (num2str(elev_float)) 'm) , et au moins 1 point voisin > 0. Probablement proche terre.']); 
                                   fprintf('%s\n',[floatname ', cycle ' num2str(cycles_sorted(id)) ' ,  -20 < altitude du flotteur < 0 (' (num2str(elev_float)) 'm) , et au moins 1 point voisin > 0. Probablement proche terre.']); 
                                   if Stat == 1
                                      alertCyc_e8 = [alertCyc_e8 cycles_sorted(id)];
                                   end
                            end
                            
                                
                            %       if(elev_float>-150 & length(alerte28(k,cycles_sorted(id-15):end))>5 & mean(voisins)>-150)  %%%si elev proche surface et plusieurs loc terre aux cycles précédents
                             %         fprintf(fid_alerte,'%s\n',[floatname ', cycle ' num2str(cycles_sorted(id)) ' ,  altitude du flotteur > -150 (' (num2str(elev_float)) 'm) , et au moins 5 cycles précédents en loc terre. Probablement proche terre.']); 
                             %         fprintf('%s\n',[floatname ', cycle ' num2str(cycles_sorted(id)) ' ,  altitude du flotteur > -150 (' (num2str(elev_float)) 'm) , et au moins 5/15 cycles précédents en loc terre. Probablement proche terre.']); 
                             %         if Stat == 1
                             %             alertCyc_e8 = [alertCyc_e8 cycles_sorted(id)];
                             %         end
                             %      end                                               
                              
                            %end
                              
                      end
                       
                   end    %%fin de la condition sur le qc
                  
                end       %%% fin de la condition sur dlat et dlon    
             
                
                
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
                    
                   
                    if  id>1 &  (dlocPosition_qc_sorted<4)&(dlocPosition_qc_sortedend<4) & ~isnan(dureeCycle) & sum(idCycprec{id})>0 
                        dlatlas=T.latitude.data(idCycprec{id}(end));
                        dlonlas=T.longitude.data(idCycprec{id}(end));
                        distanceprof=a*acos(sin(dlat*pi/180)*sin(dlatlas*pi/180)+cos(dlat*pi/180)*cos(dlatlas*pi/180)*cos((dlonlas-dlon)*pi/180));
                        distancederiveprof(istat)=distanceprof; %regroupe les valeurs de dérive en profondeur
                        % la borne varie en fonction de la duree de cycle
                        % calculee et il ne faut pas que la vitesse du
                        % profileur depasse 3m/s
                        
                          
                        %%%%EN PROFONDEUR
                        if distanceprof> dureeCycle*(cycles_sorted(id)-cycles_sorted(id-1))*(24*60*60)*3 % borne a adapter en fonction de la duree du cycle (vitesse maximale du flotteur : 3 m/s)
                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',DERIVE EN PROFONDEUR > ' num2str(dureeCycle*(cycles(id)-cycles(id-1))*(24*60*60)*3/1000) 'km :' num2str(distanceprof/1000) ' km. Loc Argos probablement erronée.']);
                            fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',DERIVE EN PROFONDEUR > ' num2str(dureeCycle*(cycles(id)-cycles(id-1))*(24*60*60)*3/1000) 'km :' num2str(distanceprof/1000) ' km. Loc Argos probablement erronée.']);
                            if Stat == 1
                                alertCyc_e9 = [alertCyc_e9 cycles_sorted(id)];
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
                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', DERIVE EN SURFACE ENTRE 1er et DERNIERE LOC D'' UN CYCLE >' num2str(bornesurf/1000) 'km :' num2str(distancesurf/1000) ' km. Loc Argos probablement erronée.']);
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
                    if(id>1 & sum(idCycprec{id})>0)     
                       fid_koba = fopen([DIR_HOME,'log_koba.txt'],'wt');  
                       locDate_prec = T.juld.data(idCycprec{id}); isok = find(locDate_prec(T.position_qc.data(idCycprec{id})<4));
              
                       if(isempty(isok)==0)  
                  
                           %%% exclu le cas ou tous les qc seraient a� 4 (qd flotteur sur le continent par ex)
                           locDate_prec_ok = locDate_prec(isok);
                           locLat_sorted(find(locPosition_qc_sorted==6))=NaN; locLon_sorted(find(locPosition_qc_sorted==6))=NaN;                                                           
                           %%test de Koba sur les vitesses 
                  
                           [o_date, o_longitude, o_latitude, o_posAcc, o_posQcIn, o_posQcOut, ...
                           o_speed,o_idBadPos, Little_Time, Little_Time_Sp, deltaTime_Sp, normal_Time_Sp] = check_argos_positions_koba_all_gh(LocDate_sorted,...
                           locLon_sorted, locLat_sorted, locQc, locPosition_qc, ...
                           locDate_prec_ok(end), dlonlas, dlatlas, fid_koba);
                           %if(o_posQcOut(find(o_speed==max(o_speed)))==1)
                           %if(isempty(find(o_speed>76.88 & o_speed<76.90))==0)
                           %keyboard
                           %end
 
                           %%Pour l'affichage des vitesse brutes et flaguees
                   
                           if(isempty(o_idBadPos)==0 & o_idBadPos~=length(o_date))
                              o_idBad_speed = [o_idBadPos, o_idBadPos+1];
                           else
                              o_idBad_speed = o_idBadPos;
                           end
                           o_speed_isok = o_speed;
                           if(isempty(o_idBad_speed)==0)
                              o_speed_isok(o_idBad_speed)=NaN;
                           end
                           speed_Cyc =  [speed_Cyc, o_speed_isok];
                           if(isempty(idDoublon_date)==0)
                              ppp = ppp+length(idDoublon_date)
                              Speed_doublon(ppp-(length(idDoublon_date)-1):ppp,:) = o_speed_isok(idDoublon_date);
                           end
                           o_speed_isok(idDoublon_date)=NaN;  %%car koba ne prend pas en compte les doublons de date si pas associés à doublons de lon/lat
                           o_speed_isokd = o_speed_isok;
                  

                           if(isempty(o_idBadPos)==0)
                              fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', LOCA ARGOS ERRONEE (d''apres test de Koba) Position ' num2str(o_idBadPos') '. FLAG POSITION MIS A 6. ']);
                              fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', LOCA ARGOS ERRONEE (d''apres test de Koba) Position '  num2str(o_idBadPos')  '. FLAG POSITION MIS A 6. ']);
                              if Stat == 1
                                  alerte30(k,cycles_sorted(id)+1)=str2double(floatname);
                                  
                              end
                              locPosition_qc(o_idBadPos)=6;   %%%flag à 6 les positions reperees comme erronees.
                              
                              if(Stat == 1 & isempty(find(o_idBadPos==1))==0)
                                  alertCyc_e5 = [alertCyc_e5 cycles_sorted(id)];
                              end
                              %T.position_qc.data(idCyc(o_idBadPos))=6;   
                           end
                        
                        end
                    end   %%fin de la condition sur id>1
               end  %%fin de la condition (sum(idCyc) > 0 & id==1) || (sum(idCyc)>0  & sum(idCycprec{id}) > 0)  
            % idCycprec=idCyc_sorted;   
             
            end   %%%Fin de la boucle sur les cycles
            
            %%%Ecriture des alertes dans un fichier txt pour STATS
            
            if(Wr_Err==1)
                if(isempty(alertCyc_e3)==0)
                   fprintf(fid3,'%s\n',[floatname ', [' num2str(alertCyc_e3) ']']);   %%% Grounded        
                end
                if(isempty(alertCyc_e4)==0)
                   fprintf(fid4,'%s\n',[floatname ', [' num2str(alertCyc_e4) ']']);   %%% alerte cycle 
                end
                if(isempty(alertCyc_e5)==0)
                   fprintf(fid5,'%s\n',[floatname ', [' num2str(alertCyc_e5) ']']);    %%% alerte locdate
                end
                if(isempty(alertCyc_e6)==0)
                   fprintf(fid6,'%s\n',[floatname ', [' num2str(alertCyc_e6) ']']);   %%% alerte launchdate                
                end
                if(isempty(alertCyc_e7)==0)
                   fprintf(fid7,'%s\n',[floatname ', [' num2str(alertCyc_e7) ']']);   %%% alerte pression
                end
                if(isempty(alertCyc_e8)==0)
                   fprintf(fid8,'%s\n',[floatname ', [' num2str(alertCyc_e8) ']']);   %%% alerte loc terre                          
                end
                if(isempty(alertCyc_e9)==0)
                   fprintf(fid9,'%s\n',[floatname ', [' num2str(alertCyc_e9) ']']);   %%% alerte loc pos                          
                end
                if(isempty(alertCyc_e10)==0)
                   fprintf(fid10,'%s\n',[floatname ', [' num2str(alertCyc_e10) ']']);   %%% Incoherence meta/traj pour durée de cycles                         
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

    end  % si fichier traj trouve
    
    
    
end  %% fin de la boucle sur les flotteurs

disp('*****FIN DES TESTS******')

fclose(fid3)
fclose(fid4);fclose(fid5); fclose(fid6);fclose(fid7); fclose(fid8);fclose(fid9);
fclose(fid10);

if Stat==1
     Diff_Cyc = Diff_Medcyc_Mcyc(Diff_Medcyc_Mcyc>0);
     Name_fileStat = [DIR_STAT 'variables_stat.mat'];
     save(Name_fileStat,'Diff_Medcyc_Mcyc','Diff_Ppark','Diff_Cyc',...
    'Diff_LaunchD','Diff_LaunchLat','Diff_LaunchLon','ecartMaxLocCycle_fl',...
    'ecartMeanLocCycle_fl','ecartFirstLastLoc','ecartLaunchDateFirstLoc',...
    'ecartLaunchDateFirstLoc','distancederiveprof',...
     'distancederivesurf','pres_drift_mes','temp_drift_mes','psal_drift_mes');

    Name_fileStat2 = [DIR_STAT 'alertes_stat.mat'];
     save(Name_fileStat2,'alerte1', 'alerte2', 'alerte3', 'alerte4', 'alerte5', 'alerte6', ...
     'alerte7', 'alerte8', 'alerte9', 'alerte10', 'alerte11', 'alerte12', 'alerte13', 'alerte14', ...
     'alerte15', 'alerte16', 'alerte17', 'alerte18','alerte19', 'alerte20', 'alerte21', ...
     'alerte21', 'alerte22', 'alerte23', 'alerte24', 'alerte25', 'alerte26', 'alerte27', ...
     'alerte28', 'alerte29', 'alerte30');
end





