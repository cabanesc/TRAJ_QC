function [o_alerte7, o_alerte8, o_alerte9, o_alertCyc_e5, o_alertCyc_e6,o_alertCyc_e9, o_Diff_LaunchD,...
           o_Diff_LaunchLat, o_Diff_LaunchLon] =Test_locs(idLoc,a_cycles_sorted, iLaunch,...
           LaunchDate,LaunchLat, LaunchLon);

%%-----------------------------------------------------------------------
%% Vérifie les locs (dates et positions) et les dates/lon/lat de lancement.
%% Recupere les flotteurs et les cycles concernes par les alertes.
%%
%%INPUT:
%% o_alerte7, o_alerte8,o_alerte9: récupération des flotteurs concernés par les alertes correspondantes
%% o_alerteCyc_e6, o_alerteCyc_e9: récupération des cycles concernés par les alertes correspondantes
%% o_Diff_LaunchD, o_Diff_LaunchLat, o_Diff_LaunchLon: différence dans la date et la position de la launchDate entre fichier méta et traj (pour stats)

%% OUTPUT:
%% idLoc: indice des locs de surfaces
%% a_cycles_sorted: numéros de cycles triés dans l'ordre croissant
%% iLaunch: find(T.cycle.number.data= -1)
%% LaunchDate: date de lancement du fichier traj.
%% LaunchLat:latitude à la date de lancement
%% LaunchLon: longitude à la date de lancement
%%--------------------------------------------------------------------

global PARAM;
global T;
global M;
global P;
global floatname;
global fid_alerte;

o_alerte7=[]; o_alerte8=[]; o_alerte9=[];
o_alertCyc_e5 = []; o_alertCyc_e6= []; o_alertCyc_e9= 0;

% VERIFICATION DES DATES DE LOC ARGOS: DATES REALISTES: later than 1 st January 1997
%----------------------------------------------------------------------------------
% 17167 <= JULD < UTC date of chec
% DATES ET POSITION de MISE A L'EAU
       
        
thedate = T.juld.data(idLoc)+datenum('01011950','ddmmyyyy');
isok = (thedate >= datenum('01011997','ddmmyyyy')) & (thedate < datenum(date));
% trouve les flags 0
isflag0 = T.juld_qc.data(idLoc)==0;
T.juld_qc.data(idLoc(isok&isflag0))=1; % on remplace par 1 les flag 0 uniquement.
T.juld_qc.data(idLoc(~isok))=6;   % A voir si on met un qc special ('6') ou '4'
if sum(~isok)>0
    fprintf(fid_alerte, [ floatname ', cycle ' num2str(sum(~isok)) ', DATES DE LOC NON REALISTES. FLAG DATE MIS A 6.']);
    fprintf('%s\n',[ floatname ', cycle ' num2str(sum(~isok)) ', DATES DE LOC NON REALISTES.FLAG DATE MIS A 6.']);
    if(P.Stat==1)
        o_alertCyc_e5 = [o_alertCyc_e5 T.cycle_number.data(idLoc(isok==0))];
        o_alerte7=str2double(floatname);
     end
            
end
        
       
% VERIFICATION DES LOC ARGOS:
%%%-------------------------------------------------------------------------        
%%%LOCS REALISTES: latitudes et longitudes corrects
        
thelon = T.longitude.data(idLoc);
thelat = T.latitude.data(idLoc);
%isoklon = (thelon>=-180&thelon<=180);    %%%attention  dans certains cas, longitude de 0 à 360. (apex incois ap9) vois si selon le type de flotteur
isoklat = (thelat>=-90&thelat<=90);
%str2double(strtrim(allfloats{1}));
% trouve les flags 0
isflag0 = T.position_qc.data(idLoc)==0;
T.position_qc.data(idLoc(isoklat&isflag0))=1; % on remplace par 1 les flag 0 uniquement.
T.position_qc.data(idLoc(~isoklat))=6;   % 
%T.position_qc.data(idLoc(isoklat&isoklon&isflag0))=1; % on remplace par 1 les flag 0 uniquement.
%T.position_qc.data(idLoc(~isoklon|~isoklat))=6;   % A voir si on met un qc special ('6') ou '4'
if(sum (~isoklat)>0)
%if sum(~isoklat)>0|sum(~isoklon)>0
   fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(sum(~isoklat|~isoklon)) ', LOC NON REALISTES. FLAG POSITION MIS A 6.']);
   fprintf('%s\n',[ floatname ', cycle ' num2str(sum(~isoklat|~isoklon)) ', LOC NON REALISTES. FLAG POSITION MIS A 6.']);
   if P.Stat==1
      o_alertCyc_e9 = [o_alertCyc_e9 T.cycle_number.data(idLoc(isoklat==0))];
      o_alerte8=str2double(floatname);
   end
            
end
        
                    
% 
% VERIFICATION DES DATES ET LOCS DE LANCEMENT (coherence entre les fichiers de trajectoires et de metadonnees)
%---------------------------------
LaunchDate_meta=datenum(M.launch_date,'yyyymmddHHMMSS');
if abs(LaunchLat-M.launch_latitude)>PARAM.LOC_LAUNCH_DIFF_M | abs(LaunchLon-M.launch_longitude)>PARAM.LOC_LAUNCH_DIFF_M |LaunchDate+datenum('01011950','ddmmyyyy')~=LaunchDate_meta
     fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(sum(~isok)) ', INCOHERENCE ENTRE META ET TRAJ POUR LA POSITION ET/OU DATE DE LANCEMENT']);
     fprintf('%s\n',[ floatname ', cycle ' num2str(sum(~isok)) ', INCOHERENCE ENTRE META ET TRAJ POUR LA POSITION ET/OU DATE DE LANCEMENT']);
     T.position_qc.data(iLaunch)=6;% flag special pour confirmer ou non le launch date ?
     if P.Stat==1
         o_alertCyc_e6 = [o_alertCyc_e6 1];
         o_alerte9=str2double(floatname);
      end            
end
       
%%%pour STATS
if P.Stat==1
     o_Diff_LaunchD = LaunchDate+datenum('01011950','ddmmyyyy') - LaunchDate_meta;
     o_Diff_LaunchLat = LaunchLat - M.launch_latitude;
     o_Diff_LaunchLon = LaunchLon - M.launch_longitude;
end
        
end
