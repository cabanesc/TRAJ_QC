function param = config()

param.DIR_TOOL = '/home/gherbert/andro/andro/soft/tool_qc_traj/';
param.DIR_SOFT = '/home/gherbert/andro/andro/soft/decArgo_20180607_022a/';
param.DIR_VISU = '/home/gherbert/andro/andro/soft/visu_tool_traj/';
param.DIR_LISTE = '/home/gherbert/andro/andro/soft/listes/';
param.DIR_STAT = '/home/gherbert/andro/andro/stat/';

%DIR_FTP='/home5/pharos/REFERENCE_DATA/ARGO/DATA/201901-ArgoData/';
param.DIR_FTP ='/net/alpha/exports/sciences/data/REFERENCE_DATA/ARGO/DATA/201901-ArgoData/';
param.DAC ={'coriolis','aoml','bodc','incois','csio','csiro','jma','nmdis','meds','kma','kordi'};
%DIR_BATHY='/home5/pharos/REFERENCE_DATA/BATHYMETRY/DATA/';
%DIR_BATHY='/home/gherbert/andro/andro/bathy/ETOPO2v2c_f4.nc';
%DIR_ISAS='/home/lpo4/ISAS_LPO/CLIM_ISAS13/';
param.DIR_ISAS = '/net/alpha/exports/sciences/data/LPO_ISAS/CLIM_ISAS13/'

param.DIR_HOME = '/home/gherbert/andro/andro/soft/solene/';

param.Wr_Err = 1;   %%0: n'Ã©crit pas les erreurs dans les fichiers txt  1: l'Ã©crit.

param.Bathy=0  ; 

if(param.Bathy==1)
    param.DIR_BATHY='/home/gherbert/andro/andro/bathy/ETOPO2v2c_f4.nc';
    % Lecture du modele bathymetrique ETOPO2_v2
    bathy_file = [param.DIR_BATHY] % 'ETOPO2v2c_f4.nc'];
    C=read_netcdf_allthefile(bathy_file);
    LONG=C.x.data;
    LAT=C.y.data;
    ELEV=C.z.data;
else
    param.DIR_BATHY ='/home/gherbert/andro/andro/data/bathy/SRTM30plus/data/';
end

param.map=0;      %%%1 affiche une carte et localise les erreurs.

param.Stat = 1;   %%% 1: récupère les cycles concernées par des erreurs données   ; 0: ne récupère pas

param.DIR_NEW_TRAJ_FILE ='/export/home1/ccabanes/data/TRAJ_QC/trajFiles/';    % add cc 26/10/2020
param.save_traj_file = 1;  %%% 1: sauve les fichiers traj   0: ne sauve pas   % add cc 26/10/2020 


% liste des flotteurs Ã  traiter : 
%  cc 02/11/2020 pour compa_altran_alertes: plusieurs listes possibles
%                pour Prog_alertes_traj: pour le moment 1 seule liste (pas encore testÃ© avec plusieurs)
% param.liste_rep_alerte  = {'AOML_RISER_apex_argos'; 'JMA_APEX_PI_JAMSTEC'; 'CSIRO_APEX_argos_apf9';'CSIRO_APEX_argos_apf8'; 'INCOIS_APEX_Argos_apf8'; 'INCOIS_APEX_Argos_apf9';'AOML_SOLOII_Dmode_all'; 'AOML_SOLOII_Dmode_all_bis'};
param.liste_rep_alerte = {'INCOIS_APEX_Argos_apf8'};
% Config pour la comparaison des alertes  % add cc 02/11/2020
param.liste_alertes = {'grounded';'pressure';'locpos';'locdate';'cycle';'locterre';'launchdate';'metadur'};


for k=1:length(param.liste_rep_alerte)
param.Liste_Float{k} = [param.DIR_LISTE 'Liste_livraison2019/' param.liste_rep_alerte{k} '/' param.liste_rep_alerte{k} '.txt'];
end

%param.Liste_Float = [param.DIR_LISTE 'liste_flotteurs_argos_apex_csiro.txt';  
%param.Liste_Float = [param.DIR_LISTE 'liste_flotteurs_argos_aoml_other.txt']  ;  
%param.Liste_Float = [param.DIR_LISTE 'liste_flotteurs_iridium_coriolis.txt'];
%param.Liste_Float = [param.DIR_LISTE 'liste_flotteurs_argos_coriolis.txt']     
%param.Liste_Float =  [param.DIR_LISTE 'Liste_livraison2019/2_172_AOML_RISER_apex_argos.txt'];
%param.Liste_Float = [param.DIR_LISTE 'Liste_livraison2019/4_114_JMA_APEX_PI_JAMSTEC.txt'];
%param.Liste_Float = [param.DIR_LISTE 'Liste_livraison2019/CSIRO_APEX_argos_apf8.txt'];
%param.Liste_Float =  [param.DIR_LISTE 'Liste_livraison2019/CSIRO_APEX_argos_apf9.txt'];
%param.Liste_Float = [param.DIR_LISTE 'Liste_livraison2019/CSIRO_APEX_argos_apf8/CSIRO_APEX_argos_apf8.txt'];
%param.Liste_Float = [param.DIR_LISTE 'Liste_livraison2019/INCOIS_APEX_Argos_apf8.txt'];
%param.Liste_Float = [param.DIR_LISTE 'Liste_livraison2019/INCOIS_APEX_Argos_apf9.txt'];
%param.Liste_Float = [param.DIR_LISTE 'Liste_livraison2019/AOML_SOLOII_Dmode_all_bis.txt'];
%param.Liste_Float = [param.DIR_LISTE 'Liste_livraison2019/AOML_SOLOII_Dmode_all.txt'];


end
