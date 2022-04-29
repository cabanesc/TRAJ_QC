function param = config()
addpath('/home/lops/users/ccabanes/matlab/libs_cc/Myfun/')
addpath('./lib/lib_ext')
param.DIR_TOOL = '/home/pharos/andro/gaelle/soft/tool_qc_traj/';
%param.DIR_SOFT = '/home/pharos/andro/gaelle/soft/decArgo_20180607_022a/';
param.DIR_SOFT = '/home/pharos/andro/livraisons/Livraison_0712021/soft/DecArgo_soft/';
param.DIR_VISU = '/home/pharos/andro/gaelle/soft/visu_tool_traj/';
param.DIR_LISTE = '/export/home1/ccabanes/data/TRAJ_QC/listes/';
%param.DIR_LISTE ='/home/pharos/andro/soft/listes/';  %  cc pour tests uniquement
param.DIR_STAT = '/export/home1/ccabanes/data/TRAJ_QC/stat/';

%DIR_FTP='/home5/pharos/REFERENCE_DATA/ARGO/DATA/201901-ArgoData/';
%param.DIR_FTP ='/home5/pharos/REFERENCE_DATA/ARGO/DATA/201901-ArgoData/';
%param.DIR_FTP ='/home5/pharos/REFERENCE_DATA/ARGO/DATA/202001-ArgoData/';
%param.DIR_FTP ='/export/home1/ccabanes/data/DOI_ARGO/201906-ArgoData/';
param.DIR_FTP ='/home5/pharos/REFERENCE_DATA/ARGO/DATA/202009-ArgoData/dac/';
param.DAC ={'coriolis','aoml','bodc','incois','csio','csiro','jma','nmdis','meds','kma','kordi'};
%DIR_BATHY='/home5/pharos/REFERENCE_DATA/BATHYMETRY/DATA/';%
%DIR_BATHY='/home/gherbert/andro/andro/bathy/ETOPO2v2c_f4.nc';
%DIR_ISAS='/home/lpo4/ISAS_LPO/CLIM_ISAS13/';
param.DIR_ISAS = '/home5/pharos/REFERENCE_DATA/OCEAN_REP/ISAS/ISAS13/ANA/CLIM/'
param.DIR_HOME = '/export/home1/ccabanes/data/TRAJ_QC/';

param.Wr_Err = 1;   %%0: n'écrit pas les erreurs dans les fichiers txt  1: l'écrit.

param.Bathy=2  ; 

if(param.Bathy==1)
    param.DIR_BATHY='/home5/pharos/REFERENCE_DATA/BATHYMETRY/DATA/ETOPO2v2c_f4.nc';
    % Lecture du modele bathymetrique ETOPO2_v2
    bathy_file = [param.DIR_BATHY] % 'ETOPO2v2c_f4.nc'];
    C=read_netcdf_allthefile(bathy_file);
    LONG=C.x.data;
    LAT=C.y.data;
    ELEV=C.z.data;
elseif param.Bathy==2  % SRTM30
    %param.DIR_BATHY ='/home/pharos/andro/data/bathy/SRTM30plus/data/';
	 param.DIR_BATHY ='/export/home1/ccabanes/data/PHAROS_BCK/andro/data/bathy/SRTM30plus/data/';
elseif param.Bathy==3  % GEBCO
     param.DIR_BATHY ='/home/pharos/andro/livraisons/Livraison_0712021/soft/DecArgo_soft/_ressources/GEBCO_2020/';
end

param.map=0;      %%%1 affiche une carte et localise les erreurs.

param.Stat = 1;   %%% 1: r�cup�re les cycles concern�es par des erreurs donn�es   ; 0: ne r�cup�re pas

param.DIR_NEW_TRAJ_FILE ='/export/home1/ccabanes/data/TRAJ_QC/trajFiles/';    % add cc 26/10/2020
param.save_traj_file = 1;  %%% 1: sauve les fichiers traj   0: ne sauve pas   % add cc 26/10/2020 

param.ATLAS_FILE ='/export/home1/ccabanes/data/TRAJ_QC/ATLAS/uvPTS_RTQC_atlas.dat';    % add cc 15/01/2021


% liste des flotteurs à traiter : plusisuers listes possibles % add cc 02/11/2020
% param.liste_rep_alerte  = {'AOML_RISER_apex_argos'; 'JMA_APEX_PI_JAMSTEC'; 'CSIRO_APEX_argos_apf9';'CSIRO_APEX_argos_apf8'; 'INCOIS_APEX_Argos_apf8'; 'INCOIS_APEX_Argos_apf9';'AOML_SOLOII_Dmode_all'; 'AOML_SOLOII_Dmode_all_bis'};
%CORIOLIS_2020_PROVOR_ARVOR_IRIDIUM.txt
%param.liste_rep_alerte = {'AOML_SOLOII_Dmode_all'};


%param.liste_rep_alerte  = {'CO_ARVOR_MOCCA_2019';'AO_APEX_ARGOS_2019'; 'AO_SOLO_IRIDIUM_2021';  'CO_PROVOR_IRIDIUM_R_2018'; 'AO_APEX_ARGOS_2021'; 'CO_APEX_ARGOS_2017'; 'CO_NEMO_IRIDIUM_2021'; 'CO_PROVOR_IRIDIUM_SBD_2018'; 'AO_APEX_ARGOS_APF8_2018'; 'CO_APEX_ARGOS_2020'; 'CO_PROVOR_ARGOS_2016'; 'CS_APEX_ARGOS_APF8_2019'; 'AO_APEX_ARGOS_APF9_2018'; 'CO_APEX_ARGOS_2021'; 'CO_PROVOR_ARGOS_2017'; 'CS_APEX_ARGOS_APF9_2019'; 'AO_APEX_IRIDIUM_2021'; 'CO_APEX_IRIDIUM_2020'; 'CO_PROVOR_ARGOS_2020'; 'IN_APEX_ARGOS_APF8_2019'; 'AO_SOLO_ARGOS_2021'; 'CO_APEX_IRIDIUM_2021'; 'CO_PROVOR_ARGOS_2021'; 'IN_APEX_ARGOS_APF9_2019'; 'AO_SOLO_ARGOS_ROEMMICH_2021'; 'CO_APEX_IRIDIUM_NAVIS_2018'; 'CO_PROVOR_IRIDIUM_2016'; 'JM_APEX_ARGOS_2019'; 'AO_SOLOII_IRIDIUM_2019'; 'CO_ARVOR_IRIDIUM_2020'; 'CO_PROVOR_IRIDIUM_2021'};

param.liste_rep_alerte  = {'CO_ARVOR_MOCCA_2019';'AO_APEX_ARGOS_2019'; 'AO_SOLO_IRIDIUM_2021';  'CO_PROVOR_IRIDIUM_R_2018'; 'AO_APEX_ARGOS_2021'; 'CO_APEX_ARGOS_2017'; 'CO_NEMO_IRIDIUM_2021'; 'CO_PROVOR_IRIDIUM_SBD_2018'; 'AO_APEX_ARGOS_APF8_2018'; 'CO_APEX_ARGOS_2020'; 'CO_PROVOR_ARGOS_2016'; 'AO_APEX_ARGOS_APF9_2018'; 'CO_APEX_ARGOS_2021'; 'CO_PROVOR_ARGOS_2017'; 'AO_APEX_IRIDIUM_2021'; 'CO_APEX_IRIDIUM_2020'; 'CO_PROVOR_ARGOS_2020'; 'AO_SOLO_ARGOS_2021'; 'CO_APEX_IRIDIUM_2021'; 'CO_PROVOR_ARGOS_2021'; 'AO_SOLO_ARGOS_ROEMMICH_2021'; 'CO_APEX_IRIDIUM_NAVIS_2018'; 'CO_PROVOR_IRIDIUM_2016';'AO_SOLOII_IRIDIUM_2019'; 'CO_ARVOR_IRIDIUM_2020'; 'CO_PROVOR_IRIDIUM_2021'};

%param.liste_rep_alerte = {'TEST'};
% Config pour la comparaison des alertes  % add cc 02/11/2020
param.liste_alertes = {'grounded';'pressure';'locpos';'locdate';'cycle';'locterre';'launchdate';'metadur'};


for k=1:length(param.liste_rep_alerte)
%param.Liste_Float{k} = [param.DIR_LISTE 'Liste_livraison2019/' param.liste_rep_alerte{k} '/' param.liste_rep_alerte{k} '.txt'];
param.Liste_Float{k} = [param.DIR_LISTE 'Liste_all_livraison/' param.liste_rep_alerte{k} '.txt'];
end

end
