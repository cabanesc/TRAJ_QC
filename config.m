% CONFIGURATION file for Prog_alertes_traj.m
%==========================================================================
function conf = config()

% Current software directory
%---------------------------
conf.DIR_HOME = '/home/lops/users/ccabanes/dvlpRD/Argo/WorkOnBase/Trajectoire/git/TRAJ_QC/';

% Directory to store data and resuts
%---------------------------
conf.DIR_DATA = '/export/home1/ccabanes/data/TRAJ_QC/';

%conf.DIR_TOOL = '/home/pharos/andro/gaelle/soft/tool_qc_traj/';
%conf.DIR_SOFT = '/home/pharos/andro/gaelle/soft/decArgo_20180607_022a/';

% For m-map
%---------------------------
conf.DIR_MMAP = '/home/pharos/andro/livraisons/Livraison_0712021/soft/DecArgo_soft/_ressources/pourMatlab/';

% Directory of ISAS climatology
%---------------------------
conf.DIR_ISAS = '/home5/pharos/REFERENCE_DATA/OCEAN_REP/ISAS/ISAS13/ANA/CLIM/';

% Directory where the Argo DOI is stored
%---------------------------
%conf.DIR_FTP ='/home5/pharos/REFERENCE_DATA/ARGO/DATA/201901-ArgoData/';
%conf.DIR_FTP ='/home5/pharos/REFERENCE_DATA/ARGO/DATA/202001-ArgoData/';
%conf.DIR_FTP ='/export/home1/ccabanes/data/DOI_ARGO/201906-ArgoData/';
%conf.DIR_FTP ='/home5/pharos/REFERENCE_DATA/ARGO/DATA/202009-ArgoData/dac/';
%conf.DIR_FTP ='/home5/pharos/REFERENCE_DATA/ARGO/DATA/202103-ArgoData/dac/';
conf.DIR_FTP ='/export/home1/ccabanes/data/DOI_ARGO/202203-ArgoData/';
conf.DAC ={'coriolis','aoml','bodc','incois','csio','csiro','jma','nmdis','meds','kma','kordi'};

% Directory of Lists of floats to be processed
%---------------------------
conf.DIR_LISTE = [conf.DIR_DATA 'listes/Liste_update2022/'];
%conf.DIR_LISTE =[conf.DIR_DATA 'listes/'];  %  cc pour tests uniquement

% conf.liste_rep_alerte =  {'AOML_RISER_apex_argos'; 'JMA_APEX_PI_JAMSTEC'; 'CSIRO_APEX_argos_apf9';'CSIRO_APEX_argos_apf8'; 'INCOIS_APEX_Argos_apf8'; 'INCOIS_APEX_Argos_apf9';'AOML_SOLOII_Dmode_all'; 'AOML_SOLOII_Dmode_all_bis'};
%conf.liste_rep_alerte  =  {'AOML_SOLOII_Dmode_all'};
%conf.liste_rep_alerte  =  {'CO_ARVOR_MOCCA_2019';'AO_APEX_ARGOS_2019'; 'AO_SOLO_IRIDIUM_2021';  'CO_PROVOR_IRIDIUM_R_2018'; 'AO_APEX_ARGOS_2021'; 'CO_APEX_ARGOS_2017'; 'CO_NEMO_IRIDIUM_2021'; 'CO_PROVOR_IRIDIUM_SBD_2018'; 'AO_APEX_ARGOS_APF8_2018'; 'CO_APEX_ARGOS_2020'; 'CO_PROVOR_ARGOS_2016'; 'CS_APEX_ARGOS_APF8_2019'; 'AO_APEX_ARGOS_APF9_2018'; 'CO_APEX_ARGOS_2021'; 'CO_PROVOR_ARGOS_2017'; 'CS_APEX_ARGOS_APF9_2019'; 'AO_APEX_IRIDIUM_2021'; 'CO_APEX_IRIDIUM_2020'; 'CO_PROVOR_ARGOS_2020'; 'IN_APEX_ARGOS_APF8_2019'; 'AO_SOLO_ARGOS_2021'; 'CO_APEX_IRIDIUM_2021'; 'CO_PROVOR_ARGOS_2021'; 'IN_APEX_ARGOS_APF9_2019'; 'AO_SOLO_ARGOS_ROEMMICH_2021'; 'CO_APEX_IRIDIUM_NAVIS_2018'; 'CO_PROVOR_IRIDIUM_2016'; 'JM_APEX_ARGOS_2019'; 'AO_SOLOII_IRIDIUM_2019'; 'CO_ARVOR_IRIDIUM_2020'; 'CO_PROVOR_IRIDIUM_2021'};
%conf.liste_rep_alerte  =  {'CO_ARVOR_MOCCA_2019';'AO_APEX_ARGOS_2019'; 'AO_SOLO_IRIDIUM_2021';  'CO_PROVOR_IRIDIUM_R_2018'; 'AO_APEX_ARGOS_2021'; 'CO_APEX_ARGOS_2017'; 'CO_NEMO_IRIDIUM_2021'; 'CO_PROVOR_IRIDIUM_SBD_2018'; 'AO_APEX_ARGOS_APF8_2018'; 'CO_APEX_ARGOS_2020'; 'CO_PROVOR_ARGOS_2016'; 'AO_APEX_ARGOS_APF9_2018'; 'CO_APEX_ARGOS_2021'; 'CO_PROVOR_ARGOS_2017'; 'AO_APEX_IRIDIUM_2021'; 'CO_APEX_IRIDIUM_2020'; 'CO_PROVOR_ARGOS_2020'; 'AO_SOLO_ARGOS_2021'; 'CO_APEX_IRIDIUM_2021'; 'CO_PROVOR_ARGOS_2021'; 'AO_SOLO_ARGOS_ROEMMICH_2021'; 'CO_APEX_IRIDIUM_NAVIS_2018'; 'CO_PROVOR_IRIDIUM_2016';'AO_SOLOII_IRIDIUM_2019'; 'CO_ARVOR_IRIDIUM_2020'; 'CO_PROVOR_IRIDIUM_2021'};
conf.liste_rep_alerte  =  {'TEST'};

for k=1:length(conf.liste_rep_alerte)
    conf.Liste_Float{k} = [conf.DIR_LISTE  conf.liste_rep_alerte{k} '.txt'];
end

% Directory of Bathymety data
%---------------------------
conf.Bathy=2  ;

if(conf.Bathy==1)  % obsolete
    conf.DIR_BATHY='/home5/pharos/REFERENCE_DATA/BATHYMETRY/DATA/ETOPO2v2c_f4.nc';
    % Lecture du modele bathymetrique ETOPO2_v2
    bathy_file = [conf.DIR_BATHY] % 'ETOPO2v2c_f4.nc'];
    C=read_netcdf_allthefile(bathy_file);
    LONG=C.x.data;
    LAT=C.y.data;
    ELEV=C.z.data;
elseif conf.Bathy==2  % SRTM30
    %conf.DIR_BATHY ='/home/pharos/andro/data/bathy/SRTM30plus/data/';
    conf.DIR_BATHY ='/export/home1/ccabanes/data/PHAROS_BCK/andro/data/bathy/SRTM30plus/data/';
elseif conf.Bathy==3  % GEBCO   %
    conf.DIR_BATHY = ['/home/pharos/andro/livraisons/Livraison_0712021/soft/DecArgo_soft/_ressources/GEBCO_2020/'];
end


conf.map = 0;      %%%1 affiche une carte et localise les erreurs.
conf.Stat = 1;   %%% 1: recupere les cycles concernees par des erreurs donnees   ; 0: ne recupere pas %  cc toujours a 1

%%%% OUTPUTS

% Directory to store new Traj files (with QC flags and additional data)
%---------------------------
conf.save_traj_file = 1;  %%% 1: save new traj files   0: do not save new traj files   
conf.DIR_NEW_TRAJ_FILE =[conf.DIR_DATA '/trajFiles/'];  

% File to store NRTQC velocities atlas 
%---------------------------
%conf.ATLAS_FILE =[conf.DIR_DATA '/ATLAS/uvPTS_RTQC_atlas.dat'];
conf.ATLAS_FILE =[conf.DIR_DATA '/ATLAS/uvPTS_RTQC_atlas_TEST.dat'];    % add cc 15/01/2021


% Directory to store mat files of errors found (to make stats and comparisons)
%---------------------------
conf.DIR_ERR = [conf.DIR_DATA 'forStats/'];
conf.Wr_Err = 1;   %%0: n'écrit pas les erreurs dans les fichiers txt  1: l'écrit.
% Types of alerts stored
conf.liste_alertes = {'grounded';'pressure';'locpos';'locdate';'cycle';'locterre';'launchdate';'metadur'};








