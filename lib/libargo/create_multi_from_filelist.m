
% C.Lagadec / Septembre 2015
% programme adapte d'un programmme de Cecile
% lecture de tous les monoprofils d'un flotteur
% reconstitution d'un multiprofils avec simplement 1er profil
% recuperation des donnes dans structure FLm

% attention : si DOXY 
%  fichiers MR (= BR + R) et fichiers MD ( = D + BD)
% on r�cup�re tous les fichiers ayant M dans le nom
% sans distinction de TEmps differe (MD) ou temps reel (MR)

% dans verif_flag : l'utilisateur donne les bornes 
%                   pour traitement temps reel
% dans plotdata : on traite tous les types de fichiers

function [Fmul,Dimul,file_list]=create_multi_from_filelist(floatname,dacname,DIRFTP,file_list,vertical_sampling_scheme,Param);

%if isempty(dacname)==0
repnc = [DIRFTP dacname '/' floatname '/profiles/'];
%else
%repnc = [DIRFTP dacname '/' floatname '/'];
%end
if isempty(dacname)&isempty(floatname)
repnc = [DIRFTP];
end

Asc_files=file_list;   

FLm=[];
DIMm=[];
%keyboard
is_in_list=ones(length(file_list),1);
for k=1:length(Asc_files) 
    if isempty(Param)
        [FL,DIM,Globatt] = read_netcdf_allthefile([repnc '/' Asc_files{k}]);
    else
        Param.vertical_sampling_scheme.name='VERTICAL_SAMPLING_SCHEME';
        [FL,DIM,Globatt] = read_netcdf_allthefile([repnc '/' Asc_files{k}],Param);
    end
    if isfield(FL,'vertical_sampling_scheme')
        is_primary=findstr_tab(FL.vertical_sampling_scheme.data,vertical_sampling_scheme);
        
        if  sum(is_primary)==1
            [FL,DIM] = extract_profile_dim(FL,DIM,'N_PROF',is_primary);
            [FLm,DIMm]=cat_profile_dim(FLm,FL,DIMm,DIM,'N_PROF') ;
        else
            is_in_list(k)=0;
        end
        
    else
        if DIM.n_prof.dimlength==1
            [FLm,DIMm]=cat_profile_dim(FLm,FL,DIMm,DIM,'N_PROF') ;
        end
    end
      
end
iilist=find(is_in_list==1);
file_list=file_list(iilist);
if isfield(FLm,'pres')
test = check_isfillval_prof(FLm,'pres');
iprof=find(test.pres==1);
iprof_garde=find(test.pres==0);
[FLm,DIMm] = remove_profile_dim(FLm,DIMm,'N_PROF',iprof);
file_list=file_list(iprof_garde);
end
Fmul=FLm;
Dimul=DIMm;
