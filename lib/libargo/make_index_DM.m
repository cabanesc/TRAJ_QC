function make_index_meta(DIR_FTP,liste_floats_files,cvs_pathfilename)
% -========================================================
%   USAGE : make_index_meta(DIR_FTP,liste_floats_files,cvs_pathfilename)
%   PURPOSE : function qui a partir d'une liste de flotteurs va lire les données meta et cree un fichier csv contenant des informations sur le flotteur
% -----------------------------------
%   INPUT :
%     DIR_FTP   (char)  -path for Argo data
%
%     liste_floats_files (cell of char) - liste des flotteurs  issue du
%     programme select_floats_from_index
%
%     cvs_pathfilename
% -----------------------------------
%   OUTPUT :
%     none
% -----------------------------------
%   HISTORY  : created (2019) ccabanes
%            : modified (yyyy) byxxx 
%   CALLED SUBROUTINES: none
% ========================================================
NcVar.format_version.name='FORMAT_VERSION';
NcVar.platform_type.name='PLATFORM_TYPE';
NcVar.trans_system.name='TRANS_SYSTEM';
NcVar.platform_model.name='PLATFORM_MODEL';
NcVar.inst_reference.name='INST_REFERENCE';
NcVar.wmo_inst_type.name='WMO_INST_TYPE';
NcVar.format_version.name='FORMAT_VERSION';
NcVar.platform_maker.name='PLATFORM_MAKER';
NcVar.firmware_version.name=upper('firmware_version');
NcVar.dac_format_id.name=upper('dac_format_id');
NcVar.pi_name.name=upper('pi_name');
NcVar.project_name.name=upper('project_name');
NcVar.operating_institution.name=upper('operating_institution');
NcVar.standard_format_id.name='STANDARD_FORMAT_ID';
NcVar.controller_board_type_primary.name=upper('controller_board_type_primary');
NcVar.sensor.name=upper('sensor');
NcVar.sensor_model.name=upper('sensor_model');
NcVar.sensor_serial_no.name=upper('sensor_serial_no');
NcVar.launch_date.name=upper('launch_date');
NcVarF.juld.name='JULD';
NcVarF.scientific_calib_date.name='SCIENTIFIC_CALIB_DATE';
NcVarF.parameter.name='PARAMETER';

clear A
fw1=fopen(cvs_pathfilename,'w')
 %fprintf(fw1,'%s\n',['DACNAME;PLATFORM;LAUNCH_DATE;FORMAT_VERSION;TRANS_SYSTEM;PLATFORM_TYPE;PLATFORM_MODEL;CONTROLLER_BOARD;FIRMWARE_VERSION;DAC_FORMAT_ID;STANDARD_FORMAT_ID;WMO_INST_TYPE;SENSOR_MODEL;SENSOR_SERIAL_NO;PI_NAME;PROJECT_NAME;OPERATING_INSTITUTION;DATE_LAST_CT;LAST_CY_R;LAST_CY_D;LAST_CY_BR;LAST_CY_BD; Actif/Not Actif']);
  fprintf(fw1,'%s\n',['DACNAME;PLATFORM;LAUNCH_DATE;TRANS_SYSTEM;PLATFORM_TYPE;SENSOR_MODEL;SENSOR_SERIAL_NO;PI_NAME;PROJECT_NAME;Actif/Not Actif;LAST_CY_R;LAST_CY_D;DATE_LAST_DM_PSAL;LAST_CY_BR;LAST_CY_BD;DATE_LAST_DM_OXY ']);
 
length(liste_floats_files)
 for k=1:length(liste_floats_files)
    disp([num2str(k) '/' num2str(length(liste_floats_files))])
%      thefields=fieldnames(A);
%      for il=1:length(thefields)
%           A.(thefields{il}){k}='';
%      end
%      

    [dacName,rest]=strtok(liste_floats_files{k},'/');
    [floatname,rest]=strtok(rest,'/');
    A.dacname{k}=dacName;
    A.platform_number{k}=floatname;
    FileNameM = [DIR_FTP '/' dacName '/' floatname '/' floatname '_meta.nc'];
	
	% recherche les fichiers R D, BR, BD et trouve le dernier cycle dispo
	[file_list]=select_float_files_on_ftp(floatname,dacName,DIR_FTP,'CD'); 
	if isempty(file_list)==0
    F=read_netcdf_allthefile([DIR_FTP '/' dacName '/' floatname '/profiles/' file_list{end}],NcVarF);
	isparam=findstr_tab(cellstr(squeeze(F.parameter.data(1,end,:,:))),'PSAL');
	if sum(isparam)==1&& isfield(F, 'scientific_calib_date')
	date_temp=squeeze(F.scientific_calib_date.data(1,end,isparam,:))';
	DATE_LAST_DM_PSAL=datestr(datenum(strtrim(date_temp),'yyyymmddHHMMSS'),'yyyy-mm-dd');
	else
	DATE_LAST_DM_PSAL='-';
	end
	else
	DATE_LAST_DM_PSAL='-';
	end
	
	[file_list]=select_float_files_on_ftp(floatname,dacName,DIR_FTP,'BD'); 
	if isempty(file_list)==0
    [F,Dim,G]=read_netcdf_allthefile([DIR_FTP '/' dacName '/' floatname '/profiles/' file_list{end}],NcVarF);
	isparam=findstr_tab(cellstr(squeeze(F.parameter.data(1,end,:,:))),'DOXY')&~(findstr_tab(cellstr(squeeze(F.parameter.data(1,end,:,:))),'_'));
	if sum(isparam)==1&& isfield(F, 'scientific_calib_date')
	date_temp=squeeze(F.scientific_calib_date.data(1,end,isparam,:))';
	DATE_LAST_DM_DOXY=datestr(datenum(strtrim(date_temp),'yyyymmddHHMMSS'),'yyyy-mm-dd');
    else
	DATE_LAST_DM_DOXY='-';
	end
	else
	DATE_LAST_DM_DOXY='-';
	end
	
	[file_list]=select_float_files_on_ftp(floatname,dacName,DIR_FTP,'C'); 
	if isempty(file_list)==0
	[F,Dim,G]=read_netcdf_allthefile([DIR_FTP '/' dacName '/' floatname '/profiles/' file_list{end}],NcVarF);
    LAST_CY_DATE=datestr(F.juld.data(1,:)+datenum('19500101','yyyymmdd'),'yyyy-mm-dd');
	else
	LAST_CY_DATE='-';
	end
	
	if (F.juld.data(1,:)+datenum('19500101','yyyymmdd'))< datenum(date)-20
	Activity='NA';
	else
	Activity='A';
	end
	
	[file_list]=select_float_files_on_ftp(floatname,dacName,DIR_FTP,'CR'); 
	last_CR=find_last_cycle(file_list);
	[file_list]=select_float_files_on_ftp(floatname,dacName,DIR_FTP,'CD');
	last_CD=find_last_cycle(file_list);
	[file_list]=select_float_files_on_ftp(floatname,dacName,DIR_FTP,'BR');
	last_BR=find_last_cycle(file_list);
	[file_list]=select_float_files_on_ftp(floatname,dacName,DIR_FTP,'BD');
    last_BD=find_last_cycle(file_list);

	
    A.platform_type{k}='X';
    A.trans_system{k}='X';
    A.platform_model{k}='X';
    A.format_versionM{k}='X';
    A.platform_maker{k}='X';
    A.firmware_version{k}='X';
    A.standard_format_id{k}='X'; 
    A.dac_format_id{k}='X'; 
    A.wmo_inst_type{k}='X';
    A.pi_name{k}='X';
    A.project_name{k}='X';
    A.operating_institution{k}='X';
    A.controller_board_type_primary{k}='X';
    A.sensor_model{k}='X';
    A.sensor_serial_no{k}='X';
    A.launch_date{k}='X';
    
    if exist(FileNameM)
        M=read_netcdf_allthefile(FileNameM,NcVar);
        if isfield(M,'launch_date')
        A.launch_date{k}=datestr(datenum(strtrim(M.launch_date.data'),'yyyymmddHHMMSS'),'yyyy-mm-dd');
        end
        
        if isfield(M,'platform_type')
        A.platform_type{k}=strtrim(M.platform_type.data');
        end
        if isfield(M,'trans_system')
            if size(M.trans_system.data,1)==16
            A.trans_system{k}=strtrim(M.trans_system.data');
            else
            A.trans_system{k}=strtrim(M.trans_system.data);
            end
        end
        
        if isfield(M,'platform_model')
        A.platform_model{k}=strtrim(M.platform_model.data');
        end
        
        if isfield(M,'format_version')
        A.format_versionM{k}=strtrim(M.format_version.data');
        end
        
        if isfield(M,'platform_maker')
        A.platform_maker{k}=strtrim(M.platform_maker.data');
        end
        
        if isfield(M,'firmware_version')
        A.firmware_version{k}=strtrim(M.firmware_version.data');
        end
        
        if isfield(M,'standard_format_id')
        A.standard_format_id{k}=strtrim(M.standard_format_id.data');
        end
        
        if isfield(M,'dac_format_id')
            A.dac_format_id{k}=strtrim(M.dac_format_id.data');
        else
            if isfield(M,'inst_reference')
                A.dac_format_id{k}=strtrim(M.inst_reference.data');
            end
        end
        
        if isfield(M,'wmo_inst_type')
        A.wmo_inst_type{k}=strtrim(M.wmo_inst_type.data');
        end
        
        if isfield(M,'pi_name')
        A.pi_name{k}=strtrim(M.pi_name.data');
        end
        
        if isfield(M,'project_name')
        A.project_name{k}=strtrim(M.project_name.data');
        end
        
        if isfield(M,'operating_institution')
        A.operating_institution{k}=strtrim(M.operating_institution.data');
        end
        
        if isfield(M,'inst_reference')
        A.inst_reference{k}=strtrim(M.inst_reference.data');
        end
        
        if isfield(M,'controller_board_type_primary')
        A.controller_board_type_primary{k}=strtrim(M.controller_board_type_primary.data');
        end
        
        if isfield(M,'sensor')
            ijk=find(findstr_tab(M.sensor.data,'CTD_CNDC'));
            if isempty(ijk)==0
            if isfield(M,'sensor_model')
            A.sensor_model{k}=strtrim(M.sensor_model.data(ijk,:));
            end
            if isfield(M,'sensor_serial_no')
            A.sensor_serial_no{k}=strtrim(M.sensor_serial_no.data(ijk,:));
            end
            end
        end
    end
	
	
    %fprintf(fw1,'%s\n',[A.dacname{k} ';'  A.platform_number{k} ';' A.launch_date{k} ';' A.format_versionM{k}  ';' A.trans_system{k} ';' A.platform_type{k} ';' A.platform_model{k} ';' A.controller_board_type_primary{k} ';'  A.firmware_version{k} ';' A.dac_format_id{k} ';' A.standard_format_id{k} ';' A.wmo_inst_type{k} ';' A.sensor_model{k} ';' A.sensor_serial_no{k} ';' A.pi_name{k} ';' A.project_name{k} ';' A.operating_institution{k} ';' LAST_CY_DATE ';' last_CR ';' last_CD ';' last_BR ';' last_BD ';' Activity]);
	fprintf(fw1,'%s\n',[A.dacname{k} ';'  A.platform_number{k} ';' A.launch_date{k} ';' A.trans_system{k} ';' A.platform_type{k}  ';' A.sensor_model{k} ';' A.sensor_serial_no{k} ';' A.pi_name{k} ';' A.project_name{k} ';' Activity  ';' last_CR ';' last_CD ';' DATE_LAST_DM_PSAL ';' last_BR ';' last_BD ';' DATE_LAST_DM_DOXY]);
end 
fclose(fw1)

function last_cy=find_last_cycle(file_list)

if isempty(file_list)==0
	[p,r]=strtok(file_list,'_');
	r=strrep(r,'.nc','');
	r=strrep(r,'_','');
	last_cy=r{end};
	
else
last_cy='_';
end

% 
% % faire des listes de flotteurs
% %for il1=1:length(ustatus)
%     %for il2=1:length(utype)
%         %for il3=1:length(utrans)
%             %for il4=1:length(uformat)
%                 %for il5=1:length(umode)
%                     for il6=1:length(udac)
%                     %name_list=['AN_' udac{il6}, '_' ,ustatus{il1},'_',utype{il2},'_',utrans{il3},'_',uformat{il4},'_',umode{il5},'.txt']
%                     %name_list=['GL_', udac{il6}, '_' ,ustatus{il1},'_' ,umode{il5},'.txt']
%                     name_list=['GL_', udac{il6},'.txt']
%                     % trouve les flotteurs de la liste
%                     ijk=find(strcmp(A.dacname,udac{il6})&~(strcmp(A.status,'D')|strcmp(A.status,'X'))&strcmp(A.format_versionTr,'3.1'));
%                     %ijk=find(strcmp(A.dacname,udac{il6})&strcmp(A.status,ustatus{il1})&strcmp(A.modeTr,umode{il5})&~(strcmp(A.status,'D')|strcmp(A.status,'X'))&strcmp(A.format_versionTr,'3.1'));
%                     %ijk=find(strcmp(type_plat,utype{il2})&strcmp(A.dacname,udac{il6})&strcmp(A.trans_system,utrans{il3})&strcmp(A.status,ustatus{il1})&strcmp(A.format_versionTr,uformat{il4})&strcmp(A.modeTr,umode{il5})&~(strcmp(A.status,'D')|strcmp(A.status,'X'))&(A.format_versionTr,'3.1'));
%                     if isempty(ijk)==0
%                     % cree la liste
% %                      if il2==1&il3==1&il1==2&il4==1&il5==1
% %                      %keyboard
% %                      end
%                     length(ijk)
%                         fw1=fopen([WORK.directory name_list],'w')
%                         for ul=1:length(ijk)
%                        % keyboard
%                             if isempty(A.Dep_last_date{ijk(ul)})
%                                strdate=' ';
%                             else
%                                strdate=datestr(A.Dep_last_date{ijk(ul)},'yyyy-mm-dd');
%                             end
%                             fprintf(fw1,'%s\n',[A.platform_number{ijk(ul)} ',' num2str(A.reste_cycle{ijk(ul)}) ',' A.status{ijk(ul)} ',' A.modeTr{ijk(ul)} ',' A.trans_system{ijk(ul)} ',' strrep( A.pi_name{ijk(ul)},',',' ') ',' strrep(A.operating_institution{ijk(ul)},',',' ') ',' strdate ',' datestr(A.Traj_last_date{ijk(ul)},'yyyy-mm-dd') ',' A.dac_format_id{ijk(ul)}]);
%                         end
%                         fclose(fw1);
%                     end
%                 %end
%             %end
%         %end
%     %end
% %end
% end
    
