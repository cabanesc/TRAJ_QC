function [M,Mf] = read_file_meta_nc(meta_fileName)

% read all the variables
[Mf,DimM] = read_netcdf_allthefile(meta_fileName);
Mf=replace_fill_bynan(Mf);


M.config_mission_number = Mf.config_mission_number.data;
launch_config_parameter = cellstr(Mf.launch_config_parameter_name.data);
config_parameter = Mf.config_parameter_name.data;

% look for configuration parameters and values
configParam = {'CycleTime','DownTime','UpTime','ParkPressure','ProfilePressure','DeepProfileFirstFloat'};

% *** CycleTime **** :For APEX and ARVOR floats this is the total duration of one cycle, usually 240 hours (10 days).
% For SOLO floats this is the total duration of one cycle, assuming that all float operations reach their full time-out
% intervals before moving to the next float stage. Typically the actual cycle time will be shorter than this value.
% *** DownTime ****: preset duration of the down time (start of descent to start of ascent) APEX, NINJA
% *** UpTime ****: preset duration of the 'up' time (start of ascent to start of
% descent  APEX, NINJA
% *** DeepProfileFirstFloat ****: Some APEX floats are set to do a "profile on deployment". This means that when the float is deployed
% it skips the drift phase at the PARKING depth and sinks directly to the PROFILE depth and starts
% ascending for the first profile. The result is that the first cycle is of shorter duration with a profile depth
% equal to PROFILE configuration pressure (Yes = 1, No = 0).



%--------------------------------------------------------------------------------------------------------------
% - first look for launch configuration parameters (this includes the subset of parameters that do not change for the life of the float.
units={'days','hours','minutes','secondes'};
ratio=[1,24,1440,86400];

for iparam = 1:3
    % initialize
    M.(configParam{iparam}) = NaN*zeros(DimM.n_missions.dimlength,1);
    
    for k=1:length(units)
        idF = find(findstr_tab(launch_config_parameter,['CONFIG_' configParam{iparam} '_' units{k}]));
        if isempty(idF)==0 & Mf.launch_config_parameter_value.data(idF)>0
            M.(configParam{iparam}) = repmat(Mf.launch_config_parameter_value.data(idF),DimM.n_missions.dimlength,1)./ratio(k);
        end
    end
end

units={'dbar','cbar','mbar','bar'};
ratio=[1,10,100,0.1];

for iparam = 4:5
    % initialize
    M.(configParam{iparam}) = NaN*zeros(DimM.n_missions.dimlength,1);
    
    for k=1:length(units)
        idF = find(findstr_tab(launch_config_parameter,['CONFIG_' configParam{iparam} '_' units{k}]));
        if isempty(idF)==0 & Mf.launch_config_parameter_value.data(idF)>0
            M.(configParam{iparam}) = repmat(Mf.launch_config_parameter_value.data(idF),DimM.n_missions.dimlength,1)./ratio(k);
        end
    end
end

iparam = 6;
% initialize
M.(configParam{iparam}) = NaN*zeros(DimM.n_missions.dimlength,1);
idF = find(findstr_tab(launch_config_parameter,['CONFIG_' configParam{iparam} '_LOGICAL']));
if isempty(idF)==0 & Mf.launch_config_parameter_value.data(idF)
    M.(configParam{iparam}) = repmat(Mf.launch_config_parameter_value.data(idF),DimM.n_missions.dimlength,1)./ratio(k);
end
   
%--------------------------------------------------------------------------------------------------------------
% - then look for configuration parameters (this includes the subset of parameters that change for each mission )
units={'days','hours','minutes','secondes'};
ratio=[1,24,1440,86400];

for iparam = 1:3
    for k=1:length(units)
        idF = find(findstr_tab(config_parameter,['CONFIG_' configParam{iparam} '_' units{k}]));
        if isempty(idF)==0 & Mf.config_parameter_value.data(:,idF)>0
            M.(configParam{iparam}) = Mf.config_parameter_value.data(:,idF)./ratio(k);
        end
    end
end

units={'dbar','cbar','mbar','bar'};
ratio=[1,10,100,0.1];

for iparam = 4:5   
    for k=1:length(units)
        idF = find(findstr_tab(config_parameter,['CONFIG_' configParam{iparam} '_' units{k}]));
        if isempty(idF)==0&Mf.config_parameter_value.data(:,idF)>0
            M.(configParam{iparam}) = Mf.config_parameter_value.data(:,idF)./ratio(k);
        end
    end
end

iparam = 6;
idF = find(findstr_tab(config_parameter,['CONFIG_' configParam{iparam} '_LOGICAL']));
if isempty(idF)==0&Mf.config_parameter_value.data(:,idF)>0
    M.(configParam{iparam}) = Mf.config_parameter_value.data(:,idF)./ratio(k);
end
   
%--------------------------------------------------------------------------------------------------------------
% Finalize the output with other useful variables
ii=isnan(M.CycleTime);
M.CycleTime(ii)=M.DownTime(ii)+M.UpTime(ii);


M.config_mission_number=M.config_mission_number';
M.CycleTime = M.CycleTime';
M.DownTime = M.DownTime';
M.UpTime = M.UpTime';
M.ParkPressure = M.ParkPressure';
M.ProfilePressure = M.ProfilePressure';
M.DeepProfileFirstFloat = M.DeepProfileFirstFloat';
M.platform_number=strtrim(Mf.platform_number.data)';
M.trans_system=strtrim(Mf.trans_system.data);
M.platform_type=strtrim(Mf.platform_type.data)';
M.launch_date=strtrim(Mf.launch_date.data)';
M.launch_latitude=Mf.launch_latitude.data;
M.launch_longitude=Mf.launch_longitude.data;
M.launch_qc=Mf.launch_qc.data;
M.pi_name=strtrim(Mf.pi_name.data)';
M.dac_format_id=strtrim(Mf.dac_format_id.data)';
M.controller_board_type_primary=strtrim(Mf.controller_board_type_primary.data)';


