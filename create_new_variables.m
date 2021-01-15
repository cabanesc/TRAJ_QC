 
function T =  create_new_variables(T, DimT,var); % cree les champs vides pour les vitesses (deep, surface), les erreurs 

varupp = upper(var);
varlow = lower(var);

% create LONGITUDE_DEEP_VELOCITY 
T.(['longitude_' varlow]).name=['LONGITUDE_' varupp];
T.(['longitude_' varlow]).dim={'N_CYCLE'};
T.(['longitude_' varlow]).data=NaN*zeros(DimT.n_cycle.dimlength,1);
T.(['longitude_' varlow]).long_name=['longitude of the location where the ' varlow ' is calculated'];
T.(['longitude_' varlow]).units='degree_east';
T.(['longitude_' varlow]).FillValue_=99999;
T.(['longitude_' varlow]).valid_min=-180;
T.(['longitude_' varlow]).valid_max=180;
T.(['longitude_' varlow]).axis='X';
T.(['longitude_' varlow]).type=6;


% create LATITUDE_DEEP_VELOCITY 
T.(['latitude_' varlow]).name=['LATITUDE_' varupp];
T.(['latitude_' varlow]).dim={'N_CYCLE'};
T.(['latitude_' varlow]).data=NaN*zeros(DimT.n_cycle.dimlength,1);
T.(['latitude_' varlow]).long_name=['latitude of the location where the ' varlow ' is calculated'];
T.(['latitude_' varlow]).units='degree_north';
T.(['latitude_' varlow]).FillValue_=99999;
T.(['latitude_' varlow]).valid_min=-90;
T.(['latitude_' varlow]).valid_max=90;
T.(['latitude_' varlow]).axis='Y';
T.(['latitude_' varlow]).type=6;

% create JULD_DEEP_VELOCITY 
T.(['juld_' varlow]).name=['JULD_' varupp];
T.(['juld_' varlow]).dim={'N_CYCLE'};
T.(['juld_' varlow]).data=NaN*zeros(DimT.n_cycle.dimlength,1);
T.(['juld_' varlow]).long_name=['Julian day (UTC) when ' varlow ' is estimated'];
T.(['juld_' varlow]).units='days since 1950-01-01 00:00:00 UTC';
T.(['juld_' varlow]).conventions= 'Relative julian days with decimal part (as parts of day)';
T.(['juld_' varlow]).resolution= 1.1574e-05;
T.(['juld_' varlow]).FillValue_=999999;
T.(['juld_' varlow]).axis='T';
T.(['juld_' varlow]).type=6;

% create U_DEEP_VELOCITY 
T.(['u_' varlow]).name=['U_' varupp];
T.(['u_' varlow]).dim={'N_CYCLE'};
T.(['u_' varlow]).data=NaN*zeros(DimT.n_cycle.dimlength,1);
T.(['u_' varlow]).long_name=['Estward component of the ' varlow];
T.(['u_' varlow]).units='cm/s';
T.(['u_' varlow]).FillValue_=99999;
T.(['u_' varlow]).type=6;

% create V_DEEP_VELOCITY 
T.(['v_' varlow]).name=['V_' varupp];
T.(['v_' varlow]).dim={'N_CYCLE'};
T.(['v_' varlow]).data=NaN*zeros(DimT.n_cycle.dimlength,1);
T.(['v_' varlow]).long_name=['Northward component of the ' varlow];
T.(['v_' varlow]).units='cm/s';
T.(['v_' varlow]).FillValue_=99999;
T.(['v_' varlow]).type=6;

% create UERR_DEEP_VELOCITY
T.(['uerr_' varlow ]).name=['UERR_' varupp ];
T.(['uerr_' varlow ]).dim={'N_CYCLE'};
T.(['uerr_' varlow ]).data=NaN*zeros(DimT.n_cycle.dimlength,1);
T.(['uerr_' varlow ]).long_name=['Error on the eastward component of the ' varlow];
T.(['uerr_' varlow ]).units='cm/s';
T.(['uerr_' varlow ]).FillValue_=99999;
T.(['uerr_' varlow ]).type=6;

% create V_DEEP_VELOCITY_ERROR
T.(['verr_' varlow ]).name=['VERR_' varupp ];
T.(['verr_' varlow ]).dim={'N_CYCLE'};
T.(['verr_' varlow ]).data=NaN*zeros(DimT.n_cycle.dimlength,1);
T.(['verr_' varlow ]).long_name=['Error on the northward component of the ' varlow];
T.(['verr_' varlow ]).units='cm/s';
T.(['verr_' varlow ]).FillValue_=99999;
T.(['verr_' varlow ]).type=6;
