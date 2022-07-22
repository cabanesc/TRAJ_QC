function parameters = param()

parameters.TIME_STD_DUREE_CYCLE_M = 100;% in day % Max value for Standard Deviation of cycle duration (for more than one cycle)  
parameters.TIME_DIFF_CYCLE = 1; % in day % Max value for difference in cycle duration between META and TRAJ file
parameters.TIME_DIFF_DOUBLE_DATE_LOC = 0; % in hour (initial 30s) % add cc 05/10/2020
parameters.PRESS_PARK_DIFF_BATHY_QC = 100; % Put QC=4 if pressure is xxx dbar below the bathymertry
parameters.PRESS_PARK_DUMB = 2200; % Threshold for dumb park pressure, put QC=4 (cc : ce parametre doit s'adapter si deep => fait ligne 255 Prog_alertes_traj.m)
parameters.PRESS_FOND = 1500; % keep pressure when bathymetry is below
parameters.PRESS_SURF = 50; % keep pressure below
parameters.PRESS_STD_MAX = 2000; % compute pressure std for pressure lower than
parameters.PRESS_STD_MIN = 200; % compute pressure std for pressure greater tan
parameters.PRESS_STD = 200; % Threshold for acceptable STD of parking pressure median estimate
parameters.PRESS_PARK_DIFF_M = 30; % difference of pressure between estimated median pressure and metadata
parameters.PRESS_PARK_DIFF_ISAS = 100; % difference of pressure between estimated pressure from ISAS temperature and estimated from measurements
parameters.PRESS_PARK_DIFF_BATH= 30; % difference of measured parking pressure and bathymetry
parameters.TIME_DIFF_2LOCS = 10; % in hour % Threshold for acceptable difference between successive location
parameters.TIME_DIFF_FLLOCS_AR = 24; % in hour % Thershold for acceptable difference between first and last location for Argos transmission system (within one cycle)
parameters.TIME_DIFF_FLLOCS_IR = 10; % in hour % Thershold for acceptable difference between first and last location for Iridium transmission system (within one cycle)
parameters.TIME_DIFF_FLLOCS_EOL = 60; % in hour % Threshold for assuming EOL (float End Of Life)
%parameters.TIME_LAUNCH_FIRST_LOC = -3/24; % in day % Difference between launch time and first loc
parameters.TIME_LAUNCH_FIRST_LOC = 3; % in hour % Difference between launch time and first loc  % correction cc 15/09/2020 mettre tous en heure
parameters.TIME_LAUNCH_FIRST_LOC_DPF = 10.5; % in hour % Threshold for assuming DPF (Deep Profile First)
parameters.TIME_LAUNCH_FIRST_LOC_DUMB = 100; % on hour % Threshold for assuming dumb difference between launch time and first location 
parameters.ALT_GROUND_NEAR_COAST = -200;  % Min depth for loching coastal grounding
parameters.SPEED_MAX = 3; %m/s % Thershold for aberant deep drift velocity

% Test_TS.m
parameters.T_N_STD = 10.; % N STD for ISAS/in drift measurement temperature comparison 
parameters.S_N_STD = 50.; % N STD for ISAS/in drift measurement salinity comparison

% Test_cycle.m
%parameters.TIME_DUREE_CYCLE_M = 0.15; % Threshold for comparison between estimate median cycle duration and metadata
parameters.TIME_DUREE_CYCLE_M = 90; % Threshold (in %) for comparison between estimate median cycle duration and the actual cycle duration  correction cc 15/09/2020
parameters.TIME_DUREE_CYCLE_JUMP = 0.5; % Threshold for assuming cycle jump/mission cycle

% Test_locs.m
parameters.LOC_LAUNCH_DIFF_M = 0.1; % Threshold for comparison between traj launch location and metadata 
parameters.DATE_LAUNCH_DIFF_M = 1e-4; % Threshold for comparison between traj launch date and metadata -minute près % add cc 18/12/2020


end
