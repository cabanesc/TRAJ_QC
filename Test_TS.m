function [o_alerte11, o_alerte12, temp_alert, temp_non_ref] = Test_TS(idCyc_drift, id, cycles_sorted,pres_long,pres_lat,numMis, missions,i_rpp, idepth_drift, ...
    ilong_drift, ilat_drift, idepth_std_drift)

%%-----------------------------------------------------------------------
%% realise des comparaisons entre la temperature et la salinite mesuree et
%% celle estimee par la climatologie et alerte si pb observ�.
%% Recupere les cycles concernes par les alertes et les param�tres temp_alert et temp_non_ref.
%%
%%INPUT:
%% idcyc_drift: indices de d�rive en profondeur du cycle n
%% id: it�ration du num�ro de cycle n
%% cycles_sorted: tableau des num�ros de cycles tri�s
%% pres_long, pres_lat: longitude et latitude
%% numMis: num�ro de misions donn�
%% missions: tableaux des num�ros de mision
%% i_rpp: indice des mesures � la profondeur de parking

%% OUTPUT:
%% o_alerte11: cycles concern�s par l'alerte sur la temperature
%% o_alerte12: cycles concernes par l'alerte sur la salinite
%% temp_alert: valeur 0 ou 1. 1 si pb d�tect� sur entre la temp�rature mesur�e et la temp�rature th�orique
%% temp_non_ref: valeur 0 ou 1. 1 si la temp�rature n'est pas r�f�renc�e dans la climatologie
%%--------------------------------------------------------------------
global PARAM;
global T;
global P;
global file_alerte;
global floatname;
global TEMP_LONG;
global TEMP_LAT;
global TEMP_TEMP_STD;
global PSAL_PSAL_STD;
global TEMP_TEMP;
global PSAL_PSAL;
global I_psal_std;
global I_temp_std;
global I_psal;
global I_temp;


o_alerte11 = []; o_alerte12 = [];

%temp_hold=50;
%psal_hold=50;
itemp=0;
ipsal=0;
% recupeation de la temperature et de la salinite a la
% premiere position du flotteur durant le cycle


temp_drift_mes(i_rpp)=mynanmean(T.temp.data(idCyc_drift));
psal_drift_mes(i_rpp)=mynanmean(T.psal.data(idCyc_drift));
%pres_drift_mes(i_rpp)=abs(mynanmean(T.pres.data(idCyc_drift)));

% conversion d'echelle de la temperature et de la salinite
% mesurees
temp=I_temp.temp.add_offset+I_temp.temp.scale_factor.*double(TEMP_TEMP(:,idepth_drift,ilat_drift,ilong_drift));
temp_std=I_temp_std.temp_std.scale_factor*double(TEMP_TEMP_STD(:,idepth_std_drift,ilat_drift,ilong_drift));
psal=I_psal.psal.add_offset+I_psal.psal.scale_factor.*double(PSAL_PSAL(:,idepth_drift,ilat_drift,ilong_drift));
psal_std=I_psal_std.psal_std.scale_factor*double(PSAL_PSAL_STD(:,idepth_std_drift,ilat_drift,ilong_drift));

temp_drift_th(i_rpp)=temp;
psal_drift_th(i_rpp)=psal;
idMis = find(double(missions)==numMis);
%temp_drift_th2(id)=temp_drift_th(i_rpp);    %%si l'on veut r�cup�rer la variable pour plotter ensuite


%%%Verification de la temp mesuree

if (temp_drift_mes(i_rpp)<temp_drift_th(i_rpp)-PARAM.T_N_STD*temp_std ...
        | temp_drift_mes(i_rpp)>temp_drift_th(i_rpp)+PARAM.T_N_STD*temp_std)&temp_drift_th(i_rpp)<50
    
    itemp=itemp+1;
    cyc_wrong_temp(itemp)=i_rpp;
    temp_wrong(itemp)=temp_drift_mes(i_rpp);
    fid_alerte=fopen(file_alerte,'a');
    fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB TEMPERATURE. TEMP MESUREE (' num2str(temp_drift_mes(i_rpp)) ') ~= TEMP THEORIQUE (' num2str(temp_drift_th(i_rpp)) ')']);
    fclose(fid_alerte);
    fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB TEMPERATURE. TEMP MESUREE (' num2str(temp_drift_mes(i_rpp)) ') ~= TEMP THEORIQUE (' num2str(temp_drift_th(i_rpp)) ')'])
    if P.Stat==1
        alerte11(cycles_sorted(id)+1)=str2double(floatname);
    end
    if(P.map==1)
        m_plot(TEMP_LONG(ilong_drift,1),TEMP_LAT(ilat_drift,1),'color',[0.3010 0.7450 0.9330],'marker','*','markersize',5)
    end
    temp_alert(i_rpp)=1;
    temp_non_ref(i_rpp)=0;
elseif temp_drift_th(i_rpp)>50
    fid_alerte=fopen(file_alerte,'a');
    fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', TEMPERATURE NON REFERENCEE']);
    fclose(fid_alerte);
    fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', TEMPERATURE NON REFERENCEE']);
    temp_non_ref(i_rpp)=1;
    temp_alert(i_rpp)=0;
else
    temp_alert(i_rpp)=0;
    temp_non_ref(i_rpp)=0;
end

% Verification de la salinite mesuree
if (psal_drift_mes(i_rpp)<psal_drift_th(i_rpp)-PARAM.S_N_STD*psal_std ...
        | psal_drift_mes(i_rpp)>psal_drift_th(i_rpp)+PARAM.S_N_STD*psal_std)&psal_drift_th(i_rpp)<60
    ipsal=ipsal+1;
    % recuperation des donnees aberrantes
    cyc_wrong_psal(ipsal)=i_rpp;
    psal_wrong(ipsal)=psal_drift_mes(i_rpp);
    fid_alerte=fopen(file_alerte,'a');
    fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB SALINITE. SAL MESUREE (' num2str(psal_drift_mes(i_rpp)) ') ~= SAL THEORIQUE (' num2str(psal_drift_th(i_rpp)) ')']);
    fclose(fid_alerte);
    fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', PB SALINITE. SAL MESUREE (' num2str(psal_drift_mes(i_rpp)) ') ~= SAL THEORIQUE (' num2str(psal_drift_th(i_rpp)) ')']);
    if P.Stat==1
        alerte12(cycles_sorted(id)+1)=str2double(floatname);
    end
    if(P.map==1)
        m_plot(TEMP_LONG(ilong_drift,1),TEMP_LAT(ilat_drift,1),'color',[0.8500 0.3250 0.0980],'marker','s','markersize',8)
    end
elseif psal_drift_th(i_rpp)>60
    fid_alerte=fopen(file_alerte,'a');
    fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', SALINITE NON REFERENCEE']);
    fclose(fid_alerte);
    fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', SALINITE NON REFERENCEE']);
end


end
