function [o_alerte11,o_alerte12,isas_alert,isas_non_ref,proche_surface] = Test_PTS_isas(idCyc_drift, id, cycles_sorted,ilong_drift, ilat_drift,i_rpp)

%%-----------------------------------------------------------------------
%% Estime la pression a partir de la temperature mesurée et de la clim ISAS
%% pour chaque mesure en dérive
%% flaggue a mauvais les pressions en alerte
%%INPUT:
%% idcyc_drift: indices de derive en profondeur du cycle n
%% id: iteration du numero de cycle n
%% cycles_sorted: tableau des numeros de cycles tries
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

isfill=((I_temp.temp.data(1,:,ilat_drift,ilong_drift))==I_temp.temp.FillValue_);
I_temp_profile = squeeze(I_temp.temp.add_offset+I_temp.temp.scale_factor.*double(I_temp.temp.data(1,:,ilat_drift,ilong_drift)));
I_temp_profile(isfill)=NaN;

isfill=(I_temp_std.temp_std.data(1,:,ilat_drift,ilong_drift)==I_temp_std.temp_std.FillValue_);
I_temp_std_profile = squeeze(I_temp_std.temp_std.add_offset+I_temp_std.temp_std.scale_factor.*double(I_temp_std.temp_std.data(1,:,ilat_drift,ilong_drift)));
I_temp_std_profile(isfill)=NaN;
id_Notnan_temp_std_profile=find(~isnan(I_temp_std_profile));


isfill=((I_psal.psal.data(1,:,ilat_drift,ilong_drift))==I_psal.psal.FillValue_);
I_psal_profile = squeeze(I_psal.psal.add_offset+I_psal.psal.scale_factor.*double(I_psal.psal.data(1,:,ilat_drift,ilong_drift)));
I_psal_profile(isfill)=NaN;

isfill=(I_psal_std.psal_std.data(1,:,ilat_drift,ilong_drift)==I_psal_std.psal_std.FillValue_);
I_psal_std_profile = squeeze(I_psal_std.psal_std.add_offset+I_psal_std.psal_std.scale_factor.*double(I_psal_std.psal_std.data(1,:,ilat_drift,ilong_drift)));
I_psal_std_profile(isfill)=NaN;
id_Notnan_psal_std_profile=find(~isnan(I_psal_std_profile));

% Definie un min et max de temperature dans la zone pour QC
% elargit la zone au cas ou le flotteur arrive sur de la bathy au cours de
% sa derive
veclat=[max(ilat_drift-1,1):min(ilat_drift+1,size(I_temp.temp.data,3))];
veclon=[max(ilong_drift-1,1):min(ilong_drift+1,size(I_temp.temp.data,4))];

isfill=((I_temp.temp.data(1,:,veclat,veclon))==I_temp.temp.FillValue_);
I_temp_profile_all = squeeze(I_temp.temp.add_offset+I_temp.temp.scale_factor.*double(I_temp.temp.data(1,:,veclat,veclon)));
I_temp_profile_all(isfill)=NaN;
I_temp_profile_all_min = min(min(min(I_temp_profile_all)));
I_temp_profile_all_max = max(max(max(I_temp_profile_all)));

isfill=((I_psal.psal.data(1,:,veclat,veclon))==I_psal.psal.FillValue_);
I_psal_profile_all = squeeze(I_psal.psal.add_offset+I_psal.psal.scale_factor.*double(I_psal.psal.data(1,:,veclat,veclon)));
I_psal_profile_all(isfill)=NaN;
I_psal_profile_all_min = min(min(min(I_psal_profile_all)));
I_psal_profile_all_max = max(max(max(I_psal_profile_all)));

[I_temp_max,idmax]=max(I_temp_profile);
[I_temp_min,idmin]=min(I_temp_profile);

ii=find(I_temp.depth.data(idmax)==I_temp_std.depth.data);
if isempty(ii); ii=length(I_temp_std.depth.data);end
I_temp_max=I_temp_profile_all_max+PARAM.T_N_STD*I_temp_std_profile(ii);
ii=find(I_temp.depth.data(idmin)==I_temp_std.depth.data);
if isempty(ii); ii=length(I_temp_std.depth.data);end
I_temp_min=max(I_temp_profile_all_min-PARAM.T_N_STD*I_temp_std_profile(ii),-2.5);

[I_psal_max,idmax]=max(I_psal_profile);
[I_psal_min,idmin]=min(I_psal_profile);

ii=find(I_psal.depth.data(idmax)==I_psal_std.depth.data);
if isempty(ii); ii=length(I_psal_std.depth.data);end
I_psal_max=I_psal_profile_all_max+PARAM.S_N_STD*I_psal_std_profile(ii);
ii=find(I_psal.depth.data(idmin)==I_psal_std.depth.data);
if isempty(ii); ii=length(I_psal_std.depth.data);end
I_psal_min=max(I_psal_profile_all_min-PARAM.S_N_STD*I_psal_std_profile(ii),-2);




pres_alert=0;
temp_alert=0;
psal_alert=0;
non_ref=1;
proche_surface=0;

temp_mes = T.temp.data(idCyc_drift);
pres_mes = T.pres.data(idCyc_drift);
psal_mes = T.psal.data(idCyc_drift);
temp_noqc4 = T.temp_qc.data(idCyc_drift)<4;
pres_noqc4 = T.pres_qc.data(idCyc_drift)<4;
psal_noqc4 = T.psal_qc.data(idCyc_drift)<4;

%  1ere verification: par rapport à des bornes
% on fait au moins 1 test donc on passe  le flag 0 a 1 pour PRES,TEMP et PSAL %
% si la verif par rapport aux bornes echoue => flag 6
% cc 14/01/2021
temp_qc0 = T.temp_qc.data(idCyc_drift)==0;
pres_qc0 = T.pres_qc.data(idCyc_drift)==0;
psal_qc0 = T.psal_qc.data(idCyc_drift)==0;
T.temp_qc.data(idCyc_drift(temp_qc0))=1;
T.pres_qc.data(idCyc_drift(pres_qc0))=1;
T.psal_qc.data(idCyc_drift(psal_qc0))=1;

bad_temp=(temp_mes<I_temp_min|temp_mes>I_temp_max);
bad_pres=(pres_mes>PARAM.PRESS_PARK_DUMB|pres_mes<-5|(pres_mes==0&temp_mes==0));
bad_psal=(psal_mes<I_psal_min|psal_mes>I_psal_max);

if sum(bad_temp&temp_noqc4)>=1
    T.temp_qc.data(idCyc_drift(bad_temp&temp_noqc4))=6;
    T.psal_qc.data(idCyc_drift(bad_temp&psal_noqc4))=6;
    fid_alerte=fopen(file_alerte,'a');
    fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',BAD TEMPERATURE DETECTED ']);
    fclose(fid_alerte);
    fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',BAD TEMPERATURE DETECTED '])
    temp_alert=1;
end
if sum(bad_pres&pres_noqc4)>=1
    T.pres_qc.data(idCyc_drift(bad_pres&pres_noqc4))=6;
    T.psal_qc.data(idCyc_drift(bad_pres&psal_noqc4))=6;
    fid_alerte=fopen(file_alerte,'a');
    fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',BAD PRESSURE DETECTED ']);
    fclose(fid_alerte);
    fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',BAD PRESSURE DETECTED '])
    pres_alert=1;
end
if sum(bad_psal&psal_noqc4)>=1
    T.psal_qc.data(idCyc_drift(bad_psal&psal_noqc4))=6;
    fid_alerte=fopen(file_alerte,'a');
    fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',BAD SALINITY DETECTED ']);
    fclose(fid_alerte);
    fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',BAD SALINITY DETECTED '])
    psal_alert=1;
end

% % trouve les min et max
ismin = (T.measurement_code.data(idCyc_drift)==297)&(~isnan(T.pres.data(idCyc_drift)));
ismax = (T.measurement_code.data(idCyc_drift)==298)&(~isnan(T.pres.data(idCyc_drift)));

isminmax = ismin|ismax;
ismean = T.measurement_code.data(idCyc_drift)==296;
if pres_mes(ismax)<pres_mes(ismin)
    T.pres_qc.data(idCyc_drift(ismin&pres_noqc4))=6;
    T.pres_qc.data(idCyc_drift(ismax&pres_noqc4))=6;
    bad_pres(isminmax)=1;
end

isstd = T.measurement_code.data==294&T.cycle_number.data==cycles_sorted(id);
presstd =T.pres.data(isstd);

% % on ne verifie la coherence P,T (istocheck=1) a partir d'ISAS:
% % Pour les variables instantannées : oui
% % Pour les variables moyennes : non sauf si la pression du flotteur reste stable pdt la dérive
% % Pour les min/max : non
istocheck(ismean)=1;

% if sum(bad_pres(isminmax))==0&sum(ismin)~=0&sum(ismax)~=0
  % on verifie que la pression moyenne est comprise entre la pression min et la pression max
    % istocheck = ismean&((pres_mes> pres_mes(ismax)|pres_mes<pres_mes(ismin))& bad_pres==0); % indicateur qui sert pour flagguer les means plus tard
% end
 if presstd<100
 istocheck(ismean)=1;
 end
  if presstd>100
 istocheck(ismean)=0;
 end
istocheck(~ismean)=1;


%on verifie les couples (P,T),(P,S) quand c'est possible
for i=1:length(T.temp.data(idCyc_drift))
    if unique(T.cycle_number.data(idCyc_drift))==10
        % keyboard
    end
    temp_mes_i = T.temp.data(idCyc_drift(i));
    pres_mes_i = T.pres.data(idCyc_drift(i));
    psal_mes_i = T.psal.data(idCyc_drift(i));
	
	% isas commence à 1m, si le flotteur est proche de la surface, on cherche la temperature de surface correspondante
	if pres_mes_i<1&pres_mes_i>=-5
	   fid_alerte=fopen(file_alerte,'a');
       fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',Proche Surface, (' num2str(pres_mes_i)  ')']);
       fclose(fid_alerte);
       fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',Proche Surface, (' num2str(pres_mes_i)  ')'])
       pres_mes_i=1;   
       proche_surface=1;	   
	end
    % recuperation temp et psal theorique et temp_std et psal_std
	
    temp_th_i = interp1(I_temp.depth.data,I_temp_profile,pres_mes_i);
	% on  regarde les profiles voisins
	temp_th_i_all=NaN*zeros(size(I_temp_profile_all,2),size(I_temp_profile_all,3));
	for ii=1:size(I_temp_profile_all,2)
	    for jj=1:size(I_temp_profile_all,3)
		    temp_th_i_all(ii,jj) = interp1(I_temp.depth.data,I_temp_profile_all(:,ii,jj),pres_mes_i);
	    end
	end
    temp_std_th_i = interp1(I_temp_std.depth.data,I_temp_std_profile,pres_mes_i);
    if isnan(temp_std_th_i)&~isempty(id_Notnan_temp_std_profile)
        
        if pres_mes_i>I_temp_std.depth.data(id_Notnan_temp_std_profile(end))
            temp_std_th_i= I_temp_std_profile(id_Notnan_temp_std_profile(end));
        else
            temp_std_th_i= I_temp_std_profile(id_Notnan_temp_std_profile(1));
        end
    end
    if  ~isnan(temp_th_i)&~(isminmax(i))
	    non_ref=0;
    end
	if isnan(temp_th_i)&sum(sum(~isnan(temp_th_i_all)))>=1
	   temp_th_i =meanoutnan(meanoutnan(temp_th_i_all));
	   temp_std_th_i = sqrt(temp_std_th_i^2+std(temp_th_i_all(~isnan(temp_th_i_all)))^2);
	   fid_alerte=fopen(file_alerte,'a');
       fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', Temp non referencee: on utilise les plus proches voisins']);
       fclose(fid_alerte);                  
	end
	
    psal_th_i = interp1(I_psal.depth.data,I_psal_profile,pres_mes_i);
	psal_th_i_all=NaN*zeros(size(I_psal_profile_all,2),size(I_psal_profile_all,3));
	for ii=1:size(I_psal_profile_all,2)
	    for jj=1:size(I_psal_profile_all,3)
		    psal_th_i_all(ii,jj) = interp1(I_psal.depth.data,I_psal_profile_all(:,ii,jj),pres_mes_i);
	    end
	end
    psal_std_th_i = interp1(I_psal_std.depth.data,I_psal_std_profile,pres_mes_i);
    if isnan(psal_std_th_i)&~isempty(id_Notnan_psal_std_profile)
        if pres_mes_i>I_psal_std.depth.data(id_Notnan_psal_std_profile(end))
            psal_std_th_i= I_psal_std_profile(id_Notnan_psal_std_profile(end));
        else
            psal_std_th_i= I_psal_std_profile(id_Notnan_psal_std_profile(1));
        end
    end
    if isnan(psal_th_i)&sum(sum(~isnan(psal_th_i_all)))>=1
	   psal_th_i =meanoutnan(meanoutnan(psal_th_i_all));
	   psal_std_th_i = sqrt(psal_std_th_i^2+std(psal_th_i_all(~isnan(psal_th_i_all)))^2);
	end
    
    if ~isnan(temp_mes_i)&~isnan(pres_mes_i)
        if ~isnan(temp_th_i) % valeur referencee dans isas
            % if ismean(i) & presstd>100& presstd<1000
             % delta= 10*PARAM.T_N_STD; % on est plus tolerant sur la coherence quand c'est la moyenne et que le flotteur n'a pas derive a pression ~constante
            % else
            delta= PARAM.T_N_STD;
            %end
            if ~(isminmax(i)) % on ne flaggue pas la non coherence pour les min, max
                if (abs(temp_mes_i-temp_th_i) > delta*temp_std_th_i)
                    
                    if bad_temp(i)==0 & bad_pres(i)==0 &istocheck(i)==1% aucune des donnees n'est vraiment mauvaise par le 1er test, sinon elles ont deja ete flagguees
                        if pres_noqc4(i)==1; T.pres_qc.data(idCyc_drift(i))=6;end;
                        if temp_noqc4(i)==1; T.temp_qc.data(idCyc_drift(i))=6;end;
                        fid_alerte=fopen(file_alerte,'a');
                        fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',P and T are not consistant :  TEMP MESUREE, (' num2str(temp_mes_i) ') ~= TEMP THEORIQUE (' num2str(temp_th_i) ')']);
                        fclose(fid_alerte);
                        fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', P and T are not consistant :  TEMP MESUREE, (' num2str(temp_mes_i) ') ~= TEMP THEORIQUE (' num2str(temp_th_i) ')'])
                        pres_alert=1;
                    end
                    
                else
                    % temp et pres sont consistant
                    if ~(ismean(i))
                        if pres_noqc4(i)==1; T.pres_qc.data(idCyc_drift(i))=1;end;
                        if temp_noqc4(i)==1; T.temp_qc.data(idCyc_drift(i))=1;end;
                        if (abs(psal_mes_i-psal_th_i) > PARAM.S_N_STD*psal_std_th_i)
                            if psal_noqc4(i)==1; T.psal_qc.data(idCyc_drift(i))=6;end;
                            fid_alerte=fopen(file_alerte,'a');
                            fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',P and S are not consistant :  PSAL MESUREE, (' num2str(psal_mes_i) ') ~= PSAL THEORIQUE (' num2str(psal_th_i) ')']);
                            fclose(fid_alerte);
                            fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', P and S are not consistant :  PSAL MESUREE, (' num2str(psal_mes_i) ') ~= PSAL THEORIQUE (' num2str(psal_th_i) ')'])
                            psal_alert=1;
						else
                            if psal_noqc4(i)==1; T.psal_qc.data(idCyc_drift(i))=1;end;
                            
                        end
                    end
                end
            end
        else
            % valeur temp non referencee dans ISAS : soit le float arrive sur de la bathy, resolution grille isas limite, soit pres_mes est dans les choux
            % on va chercher une pression theorique en s'appuyant sur la
            % valeur mesuree de T
            %if ~(ismean(i)|isminmax(i))
            
            
            if ~(isminmax(i))
                fid_alerte=fopen(file_alerte,'a');
                fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', valeur temp non referencee dans ISAS' ]);
                fclose(fid_alerte);
                fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', valeur temp non referencee dans ISAS ']);
                %keyboard
                if bad_temp(i)==0 & bad_pres(i)==0
                    diff_temp = I_temp_profile-temp_mes_i;
                    id_zero = find((sign(diff_temp(1:end-1)).*sign(diff_temp(2:end)))==-1);  % trouve les zeros (-1)
                    
                    dT_dP=diff(I_temp_profile)'./diff(I_temp.depth.data);
					dT_dP(isnan(dT_dP))=meanoutnan(dT_dP);
                    if isempty(id_zero)==0 % temperature referencee
                        int_pres_alert=1; % =1 si une des pressions est mauvaise
                        % trouve les pressions theoriques (plusieurs si inversion de T)
                        for k=1:length(id_zero)
                            pres_th_stdm(k)= interp1(I_temp_profile(id_zero(k):id_zero(k)+1),I_temp.depth.data(id_zero(k):id_zero(k)+1),temp_mes_i);
                            DT(k)= interp1(I_temp_std.depth.data,I_temp_std_profile,pres_th_stdm(k));
                            dT_dPi(k) = interp1((I_temp.depth.data(1:end-1)+I_temp.depth.data(2:end))/2,dT_dP,pres_th_stdm(k));
                            DP(k)=abs(DT(k)./dT_dPi(k));
							
							
                            if abs(pres_mes_i-pres_th_stdm(k))<10*DP(k)
                                int_pres_alert=0;
                            end
                        end
                        if int_pres_alert==1
%                             if pres_noqc4(i)==1; T.pres.data(idCyc_drift(i))=6;end;
%                             fid_alerte=fopen(file_alerte,'a');
%                             fprintf(fid_alerte,'%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ',P and T are not consistant :  PRES MESUREE, (' num2str(pres_mes_i) ')' ]);
%                             fclose(fid_alerte);
%                             fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', P and T are not consistant :  PRES MESUREE, (' num2str(pres_mes_i) ')'])
%                             pres_alert=1;
                          
                        else
                            if pres_noqc4(i)==1; T.pres.data(idCyc_drift(i))=1;end;
                            fprintf('%s\n',[ floatname ', cycle ' num2str(cycles_sorted(id)) ', P and T are  consistant :  PRES MESUREE, (' num2str(pres_mes_i) ')'])
                            
                        end
                    end
                end
            end
        end
        
    end
end


    %keyboard

if pres_alert==1
    if P.Stat==1
        o_alerte11(cycles_sorted(id)+1)=str2double(floatname);
    end
    if(P.map==1)
        m_plot(TEMP_LONG(ilong_drift,1),TEMP_LAT(ilat_drift,1),'color',[0.3010 0.7450 0.9330],'marker','*','markersize',5)
    end
    isas_alert=1;
else
    isas_alert=0;
end
if psal_alert==1
	if P.Stat==1
			o_alerte12(cycles_sorted(id)+1)=str2double(floatname);
		end
		if(P.map==1)
			m_plot(TEMP_LONG(ilong_drift,1),TEMP_LAT(ilat_drift,1),'color',[0.8500 0.3250 0.0980],'marker','s','markersize',8)
	end
end
isas_non_ref=non_ref;
    
   