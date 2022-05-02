 

function [o_alerte3, o_alerte4, o_alerte5, o_alerte6, o_alertCyc_e4,o_alertCyc_e5, idCycprec] = ...
    Test_cycles(a_idLoc, a_cycles, a_cycles_1, a_cycles_sorted, a_duree_cycle,a_mission, a_dureeMedianCycle)


%%------------------------------------------------------------------------------------------------------------
% BOUCLE SUR CHAQUE CYCLE POUR VERIFIER:
% Vérifie: - les double de cycles
%          - pb de numéros de cycle non croissants
%          - mauvaise duree du cycle
%%INPUT:
%% a_idLoc: indice des locs de surface
%% a_cycles: tableau des numéros de cycles
%% a_cycles_1: tableau des cycles théoriques
%% a_cycles_sorted: tableau des cycles tries
%% a_duree_cycle: durée des cycles
%% a_mission: numéro de mission donné
%% a_dureeMedianCycle: mediane de la duree de cycle calculee pour le numero
%% de mission donné.
%% k: indice du flotteur
%%OUTPUT:
%%% o_alerte3,o_alerte4,o_alerte5,o_alerte6: recuperation des flotteurs avec alertes sur les cycles
%%% o_alertCyc_e4: recuperations des cycles concernes par l'alerte 'cycle'.
%%%idCycprec: indices du cycle n-1
%%---------------------------------------------------------------------------------------------------------------

global PARAM;
global T;
global M;
global file_alerte;
global floatname;

cyc=(1);
idCyc=(1);
idCycprec{1}=[];
o_alertCyc_e4 = [];
o_alertCyc_e5 = [];
o_alerte3=[]; o_alerte4=[]; o_alerte5=[];o_alerte6=[];

unique_cyc = unique(T.cycle_number.data); 
[ucyc, occcyc]=unique_withocc(T.cycle_number.data); % add cc 26/10/2020 pour initialiser cyc
maxocc=max(occcyc);
cyc=NaN*zeros(length(a_cycles),maxocc);



numMis = T.config_mission_number.data(1);
%isMis = find(M.config_mission_number == numMis);
isCyc = (T.cycle_number.data(a_idLoc) == a_cycles(1));
idCyc = a_idLoc(isCyc);
cyc(1,1:length(T.juld.data(idCyc)))=T.juld.data(idCyc);  % add cc 18/12/2020 pour pouvoir tester duplication des donnees cycle 0-1 ou 1-2 

% VERIFICATION de la CROISSANCE DES NUMEROS DE CYCLE
 if isequal(a_cycles_sorted,a_cycles)==0   % cc les cycles ne sont pas ranges dans l'ordre croissant            
	o_alerte4=str2double(floatname);
	%o_alertCyc_e4 = [o_alertCyc_e4 a_cycles_sorted(id)];  
	% cc26/10/2020 ne pas faire une alerte sur les cycles
	fid_alerte=fopen(file_alerte,'a');
	fprintf(fid_alerte,'%s\n',[ floatname  ',, warning,CYCLES NUMBERS are not increasing']);
	fclose(fid_alerte);
	fprintf('%s\n',[ floatname  ',,warning, CYCLES NUMBERS are not increasing']);
end

% VERIFICATION DE DUPLICATION des LOCATIONS DATES DANS 2 CYCLES DIFFERENTS
for id=2:length(a_cycles)
    %idCyc_prec=idCyc;
    numCycle = a_cycles(id);
    numCycle_prec = a_cycles(id-1);
    isCyc = (T.cycle_number.data(a_idLoc) == numCycle);
    idCyc = a_idLoc(isCyc);
    numMis = T.config_mission_number.data(id);
    %isMis = find(M.config_mission_number == numMis);
    %idMis = find(double(a_mission)==numMis);
    if(id>1)
        isCycprec = (T.cycle_number.data(a_idLoc)==(numCycle-1));
        idCycprec{id} = [a_idLoc(isCycprec)];
    end
    
    if (id>1 & sum(idCyc)>0 & sum(idCycprec{id})>0)
        
        cyc(id,1:length(T.juld.data(idCyc)))=T.juld.data(idCyc);
        
%        remove cc 26/10/1=2020 car initialisation plus haut        
%         for i = 1:size(cyc,1)
%             a = cyc(i,:);
%             a(find(ismember(cyc(i,:),0)==1))=NaN;    %%%remplace les 0 en Nan pour utilisation ismember.
%             cyc(i,:) = a;
%         end
        
        if(isempty(find(ismember(cyc(1:id-1,:),cyc(id,:))==1))==0 & id>1)
            
            %if ismember(cyc(1:id-1,:),cyc(id,:),'rows')>1 & id >1 % NC
            o_alerte3(a_cycles(id)+1)=str2double(floatname);
            fid_alerte=fopen(file_alerte,'a');
            fprintf(fid_alerte,'%s\n',[floatname ', cycle ' num2str(a_cycles(id)) ',discarded(for u&v), (LOCATION DATES are duplicated in two different cycles']);
            fclose(fid_alerte);
            fprintf('%s\n',[floatname ', cycle ' num2str(a_cycles(id)) ',discarded(for u&v), LOCATION DATES are duplicated in two different cycles']);    %%parfois dû à erreur en fin de batterie, float qui reste en surface
            o_alertCyc_e4 = [o_alertCyc_e4 a_cycles_sorted(id)];
            %T.juld_qc.data(idCyc)=4;
            %T.position_qc.data(idCyc)=4;
            
        end
        
		
        
    end
end           %%%fin de la boucle sur les cycles



% VERIFICATION DU DECOUPAGE DES CYCLES
% modif cc 12/11/2020
% On calcule une duree cumulee (dureeCumulee)  en multipliant la duree du cycle par le (numero de cycle -1)
% Pour chaque cycle on tient compte du numero de mission pour avoir la  duree de cycle correspondante => dureeCycleTh
% On trouve une date "modulo le nombre de cycle"  => dateMod  ie JULD -dureeCumulee
% Sur le diagramme nc_trace_times les dateMod doivent être alignés
% Pour chaque cycle on calcule d'abord une valeur mediane de dateMod => dateModMed, on verifie alors qu'il n'y a pas de date erronée au sein du cycle
% Si c'est le cas elle est flagguée à 6 (o_alertCyc_e5)
% Ensuite on verifie que les dateModMed sont alignés (ou presque ex: derive horloge ). Si ce n'est pas le cas, on identifie les cycles avec un pb de numerotation de cycle

dureeCumulee= zeros(length(a_cycles),1);
dateModMed= NaN*zeros(length(a_cycles),1);
dureeCycleTh= NaN*zeros(length(a_cycles),1);


% calcul en utilisant M.CycleTime comme duree de cycle
a_config_mission_number=M.config_mission_number;
[dureeCumuleeMeta,dateModMedMeta,dureeCycleThMeta,o_alertCyc_e5Meta] = calcule_dureeCumulee(a_cycles,a_cycles_sorted,M.CycleTime,a_config_mission_number,a_idLoc,o_alertCyc_e5);

% calcul en utilisant a_dureeMedianCycle comme duree de cycle
a_config_mission_number=unique(T.config_mission_number.data(~isnan(T.config_mission_number.data)),'stable')';
[dureeCumulee,dateModMed,dureeCycleTh,o_alertCyc_e5] = calcule_dureeCumulee(a_cycles,a_cycles_sorted,a_dureeMedianCycle,a_config_mission_number,a_idLoc,o_alertCyc_e5);

r=floor(length(dateModMedMeta)/2);

% indicateur pour determiner si on va se baser sur a_dureeMedianCycle ou M.CycleTime 
drift=abs(meanoutnan(dateModMed(1:r))-meanoutnan(dateModMed(end-r:end)));
driftMeta=abs(meanoutnan(dateModMedMeta(1:r))-meanoutnan(dateModMedMeta(end-r:end)));

if (driftMeta<drift)|isnan(drift)
dureeCumulee=dureeCumuleeMeta;
dateModMed=dateModMedMeta;
dureeCycleTh=dureeCycleThMeta;
o_alertCyc_e5=o_alertCyc_e5Meta;
end

% Date de loc erronées au sein du cycle
 if isempty(o_alertCyc_e5)==0
 
		     fid_alerte=fopen(file_alerte,'a');
			 fprintf(fid_alerte, '%s\n',[ floatname ', ' num2str(unique(o_alertCyc_e5)) ',warning, SOME LOCATION DATES are suspicious']);
			 fclose(fid_alerte);
			 fprintf('%s\n',[ floatname ', ' num2str(unique(o_alertCyc_e5)) ',warning, SOME LOCATION DATES are suspicious']);

end

% Pb de duree de cycle:
diffdateModMed(1)=NaN;
diffdateModMed(2:length(dateModMed))=dateModMed(2:end)-dateModMed(1:end-1);

isok=abs(diffdateModMed)'<dureeCycleTh*0.5;
isok(a_cycles_sorted==0|a_cycles_sorted==1)=1; % ne tient pas compte des duree de cycle 0 ou 1

 if sum(isok==0)>1
		     fid_alerte=fopen(file_alerte,'a');
			 fprintf(fid_alerte, '%s\n',[ floatname ', ' num2str(a_cycles_sorted(~isok)) ',warning, CYCLE TIME is not consistent with the expected cycle time ']);
			 fclose(fid_alerte);
			 fprintf('%s\n', [ floatname ', ' num2str(a_cycles_sorted(~isok)) ',warning, CYCLE TIME is not consistent with the expected cycle time']);
			 o_alerte6= [o_alerte6 a_cycles_sorted(~isok)];
end

%keyboard

% Pb numerotation de cycle : on cherche la duplication de cycle et les sauts de cycle.
%keyboard
% cycle duplique avec incrementation erronee du numero de cycle
iserrcycl=(a_duree_cycle==0&~isok');
if ~isempty(iserrcycl)&sum(iserrcycl)>0 
		fid_alerte=fopen(file_alerte,'a');
		fprintf(fid_alerte, '%s\n',[ floatname ', ' num2str(a_cycles_sorted(iserrcycl)) ',discarded(for u&v), CYCLE NUMBER could be wrong (duplicate)']);
		fclose(fid_alerte);
		fprintf('%s\n',[ floatname ', ' num2str(a_cycles_sorted(iserrcycl)) ',discarded(for u&v), CYCLE NUMBER could be wrong (duplicate)']);
		o_alertCyc_e4 = [o_alertCyc_e4 a_cycles_sorted(iserrcycl)];
end

% Pour les flotteurs qui n'ont pas plus de deux missions avec des duree de cycle differentes 
% on trouve les cycles dont la duree n'est pas celle attendue: on
% calcule la répartition (hist)
% Si la duree de cycle est peu frequente (occurence <=2) et si elle est multiple de la longueur theorique ; le cycle est eliminé du calcul u &v
if length(unique(M.CycleTime))<=2 
	[m,h]=hist(a_duree_cycle(~isok),[0:1:100]); % repartition
	ll=find(m>=1);                               % indice des duree qui apparaissent au moins une fois
    if isempty(ll)==0                            
        dd=diff(ll)>1;                           % on ne prends pas en compte les durees de cycle avec fortes occurences(>2) et les durees de cycle proche 
        if isempty(dd)==0
            dd=[dd(1) dd]; %
            ik=m(ll)<=2&dd;
            duree_prob=h(ll(ik));                % les durees problematiques sont celles avec une faible occurence et isolée.
        else
            ik=m(ll)<=2;
            duree_prob=h(ll(ik));
        end
        for ill=1:length(duree_prob)
            if duree_prob(ill)>0                  % si les durees problemeatique sont multiples de la duree theorique, on elimine le cycle du calcul de u&v
                iserrcycl= find((a_duree_cycle>=duree_prob(ill)-0.5)&(a_duree_cycle<duree_prob(ill)+0.5));
                isbad=(mod(duree_prob(ill),round(dureeCycleTh(iserrcycl)))==0&duree_prob(ill)~=round(dureeCycleTh(iserrcycl)));
                iserrcycl=iserrcycl(isbad);
                if ~isempty(iserrcycl)&sum(iserrcycl)>0
                    fid_alerte=fopen(file_alerte,'a');
                    fprintf(fid_alerte, '%s\n',[ floatname ', ' num2str(a_cycles_sorted(iserrcycl)) ',discarded(for u&v), CYCLE NUMBER could be wrong ']);
                    fclose(fid_alerte);
                    fprintf('%s\n',[ floatname ', ' num2str(a_cycles_sorted(iserrcycl)) ',discarded(for u&v), CYCLE NUMBER could be wrong']);
                    o_alertCyc_e4 = [o_alertCyc_e4 a_cycles_sorted(iserrcycl)];
                end
            end
        end
    end
end



%uniquement pour des flotteurs mono-mission

% if max(M.CycleTime)<=2   % TODO test à ameliorer
	% if sum(isok)>=1
	% iserrcycl=abs(dateModMed-medianoutnan(dateModMed(isok)))>dureeCycleTh*0.5;
	% else
	% iserrcycl=[];
	% end

	% iserrcycl(1)=0;


	% if ~isempty(iserrcycl)&sum(iserrcycl)>0 
		% fid_alerte=fopen(file_alerte,'a');
		% fprintf(fid_alerte, '%s\n',[ floatname ', ' num2str(a_cycles_sorted(iserrcycl)) ',discarded(for u&v), CYCLE NUMBER could be wrong ']);
		% fclose(fid_alerte);
		% fprintf('%s\n',[ floatname ', ' num2str(a_cycles_sorted(iserrcycl)) ',discarded(for u&v), CYCLE NUMBER could be wrong']);
		% o_alertCyc_e4 = [o_alertCyc_e4 a_cycles_sorted(iserrcycl)];
		%o_alertCyc_e5 = [o_alertCyc_e5 a_cycles_sorted(iserrcycl)];
	% elseif isempty(iserrcycl)
		% fid_alerte=fopen(file_alerte,'a');
		% fprintf(fid_alerte, '%s\n',[ floatname ',warning, Can''t check cycle number ']);
		% fclose(fid_alerte);
		% fprintf('%s\n',[ floatname ',warning, Can''t check cycle number ']);
	% end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dureeCumulee,dateModMed,dureeCycleTh,o_alertCyc_e5] = calcule_dureeCumulee(a_cycles,a_cycles_sorted,a_dureeMedianCycle,a_config_mission_number,a_idLoc,o_alertCyc_e5)
global PARAM;
global T;
global M;
global file_alerte;
global floatname;

dureeCumulee= zeros(length(a_cycles),1);
dateModMed= NaN*zeros(length(a_cycles),1);
dureeCycleTh= NaN*zeros(length(a_cycles),1);
isCyc = find(a_cycles == a_cycles_sorted(1));
numMis = unique(T.config_mission_number.data(isCyc));

if ~isnan(numMis)
    isMis = find(a_config_mission_number == numMis);
    if isempty(isMis)==0
        %dureeCumulee(1)=a_dureeMedianCycle(isMis)*(a_cycles_sorted(1)-1);
        idCyc = a_idLoc(T.cycle_number.data(a_idLoc) == a_cycles_sorted(1));
        
        dureeCumulee(1)=T.juld.data(idCyc(1));
        dureeMod=T.juld.data(idCyc)-dureeCumulee(1);
        % detection de date erronnée au sein du cycle
        iserr=abs(dureeMod-median(dureeMod))>a_dureeMedianCycle(isMis);
        if sum(iserr)>0
            % fid_alerte=fopen(file_alerte,'a');
            % fprintf(fid_alerte, '%s\n',[ floatname ', ' num2str(a_cycles_sorted(1)) ', DATE DE LOC ERRONNEE AU SEIN DU CYCLE ']);
            % fclose(fid_alerte);
            % fprintf('%s\n',[ floatname ', ' num2str(a_cycles_sorted(1)) ', DATE DE LOC ERRONNEE AU SEIN DU CYCLE ']);
            o_alertCyc_e5 = [o_alertCyc_e5 a_cycles_sorted(1)];
        end
        
        
        dateModMed(1)=median(dureeMod);
        dureeCycleTh(1)=a_dureeMedianCycle(isMis);
    else
        dureeCumulee(1)=0;
    end
else
    dureeCumulee(1)=0;
end

for id=2:length(a_cycles_sorted)
	isCyc = find(a_cycles == a_cycles_sorted(id));
	numMis = unique(T.config_mission_number.data(isCyc));
	isMis = find(a_config_mission_number == numMis);
	if ~isnan(numMis)&~isempty(isMis)&length(isMis)==1		
		dureeCumulee(id)=dureeCumulee(id-1)+a_dureeMedianCycle(isMis)*(a_cycles_sorted(id)-a_cycles_sorted(id-1));
		idCyc = a_idLoc(T.cycle_number.data(a_idLoc) == a_cycles_sorted(id));
		dureeMod=T.juld.data(idCyc)-dureeCumulee(id);
		stdDureeMod=std(dureeMod);
		% detection de date erronnée au sein du cycle
		iserr=abs(dureeMod-median(dureeMod))>a_dureeMedianCycle(isMis);
		 if sum(iserr)>0 & stdDureeMod<a_dureeMedianCycle(isMis)
		   % fid_alerte=fopen(file_alerte,'a');
			% fprintf(fid_alerte, '%s\n',[ floatname ', ' num2str(a_cycles_sorted(id)) ', DATE DE LOC ERRONNEE AU SEIN DU CYCLE ']);
			% fclose(fid_alerte);
			% fprintf('%s\n',[ floatname ', ' num2str(a_cycles_sorted(id)) ',DATE DE LOC ERRONNEE AU SEIN DU CYCLE ']);
            o_alertCyc_e5 = [o_alertCyc_e5 a_cycles_sorted(id)];

		 end
		dateModMed(id)=median(dureeMod);	
		dureeCycleTh(id)=a_dureeMedianCycle(isMis);
	end
end

dateModMed(1)=0;
