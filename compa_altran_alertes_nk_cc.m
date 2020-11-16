%%%%COMPARAISON alerteS ALTRAN ET alerteS TEST

clear all
close all

id=0;

param=config;

addpath(param.DIR_TOOL);
addpath([param.DIR_TOOL '/Lib_Argo/']);
addpath([param.DIR_TOOL '/Lib_Argo/Plots']);
addpath([param.DIR_TOOL 'Lib_Argo/Tests_TR']);
addpath([param.DIR_TOOL '/Lib_Argo/RWnetcdf/R2008b']);
liste_rep_alerte = param.liste_rep_alerte;

liste_alertes = param.liste_alertes;

%%%%Fichiers Altran
DIR = [param.DIR_LISTE 'Liste_livraison2019/'];

%liste_al = dir([DIR '*.txt']);   %%%liste flotteurs altran

%%% Fichiers test
%DIR2 = '/home/gherbert/andro/andro/soft/solene/Erreurs/Erreurs_b/';
DIR2 = [param.DIR_HOME 'alerts/'];
%DIR2ok = [DIR2 liste_rep_alerte{i}]
%liste2 = dir([DIR2ok '/alertes_*']);    %%liste résultats alertes test


%%% Récupération des fichiers (liste initiale, liste alerte altran, liste
%%% alerte test)

%for i = 1:length(liste_rep_alerte)

nb_cyc_id_all = 0; nb_cyc_sup_all = 0; nb_cyc_miss_all = 0;
nb_flo_id_all =0; nb_flo_sup_all=0; nb_flo_miss_all= 0;
nb_flo_al_all = 0; nb_flo_test_all = 0;
flo_total_all = 0; flo_without_err_al_all = 0; flo_without_err_test_all = 0;


for l = 1:length(liste_rep_alerte)     %%%  boucle sur les listes
    
    fich_allfloat = [DIR liste_rep_alerte{l} '/' liste_rep_alerte{l} '.txt'];
    
    alerte={};
    idi=0;
    float_test=[];float_testok=[];
    float_al=[];float_alok=[];
    float_all=[];
    cyc_al = [];
    cyc_test = [];
    
	
	% Recuperation de la liste de tous les flotteurs de la liste
	mat_all = importdata(fich_allfloat);
	for k = 1:length(mat_all)
		if isa(mat_all,'cell')
			float_all(k) = str2num(mat_all{k}(1:7));  %%%récupération flotteur initiaux
		elseif isa(mat_all,'double')
			float_all(k) = mat_all(k);
		end
		
	end
	
    for j = 1:length(liste_alertes)     %%% boucle sur les types d'erreurs
        fich_err = [DIR liste_rep_alerte{l} '/' liste_rep_alerte{l} '_erreur_' liste_alertes{j} '.txt'];
		
		
        if exist(fich_err,'file')==2
		    idi = idi+1;
			alerte{idi} = liste_alertes{j};
			
		    err_al = importdata(fich_err,',');    %%% Erreurs altran
            
            fich_err_t = dir([DIR2 liste_rep_alerte{l} '/*' liste_alertes{j} '*b.txt']);
            fich_err_test = [DIR2 liste_rep_alerte{l} '/' fich_err_t.name];
            
            err_test = importdata(fich_err_test,',');      %%%Erreurs du test
			
			
			% récupération des flotteurs/cycles concernes par l'alerte Altran
            id=0;
            for k = 1:length(err_al)      %%% k: flotteurs concernés par l'erreur d'apres altran
                if isa(err_al,'cell')
                    float_al(idi,k) = str2num(err_al{k}(1:7));  %%%récupération flotteur alertes altran  idi: type erreur   k: flotteurs
                elseif isa(err_al,'double')
                    float_al(idi,k) = err_al(k);
                end
				
				% verification que ces cycles  existent bien dans les traj  add cc 03/11/2020
				NcVar.cycle_number.name='CYCLE_NUMBER'; 
				doneT=0;
				ii=0;
				while ~doneT&& ii<=length(param.DAC)-1
					% lecture des traj
					ii=ii+1;
					floatname =num2str(float_al(idi,k));
					traj_fileName_R = [param.DIR_FTP param.DAC{ii} '/' floatname '/' floatname '_Rtraj.nc'];
					traj_fileName_D = [param.DIR_FTP param.DAC{ii} '/' floatname '/' floatname '_Dtraj.nc'];
					if exist(traj_fileName_D,'file');
						[T,DimT,GlobT]=read_netcdf_allthefile(traj_fileName_D,NcVar);
						traj_fileName_final=[floatname '_Dtraj.nc'];
						doneT=1;
					elseif exist(traj_fileName_R,'file');
						[T,DimT,GlobT]=read_netcdf_allthefile(traj_fileName_R,NcVar);
						doneT=1;
						traj_fileName_final=[floatname '_Rtraj.nc'];
					end
				end 
				cycle_traj=unique(T.cycle_number.data); 
				
                id = id+1;
                if(j~=7)      %%%fichier launchdate où pas de cycle indiqué
                    if isa(err_al,'cell')
                        a=strfind(err_al{k},'['); b = strfind(err_al{k},']');
						
                        thecycles = [str2num(err_al{k}(a+1:b-1))];
						isintraj = ismember(thecycles,cycle_traj);
						cyc_al{idi,id} = thecycles(isintraj);
						
                    elseif isa(err_al, 'double')
						thecycles = err_al(length(err_al)+k); 
						isintraj = ismember(thecycles,cycle_traj);
						cyc_al{idi,id} = thecycles(isintraj);
                    end
                else
                    cyc_al{idi,id} = 1;
                end
				
            end
            
            % récupération des flotteurs/cycles concernes par l'alerte Test
            id=0;
            %  keyboard
            for k = 1:length(err_test)
                if isa(err_test,'cell')
                    float_test(idi,k) = str2num(err_test{k}(1:7));  %%%récupération flotteur alertes test
                elseif isa(err_test,'double')
                    float_test(idi,k) = err_test(k);
                end
                
                id = id+1;
                if(j~=7)      %%%fichier launchdate où pas de cycle indiqué
                    if isa(err_test,'cell')
                        a=strfind(err_test{k},'['); b = strfind(err_test{k},']');
                        cyc_test{idi,id} = [str2num(err_test{k}(a+1:b-1))];
                    elseif isa(err_test, 'double')
                        cyc_test{idi,id} = err_test(length(err_test)+k);
                    end
                else
                    cyc_test{idi,id} = 1;
                end
            end
            if isempty(err_test)     %%cas où le fichier erreur test existe mais il est vide
                float_test(idi,k)  = err_test(k);
            end
            
        end
    end
    
    
    % keyboard
    %%%% Récupération des cycles pour un flotteur donné (car répétition de meme
    %%%% flotteur dans la liste altran).
    
    
    cyc_al_ok=[];
    cycalok=[];
    
    for i = 1:size(float_al,1)
        id=1;
        for m = 2:size(float_al,2)
            
            if(float_al(i,m)~=float_al(i,m-1))
                if(m==2)
                    cycalok{i,1}=cyc_al{i,1};
                end
                id = id+1;
                cycalok{i,id} = cyc_al{i,m};
            end
            
            if(float_al(i,m)==float_al(i,m-1))
                if(m==2)
                    cycalok{i,1} = [cyc_al{i,1} cyc_al{i,2}];
                else
                    cyc_al_ok = [cycalok{i,id} cyc_al{i,m}];
                    cycalok{i,id} = cyc_al_ok;
                end
            end
            %keyboard
            
        end
    end
    
    %%%idem test
    
    %cyc_testok=[];
    %cyctestok=[];
    
    %for i = 1:size(float_test,1)
    %   id=1;
    
    
    %   for j = 2:size(float_test,2)
    
    %       if(float_test(i,j)~=float_test(i,j-1))
    %           if(j==2)
    %              cyctestok{i,1}=cyc_test(i,1);
    %           end
    %           id = id+1;
    %           cyctestok{i,id} = cyc_test(i,j);
    %       end
    
    %       if(float_test(i,j)==float_test(i,j-1))
    
    %           if(j==2)
    %              cyctestok{i,1} = [cyc_test(i,1) cyc_test(i,2)];
    %           else
    %              cyc_testok = [cyctestok{i,id} cyc_test(i,j)];
    %              cyctestok{i,id} = cyc_testok;
    %           end
    %       end
    %
    %   end
    %end
    
    
    
    
    %%% met  à NaN les 0
    
    float_al(find(float_al==0))=NaN;
    float_test(find(float_test==0))=NaN;
    %float_all(find(float_all==0))=NaN;
    for i = 1:size(float_test,1)
        float_testok{i} = unique(float_test(i,:),'stable');
        float_alok{i} = unique(float_al(i,:),'stable');
    end
    
    float_id = [];
    float_testsup = [];
    nb_idem=[];
    nb_diff=[];
    nb_sup = [];
    flo_identique = {};
    flo_supp = {};
    flo_miss = {};
    nb_flo_al = [];
    nb_flo_test = [];
    nb_flo_supp = [];
    nb_flo_miss = [];
    nb_flo_id = [];
    Prcent_flo_id = [];
    Prcent_flo_supp = [];
    Prcent_flo_miss = [];
	
    cyc_test_identique =cell(1,1);
    cyc_test_miss = cell(1,1);;
    cyc_test_supp = cell(1,1);;
    nb_cyc_id= zeros(1,length(float_testok));
    nb_cyc_supp= zeros(1,length(float_testok));
    nb_cyc_miss= zeros(1,length(float_testok));
    nb_flo_id= zeros(1,length(float_testok));
    nb_flo_supp= zeros(1,length(float_testok));
    nb_flo_miss= zeros(1,length(float_testok));
    Prcent_cyc_id= zeros(1,length(float_testok));
    Prcent_cyc_supp(i) =zeros(1,length(float_testok));
	Prcent_cyc_miss(i)=zeros(1,length(float_testok));
   
    for i = 1:length(float_testok)   %%%pour chaque alerte
 
        ii=0;
        iii=0;
        it=0;
        it2=0;
        flo_al = double(float_alok{i}(:));    %%%flotteurs altran pour l'alerte donnée
        flo_test = double(float_testok{i}(:));  %%%flotteurs test pour l'alerte donnée
        
        flo_al= flo_al(~isnan(flo_al)); flo_test=flo_test(~isnan(flo_test));
        [is_identique] = ismember(flo_test,flo_al);
        is_identique_al = ismember(flo_al,flo_test);
        
        flo_identique{i} = flo_test(is_identique);   %%%flotteur identiques dans Altran et test
        flo_supp{i} = flo_test(~is_identique);      %%flotteurs supplémentaire dans test
        flo_miss{i} = flo_al(~is_identique_al);      %%%flotteurs manquants dans test
        
        nb_flo_total = length(float_all);
        nb_flo_id(i) = length(find(~isnan(flo_identique{i})));
        Prcent_flo_id(i) = nb_flo_id(i)*100/nb_flo_total;
        nb_flo_supp(i) = length(find(~isnan(flo_supp{i})));
        Prcent_flo_supp(i) = nb_flo_supp(i)*100/nb_flo_total;
        nb_flo_miss(i) = length(find(~isnan(flo_miss{i})));
        Prcent_flo_miss(i) = nb_flo_miss(i)*100/nb_flo_total;
        
        nb_flo_al(i) = length(find(~isnan(float_al(i,:))));
        nb_flo_test(i) = length(find(~isnan(float_test(i,:))));
        
        nb_without_err_al = length(find(ismember(float_all,float_al(:,:))==0));
        nb_without_err_test = length(find(ismember(float_all,float_test(:,:))==0));
        
        %%%%Cycles identiques, supplémentaires et manquants pour un même
        %%%%flotteur  % correction cc 03/11/2020 : les cycles supp/missing pour des flotteurs identiques dans les listes altran et TRAJQC n'etaient pas comptés
		%%%%%
        float_al_test = unique([flo_al ;flo_test]); % tous les flotteurs altran/test concernes par l'alerte
		mm=0;
		for ifloat=1:length(float_al_test)
		    % on verifie si flotteur present dans flo_al
			isal=ismember(flo_al,float_al_test(ifloat));
		    istest=ismember(flo_test,float_al_test(ifloat));
			 mm=mm+1;
			if sum(isal)==1&&sum(istest)==1 % flotteur dans les deux listes
			
			   cyc_test_i_float = unique(cyc_test{i,istest});
               cyc_al_i_float = unique(cycalok{i,isal});
			   cyc_test_identique{i,mm}= [cyc_test_i_float(ismember(cyc_test_i_float,cyc_al_i_float))];
			   nb_cyc_id(i) = nb_cyc_id(i) + sum(ismember(cyc_test_i_float,cyc_al_i_float));
			   cyc_test_supp{i,mm}= [cyc_test_i_float(~ismember(cyc_test_i_float,cyc_al_i_float))];
			   nb_cyc_supp(i) = nb_cyc_supp(i) + sum(~ismember(cyc_test_i_float,cyc_al_i_float));
			   cyc_test_miss{i,mm}= [ cyc_al_i_float(~ismember(cyc_al_i_float,cyc_test_i_float))];
			   nb_cyc_miss(i) = nb_cyc_miss(i) + sum(~ismember(cyc_al_i_float,cyc_test_i_float));			   
           
            else
			cyc_test_supp{i,mm}=[];
			cyc_test_miss{i,mm}=[];
			end
            if sum(isal)==1&&sum(istest)==0 % flotteur seulement dans alertes altran
               cyc_al_i_float = unique(cycalok{i,isal});
			   cyc_test_miss{i,mm}= [cyc_test_miss{i,mm} cyc_al_i_float];
			   nb_cyc_miss(i) = nb_cyc_miss(i) + length(cyc_al_i_float);			   
            end	
            if sum(isal)==0&&sum(istest)==1 % flotteur seulement dans alertes test
			   cyc_test_i_float = unique(cyc_test{i,istest});
			   cyc_test_supp{i,mm}= [cyc_test_supp{i,mm} cyc_test_i_float];
			   nb_cyc_supp(i) = nb_cyc_supp(i) + length(cyc_test_i_float);			   
            end
        end	
		Prcent_cyc_id(i) = nb_cyc_id(i)./(nb_cyc_id(i)+nb_cyc_miss(i))*100;
		Prcent_cyc_supp(i) = nb_cyc_supp(i)./(nb_cyc_id(i)+nb_cyc_miss(i))*100;
		Prcent_cyc_miss(i) = nb_cyc_miss(i)./(nb_cyc_id(i)+nb_cyc_miss(i))*100;
		
        % %%%cycles identiques dans test et altran et comptage
        % if(isempty(find(is_identique==1))==0)
            
            % ident = find(is_identique==1); % trouve les flotteurs identiques dans les deux listes, index test
            % ident_al = find(is_identique_al==1); % trouve les flotteurs identiques dans les deux listes, index altran
			
            % for mm = 1:length(ident)% pour tous les flotteurs identique, il faut trouver les cycles identiques entre altran et test% cc 03/11/2020
			    % cyc_test_i_float = cyc_test{i,ident(mm)};
				% cyc_al_i_float=cycalok{i,ident_al(mm)};
                % cyc_test_identique{i,mm}= cyc_test{i,ident(mm)};
                % nb_cyc_id(i) = nb_cyc_id(i) + length(double(cyc_test_identique{i,mm}));
            % end
        % else
            % cyc_test_identique{i} = [];
        % end
        
        % %%%cycle supplémentaires dans test (et comptage)
        % if(isempty(find(is_identique==0))==0)
            
            % non_ident = find(is_identique==0);
            
            % for mm = 1:length(non_ident)
                % cyc_test_supp{i,mm}= cyc_test{i,non_ident(mm)};
                % nb_cyc_supp(i) = nb_cyc_supp(i)+length(double(cyc_test_supp{i,mm}));
            % end
        % else
            % cyc_test_supp{i} = [];
        % end
        
        % %%cyles manquants dans test et comptage
        % %  keyboard
        % if(isempty(find(is_identique_al==0))==0)
            % miss = find(is_identique_al==0);
            % for u = 1:length(miss)
                % cyc_test_miss{i,u} = cell2mat(cycalok{i,miss(u)});
                % nb_cyc_miss(i) = nb_cyc_miss(i) + length(double(cyc_test_miss{i,u}));
            % end
            
            % %else
            % %  cyc_test_miss{i} = [];
        % end
        
        
        
        
    end
    
    %%%figures
    nb_flo=zeros(length(Prcent_flo_id),2);
    Prcent=zeros(length(Prcent_flo_id),3);
    nb_cyc=zeros(length(Prcent_flo_id),3);
    nb_flot=zeros(length(Prcent_flo_id),3);
    
    
    
    for ii = 1:length(Prcent_flo_id)  %%pour chaque erreur
        Prcent(ii,:) = [Prcent_flo_id(ii) Prcent_flo_supp(ii) Prcent_flo_miss(ii)];
    end
    
    for ii = 1:length(nb_flo_id)   %%pour chaque erreur
        nb_flot(ii,:) = [nb_flo_id(ii) nb_flo_supp(ii) nb_flo_miss(ii)];
    end
    
    
    for ii = 1:length(nb_cyc_id)   %%pour chaque erreur
        nb_cyc(ii,:) = [nb_cyc_id(ii) nb_cyc_supp(ii)  nb_cyc_miss(ii)];
    end
    
    for ii = 1:length(nb_flo_al)    %%pour chaque erreur
        nb_flo(ii,:) = [nb_flo_al(ii) nb_flo_test(ii)]
    end
    
    flo_tot = [nb_flo_total nb_without_err_al nb_without_err_test];
    
    
    nb_cyc_id_all = nb_cyc_id_all + nb_cyc(1,1);      %% nbre de cycles identiques total pour grounded
    nb_cyc_sup_all = nb_cyc_sup_all + nb_cyc(1,2);     %% nbre de cycles supplémentaires total pour grounded
    nb_cyc_miss_all = nb_cyc_miss_all + nb_cyc(1,3);    %% nbre de cycles manquants total pour grounded
    
    nb_flo_id_all = nb_flo_id_all + nb_flot(1,1);     %%nbre de flotteurs id total pour grounded
    nb_flo_sup_all = nb_flo_sup_all + nb_flot(1,2);   %%nbre de flotteurs supp total pour grounded
    nb_flo_miss_all = nb_flo_miss_all + nb_flot(1,3);    %%nbre de flotteurs manquants total pour grounded
    
    nb_flo_al_all = nb_flo_al_all + nb_flo(1);          %%% nbre de flotteurs altran
    nb_flo_test_all = nb_flo_test_all + nb_flo(2);      %%% nbre de flotteurs test
    
    flo_total_all = flo_total_all + flo_tot(1);           %%% nbre total de flotteurs à tester
    flo_without_err_al_all = flo_without_err_al_all + flo_tot(2);     %%% nbre de flotteurs sans erreurs d'apres altran
    flo_without_err_test_all = flo_without_err_test_all + flo_tot(3);   %%% nbre de flotteurs sans erreurs d'apres test
    
    
    
    % Figure 1 Nombre de flotteurs identiques/sup/missing
	%====================================================
    
    figure
    %subplot(2,1,1)
    %bar(Prcent)
	b=bar(nb_flot);

    %pour matlab2020
    % if size(nb_flot,1)==1
	% b.FaceColor = 'flat';
	 % b.CData(1,:) = [0 0.4470 0.7410];
    % b.CData(2,:) = [0.4660    0.6740    0.1880];
	 % b.CData(3,:) = [0.9290 0.6940 0.1250]
	% else
	% b(2).FaceColor=[0.4660    0.6740    0.1880];
	% end
	
	
	xticks([1:1:length(nb_flot)])
    alerte_str = cellstr(alerte);
    xticklabels(alerte_str)
    xtickangle(45)
    title({'Pourcentage de Flotteurs identiques à Altran' ...
        liste_rep_alerte{l}},'interp','none');
    titre1 = ['Nbre de Flotteurs identiques(b), supplementaires(g), manquants(y)']
    titre2 = char(liste_rep_alerte(l));
    title({titre1  titre2});
    % ylim([0 110])
    grid on
    %     Namefig = ['/home/gherbert/Images/compa_test_altran/compa_par_type/flot_id_supp_miss_' liste_rep_alerte{l} 'pourverif2.png'];
	if exist([param.DIR_HOME '/figures/' liste_rep_alerte{l}],'dir')==0
	mkdir([param.DIR_HOME '/figures/' liste_rep_alerte{l}])
	end
	
    Namefig = [param.DIR_HOME '/figures/' liste_rep_alerte{l} '/flot_id_supp_miss_' liste_rep_alerte{l} 'pourverif2.png']; % cc modif chemin sauvegarde fgiure
    saveas(gcf,Namefig)
    
    
    % Figure2 Nombre de cycles identiques/sup/missing
	%====================================================
    figure
    alerte_str = cellstr(alerte);
    b=bar(nb_cyc);
	% pour matlab2020
	% if size(nb_flot,1)==1
	% b.FaceColor = 'flat';
	 % b.CData(1,:) = [0 0.4470 0.7410];
    % b.CData(2,:) = [0.4660    0.6740    0.1880];
	 % b.CData(3,:) = [0.9290 0.6940 0.1250]
	% else
	% b(2).FaceColor=[0.4660    0.6740    0.1880];
	% end
	
    xticks([1:1:length(nb_cyc)])
	xticklabels(alerte_str)
    xtickangle(45)
    titre1 = ['Nbre de cycles identiques(b), supplementaires(g), manquants(y)']
    %titre2 = ['pour un flotteur present dans les deux cas'];
    titre3 = char(liste_rep_alerte(l));
    %title({titre1  titre2  titre3},'interp','none');
	title({titre1    titre3},'interp','none');
    %title({'Nbre de cycles identiques, supplémentaires, manquants' ...
    %    'pour un flotteur présent dans les deux cas'});
    grid on
    %     Namefig = ['/home/gherbert/Images/compa_test_altran/compa_par_type/Cyc_id_supp_miss_' liste_rep_alerte{l} 'pourverif2.png'];
    Namefig = [param.DIR_HOME '/figures/' liste_rep_alerte{l} '/Cyc_id_supp_miss_' liste_rep_alerte{l} 'pourverif2.png']; % cc modif chemin sauvegarde fgiure
    saveas(gcf,Namefig)
    %keyboard
    
	% Figure3 Nombre flotteur ALTRAN/TRAJ_QC pour chaque alerte
	%====================================================
    figure
    positionVector1 = [0.1, 0.25, 0.4, 0.5];
    subplot('Position',positionVector1)
    %subplot(2,1,1)
    bar(nb_flo)
	xticks([1:1:length(nb_cyc)])
    xticklabels(alerte_str)
    xtickangle(45)
    grid on
    title({'Nb de flotteurs par alerte'...
        '(r:altran b:test)'})
    positionVector2 = [0.6, 0.25, 0.3, 0.3];
    subplot('Position',positionVector2)
    %subplot(2,1,2)
    set(gcf, 'colormap', [1 0 0 ; 0 0 1])
    bar(flo_tot)
    xticks([1:1:3])
    xticklabels({'nbr total'; 'sans alerte - Altran'; 'sans alerte - test'})
    xtickangle(45)
    title({'Nb total de flotteurs' ...
        'et flotteurs sans alertes'});
    grid on
    %     Namefig = ['/home/gherbert/Images/compa_test_altran/compa_par_type/Nb_flot_' liste_rep_alerte{l} 'pourverif2.png'];
    Namefig = [param.DIR_HOME '/figures/' liste_rep_alerte{l} '/Nb_flot_' liste_rep_alerte{l} 'pourverif2.png']; % cc modif chemin sauvegarde fgiure
    saveas(gcf,Namefig)
    
    
    %%Sauvegarde des alertes par type de flotteurs pour faire des stats
    %%globales
    file_save = ['alertes' liste_rep_alerte{l}];
    save(file_save,'nb_cyc','nb_flo','flo_tot');
    
    
    %%%Récupération des n° de flotteurs et cycles manquant pour chaque type
    %%%d'alerte. Ecriture dans un fichier txt.
    
    fch_1 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/missing_alerte_grounded.txt'];
    fch_2 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/missing_alerte_pressure.txt'];
    fch_3 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/missing_alerte_locpos.txt'];
    fch_4 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/missing_alerte_locdate.txt'];
    fch_5 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/missing_alerte_cycle.txt'];
    fch_6 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/missing_alerte_locterre.txt'];
    fch_7 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/missing_alerte_launchdate.txt'];
    fch_8 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/missing_alerte_metadur.txt'];
    
    
    if exist(fch_1,'file')==2; delete fch_1; end;
    if exist(fch_2,'file')==2; delete fch_2; end;
    if exist(fch_3,'file')==2; delete fch_3; end;
    if exist(fch_4,'file')==2; delete fch_4; end;
    if exist(fch_5,'file')==2; delete fch_5; end;
    if exist(fch_6,'file')==2; delete fch_6; end;
    if exist(fch_7,'file')==2; delete fch_7; end;
    if exist(fch_8,'file')==2; delete fch_8; end;
    
    fid_1 = fopen(fch_1,'w');
    fid_2 = fopen(fch_2,'w');
    fid_3 = fopen(fch_3,'w');
    fid_4 = fopen(fch_4,'w');
    fid_5 = fopen(fch_5,'w');
    fid_6 = fopen(fch_6,'w');
    fid_7 = fopen(fch_7,'w');
    fid_8 = fopen(fch_8,'w');
    
    n_err=0;
    
    for i = 1:size(float_al,1)
        flo_al = double(float_alok{i}(:));
        tp = flo_al;
        tps = cell2mat(flo_miss(1,i));
        
        %keyboard
        if(isempty(find(~isnan(tps)))==0)
            for j = 1:length(tps)
                if(isnan(tps(j))==0)
                    idfloat_miss = find(tp==tps(j));
                    n_err = n_err+1;
                    % keyboard
                    if(i==1)
                        %fprintf(fid_1,'%s\n',[num2str(tps(j)) ',[' num2str(cell2mat(cycalok{i,idfloat_miss})) ']']);
						fprintf(fid_1,'%s\n',[num2str(tps(j)) ',[' num2str((cycalok{i,idfloat_miss})) ']']);
                    elseif(i==2)%
                        %fprintf(fid_2,'%s\n',[num2str(tps(j)) ',[' num2str(cell2mat(cycalok{i,idfloat_miss})) ']']);
                        fprintf(fid_2,'%s\n',[num2str(tps(j)) ',[' num2str((cycalok{i,idfloat_miss})) ']']);
                    elseif(i==3)
                        %fprintf(fid_3,'%s\n',[num2str(tps(j)) ',[' num2str(cell2mat(cycalok{i,idfloat_miss})) ']']);
                        fprintf(fid_3,'%s\n',[num2str(tps(j)) ',[' num2str((cycalok{i,idfloat_miss})) ']']);
                    elseif(i==4)
                        %fprintf(fid_4,'%s\n',[num2str(tps(j)) ',[' num2str(cell2mat(cycalok{i,idfloat_miss})) ']']);
                        fprintf(fid_4,'%s\n',[num2str(tps(j)) ',[' num2str((cycalok{i,idfloat_miss})) ']']);
                    elseif(i==5)
                        %fprintf(fid_5,'%s\n',[num2str(tps(j)) ',[' num2str(cell2mat(cycalok{i,idfloat_miss})) ']']);
                        fprintf(fid_5,'%s\n',[num2str(tps(j)) ',[' num2str((cycalok{i,idfloat_miss})) ']']);
                    elseif(i==6)
                        %fprintf(fid_6,'%s\n',[num2str(tps(j)) ',[' num2str(cell2mat(cycalok{i,idfloat_miss})) ']']);
                        fprintf(fid_6,'%s\n',[num2str(tps(j)) ',[' num2str((cycalok{i,idfloat_miss})) ']']);
                    elseif(i==7)
                        %fprintf(fid_7,'%s\n',[num2str(tps(j)) ',[' num2str(cell2mat(cyc_alok{i,idfloat_miss})) ']']);
                        fprintf(fid_7,'%s\n',[num2str(tps(j)) ',[' num2str((cyc_alok{i,idfloat_miss})) ']']);
                    elseif(i==8)
                        %fprintf(fid_8,'%s\n',[num2str(tps(j)) ',[' num2str(cell2mat(cyc_alok{i,idfloat_miss})) ']']);
                        fprintf(fid_8,'%s\n',[num2str(tps(j)) ',[' num2str((cyc_alok{i,idfloat_miss})) ']']);
                    end
                end
            end
        end
    end
    
    
    fclose(fid_1);fclose(fid_2); fclose(fid_3);fclose(fid_4);
    fclose(fid_5);fclose(fid_6);fclose(fid_7);
    
    
    %%%Ecriture des cycles manquants dans test et flotteurs correspondants
    
    %keyboard
	fch_11 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/cycmissing_alerte_grounded.txt'];
    fch_22 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/cycmissing_alerte_pressure.txt'];
    fch_33 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/cycmissing_alerte_locpos.txt'];
    fch_44 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/cycmissing_alerte_locdate.txt'];
    fch_55 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/cycmissing_alerte_cycle.txt'];
    fch_66 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/cycmissing_alerte_locterre.txt'];
    fch_77 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/cycmissing_alerte_launchdate.txt'];
    fch_88 = [param.DIR_HOME '/alerts/' liste_rep_alerte{l} '/cycmissing_alerte_metadur.txt'];
    
    % fch_11 = ['./alerts/alerte_grounded_cycmissing_2_' liste_rep_alerte{l} '.txt'];
    % fch_22 = ['./alerts/alerte_press_cycmissing_2_' liste_rep_alerte{l} '.txt'];
    % fch_33 = ['./alerts/alerte_locpos_cycmissing_2_' liste_rep_alerte{l} '.txt'];
    % fch_44 = ['./alerts/alerte_locdate_cycmissing_2_' liste_rep_alerte{l} '.txt'];
    % fch_55 = ['./alerts/alerte_cycle_cycmissing_2_' liste_rep_alerte{l} '.txt'];
    % fch_66 = ['./alerts/alerte_locterre_cycmissing_2_' liste_rep_alerte{l} '.txt'];
    % fch_77 = ['./alerts/alerte_launchdate_mcycissing_2_' liste_rep_alerte{l} '.txt'];
    % fch_88 = ['./alerts/alerte_metadur_mcycissing_2_' liste_rep_alerte{l} '.txt'];
    
    
    if exist(fch_11,'file')==2; delete fch_11; end;
    if exist(fch_22,'file')==2; delete fch_22; end;
    if exist(fch_33,'file')==2; delete fch_33; end;
    if exist(fch_44,'file')==2; delete fch_44; end;
    if exist(fch_55,'file')==2; delete fch_55; end;
    if exist(fch_66,'file')==2; delete fch_66; end;
    if exist(fch_77,'file')==2; delete fch_77; end;
    if exist(fch_88,'file')==2; delete fch_88; end;
    
    fid_11 = fopen(fch_11,'w');
    fid_22 = fopen(fch_22,'w');
    fid_33 = fopen(fch_33,'w');
    fid_44 = fopen(fch_44,'w');
    fid_55 = fopen(fch_55,'w');
    fid_66 = fopen(fch_66,'w');
    fid_77 = fopen(fch_77,'w');
    fid_88 = fopen(fch_88,'w');
    
    
    
    for i = 1:size(float_al,1)
        
        for j = 1:size(cyc_test_miss,2)
            % cyc_miss = cell2mat(cyc_test_miss{i,j});
            if(i==1 & size(cyc_test_miss,1)>=1)
                fprintf(fid_11,'%s\n',num2str(cyc_test_miss{i,j}));
            elseif(i==2 & size(cyc_test_miss,1)>=2)
                fprintf(fid_22,'%s\n',[num2str(cyc_test_miss{i,j})]);
            elseif(i==3 & size(cyc_test_miss,1)>=3)
                fprintf(fid_33,'%s\n',[num2str(cyc_test_miss{i,j})]);
            elseif(i==4 & size(cyc_test_miss,1)>=4)
                fprintf(fid_44,'%s\n',[num2str(cyc_test_miss{i,j})]);
            elseif(i==5 & size(cyc_test_miss,1)>=5)
                fprintf(fid_55,'%s\n',[num2str(cyc_test_miss{i,j})]);
            elseif(i==6 & size(cyc_test_miss,1)>=6)
                fprintf(fid_66,'%s\n',[num2str(cyc_test_miss{i,j})]);
            elseif(i==7 & size(cyc_test_miss,1)>=7)
                fprintf(fid_77,'%s\n',[num2str(cyc_test_miss{i,j})]);
            elseif(i==8 & size(cyc_test_miss,1)>=8)
                fprintf(fid_88,'%s\n',[num2str(cyc_test_miss{i,j})]);
                
            end
        end
    end
    fclose(fid_11);fclose(fid_22); fclose(fid_33);fclose(fid_44);
    fclose(fid_55);fclose(fid_66);fclose(fid_77);fclose(fid_88);
end   %% fin de la boucle sur les listes



%%%%   ------------------------------------------------------------------
%%%% Figure alertes pour l'ensemble des types de floteurs. (attention pour la première erreur de la liste seulement)

nb_cyc_all = [nb_cyc_id_all nb_cyc_sup_all nb_cyc_miss_all];
nb_flo_all = [nb_flo_id_all*100/flo_total_all nb_flo_sup_all*100/flo_total_all nb_flo_miss_all*100/flo_total_all];
nb_flo = [nb_flo_al_all nb_flo_test_all];
flo_all = [flo_total_all flo_without_err_al_all flo_without_err_test_all];


figure
bar(nb_cyc_all)
xticks([1:1:length(nb_cyc_all)])
xticklabels({'Id' 'Supp' 'Miss'})
xtickangle(45)
titre1 = ['Nbre de cycles identiques(b), supplementaires(g), manquants(y)']
titre2 = ['pour un flotteur present dans les deux cas'];
titre3 = char(liste_rep_alerte(l));
title({titre1  titre2  titre3});
%title({'Nbre de cycles identiques, supplémentaires, manquants' ...
%    'pour un flotteur présent dans les deux cas'});
grid on
if exist([param.DIR_HOME 'figures/compa_test_altran'],'dir')==0
mkdir([param.DIR_HOME 'figures/compa_test_altran'])
end
Namefig = [param.DIR_HOME 'figures/compa_test_altran/Cyc_id_supp_miss_all2.png'];
saveas(gcf,Namefig)

figure
bar(nb_flo_all)
xticks([1:1:length(nb_flo_all)])
xticklabels({'Id' 'Supp' 'Miss'})
xtickangle(45)
titre1 = ['% de flotteurs identiques, supplementaires, manquants']
title({titre1});
grid on
%Namefig = ['/home/gherbert/Images/compa_test_altran/compa_par_type/flo_id_supp_miss_all2.png'];
Namefig = [param.DIR_HOME 'figures/compa_test_altran/flo_id_supp_miss_all2.png'];
saveas(gcf,Namefig)


figure
positionVector1 = [0.1, 0.25, 0.4, 0.5];
subplot('Position',positionVector1)
%subplot(2,1,1)
bar(nb_flo)
xticks([1:1:length(nb_flo)])
xticklabels(alerte_str)
xtickangle(45)
grid on
title({'Nb de flotteurs par alerte'...
    '(r:altran b:test)'})
positionVector2 = [0.6, 0.25, 0.3, 0.3];
subplot('Position',positionVector2)
%subplot(2,1,2)
set(gcf, 'colormap', [1 0 0 ; 0 0 1])
bar(flo_all)
xticks([1:1:3])
xticklabels({'nbr total'; 'sans alerte - Altran'; 'sans alerte - test'})
xtickangle(45)
title({'Nb total de flotteurs' ...
    'et flotteurs sans alertes'});
grid on
%Namefig = ['/home/gherbert/Images/compa_test_altran/compa_par_type/Nb_flot_all2.png'];
Namefig = [param.DIR_HOME 'figures/compa_test_altran/Nb_flot_all2.png'];
saveas(gcf,Namefig)







%%% Récupération des flotteurs/cycles supplémentaires par rapport à Altran



%keyboard
%fch_11 = '/home/gherbert/andro/andro/soft/solene/alerte_cycle_cycsup.txt';
%fch_22 = '/home/gherbert/andro/andro/soft/solene/alerte_grounded_cycsup_2.txt';
%fch_33 = '/home/gherbert/andro/andro/soft/solene/alerte_launchdate_cycsup_2.txt';
%fch_44 = '/home/gherbert/andro/andro/soft/solene/alerte_locpos_cycsup_2.txt';
%fch_55 = '/home/gherbert/andro/andro/soft/solene/alerte_locterre_cycsup_2.txt';
%fch_66 = '/home/gherbert/andro/andro/soft/solene/alerte_metadur_cycsup_2.txt';
%fch_77 = '/home/gherbert/andro/andro/soft/solene/alerte_press_mcycsup_2.txt';

%fid_11 = fopen(fch_11,'w');
%fid_22 = fopen(fch_22,'w');
%fid_33 = fopen(fch_33,'w');
%fid_44 = fopen(fch_44,'w');
%fid_55 = fopen(fch_55,'w');
%fid_66 = fopen(fch_66,'w');
%fid_77 = fopen(fch_77,'w');


%for i = 1:size(float,1)

%    for j = 1:size(recup_cycsup,2)
%        if(i==1 & isempty(recup_cycsup{i,j})==0)
%           fprintf(fid_11,'%s\n',num2str(recup_cycsup{i,j}));
%        elseif(i==2 & isempty(recup_cycsup{i,j})==0)
%          fprintf(fid_22,'%s\n',[num2str(recup_cycsup{i,j})]);
%       elseif(i==3 & isempty(recup_cycsup{i,j})==0)
%          fprintf(fid_33,'%s\n',[num2str(recup_cycsup{i,j})]);
%       elseif(i==4 & isempty(recup_cycsup{i,j})==0)
%          fprintf(fid_44,'%s\n',[num2str(recup_cycsup{i,j})]);
%       elseif(i==5 & isempty(recup_cycsup{i,j})==0)
%          fprintf(fid_55,'%s\n',[num2str(recup_cycsup{i,j})]);
%       elseif(i==6 & isempty(recup_cycsup{i,j})==0)
%          fprintf(fid_66,'%s\n',[num2str(recup_cycsup{i,j})]);
%       elseif(i==7 & isempty(recup_cycsup{i,j})==0)
%          fprintf(fid_77,'%s\n',[num2str(recup_cycsup{i,j})]);
%       end
%   end
%end
%fclose(fid_11);fclose(fid_22); fclose(fid_33);fclose(fid_44);
%fclose(fid_55);fclose(fid_66);fclose(fid_77);



