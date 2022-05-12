function Pr = format_flags_num2char(Co);
% -========================================================
%   USAGE : Pr = format_flags_num2char(Co);
%   PURPOSE : change flag char numerical vectors to string
%             ex [1 1 1 1 1 4 4 1 1 1 1 999] _> '11111441111 ' 
% -----------------------------------
%   INPUT :
%     IN1   (class)  -comments-
%             additional description
%     IN2   (class)  -comments-
%
%   OPTIONNAL INPUT :
%    OPTION1  (class)  -comments-
% -----------------------------------
%   OUTPUT :
%     OUT1   (class)  -comments-
%             additional description
%     OUT2   (class)  -comments-
%             additional description
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%            : modified (yyyy) byxxx
%   CALLED SUBROUTINES: none
% ========================================================


fillval='FillValue_';

champs = fieldnames(Co);    %champs={'psal','psalqc','psalad',....}
Nbfields = length(champs);


if isfield(Co,'fillisnan')
   original_fill=Co.fillisnan;
end
Co=replace_fill_bynan(Co); % on remplace les fillValues par des NaNs. 
Pr=Co;                           % pour les qc (si tableau numerique), cela veut dire qu'on remplace les 999 par des nan

for k=1:Nbfields            % boucle sur toutes les variables
    oneChamp=champs{k};
    
    if isempty(findstr(oneChamp,'_qc'))==0            % test si il y a 'qc' dans le nom du champ
        
        if isempty(Pr.(oneChamp).data)==0            % test si le champ est rempli
            
            if isfield (Pr.(oneChamp), 'ischar2num')&&Pr.(oneChamp).ischar2num==1;  % variable logique qui identifie
            % si le tableau de flag a ete transforme char-> num (=1) ou non (=0)
            % Le tableau a de flag charactere a prealablement ete transforme en vecteur numerique
            
                if length(size(Pr.(oneChamp).data))>2       % test si c'est un tableau de 2 dim
                   % warning('Does not accept array of dim >2')
                else
                %keyboard
                
                    Pr.(oneChamp)=rmfield(Pr.(oneChamp),'data');
                    Pr.(oneChamp).type=2; % single precision
                    Pr.(oneChamp).(fillval)=blanks(1);
                    Pr.(oneChamp).data=repmat(blanks(1),[(size(Co.(oneChamp).data))]); % fillval 
					sizePr=size(Pr.(oneChamp).data);
                    size1=sizePr(1);
                    size2=sizePr(2);
                    if size1<=size2  % cc 15/02/2021 accelère la fonction si size1 tres grand
						for i=1:size1                         % boucle sur chaque profil
							a = Co.(oneChamp).data(i,:);  % tableau temporaire
							b=blanks(length(a));
							ia=~isnan(a);
							c=num2str(a(ia));
							ic=~isspace(c);
							b(ia)=c(ic);
							Pr.(oneChamp).data(i,:)=b; % fillval
						end
					
					else
					
                    for i=1:size2                         % boucle sur chaque niveaux
							a = Co.(oneChamp).data(:,i);  % tableau temporaire
							b=blanks(length(a));
							ia=~isnan(a);
							c=num2str(a(ia));
							ic=~isspace(c);
							b(ia)=c(ic);
							Pr.(oneChamp).data(:,i)=b'; % fillval
						end
					end
                  
                    Pr.(oneChamp).ischar2num=0;
                    Pr.(oneChamp) = orderfields(Pr.(oneChamp),Co.(oneChamp)); % ordonne les champs comme dans la structure initiale
                end
            end
               
        end
    end
end
if original_fill==0
Pr=replace_nan_byfill(Pr);
end