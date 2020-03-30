%remove all the first line of the code
%load('DNn_mt.mat') % load .mat file
%metabolitesList=mt([1:29649]);
%compoundList=mt([1 2 8 10 25 26 29 31 112 458 559 1232 18610:18613]);
function Metabolites=nc_SortDataCompounds(metabolitesList)
tic;
Metabolites=struct([]);
numCompounds=length(metabolitesList);
handleWaitbar=waitbar(0,'Please wait...');
for i=1:numCompounds
% for i=5:numCompounds
    waitbar(i/numCompounds,handleWaitbar,['Sorting compounds: ' num2str(i) ' of ' num2str(numCompounds) 'for' num2str(toc)]);
    raw_data=metabolitesList(i).txt; % get specific entry based on ind (index)
    new_lines=regexp(raw_data,'[\n]'); % find all new lines
    new_lines_with_space=regexp(raw_data,'[\n] '); % find all new lines where the next line starts with a space
    new_lines(ismember(new_lines,new_lines_with_space))=[]; % remove the ones that start with a space from the list in f
    new_lines=new_lines+1;
    new_lines=[1 new_lines];  %#ok % add  first text character
    Metabolites(i).KEGG_ID=metabolitesList(i).KEGG_ID;
    %find entry
    index_position=regexp(raw_data,'ENTRY');
    k=find(new_lines==index_position);
    Metabolites(i).ENTRY=regexprep(raw_data(new_lines(k):(new_lines(k)+18)),'ENTRY','');
    Metabolites(i).ENTRY=regexprep(Metabolites(i).ENTRY,'\s+',' ');
    Metabolites(i).ENTRY=strtrim(Metabolites(i).ENTRY);
    %names
    index_position=regexp(raw_data,'NAME');
    if ~isempty(index_position)
        k=find(new_lines==index_position(1)); 
        if isempty(k)%==1
            k=find(new_lines<index_position);
            k=k(end); 
        else
            
        end
        Metabolites(i).NAME=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'NAME','');
        Metabolites(i).NAME=regexprep(Metabolites(i).NAME,'\s+',' ');
        Metabolites(i).NAME=split(Metabolites(i).NAME,{'; '});
        Metabolites(i).NAME=strtrim(Metabolites(i).NAME);
        Length_of_name=numel(Metabolites(i).NAME);
        str_to_find='L-';
        for n=1:Length_of_name
            if length(Metabolites(i).NAME{n,1})>=length(str_to_find) && strcmp(Metabolites(i).NAME{n,1}(1:length(str_to_find)),str_to_find)==1
                Metabolites(i).NAME{end+1,1}=Metabolites(i).NAME{n,1}(numel(str_to_find)+1:numel(Metabolites(i).NAME{n,1}));
            end
        end
        str_to_find='D-';
        for n=1:Length_of_name
            if length(Metabolites(i).NAME{n,1})>=length(str_to_find) && strcmp(Metabolites(i).NAME{n,1}(1:length(str_to_find)),str_to_find)==1
                Metabolites(i).NAME{end+1,1}=Metabolites(i).NAME{n,1}(numel(str_to_find)+1:numel(Metabolites(i).NAME{n,1}));
            end
        end
        str_to_find='(R)-';
        for n=1:Length_of_name
            if length(Metabolites(i).NAME{n,1})>=length(str_to_find) && strcmp(Metabolites(i).NAME{n,1}(1:length(str_to_find)),str_to_find)==1
               Metabolites(i).NAME{end+1,1}=Metabolites(i).NAME{n,1}(numel(str_to_find)+1:numel(Metabolites(i).NAME{n,1}));
            end
        end
        str_to_find='(S)-';
        for n=1:Length_of_name
            if length(Metabolites(i).NAME{n,1})>=length(str_to_find) && strcmp(Metabolites(i).NAME{n,1}(1:length(str_to_find)),str_to_find)==1
               Metabolites(i).NAME{end+1,1}=Metabolites(i).NAME{n,1}(numel(str_to_find)+1:numel(Metabolites(i).NAME{n,1}));
            end
        end
        str_to_find='(Z)-';
        for n=1:Length_of_name
             if length(Metabolites(i).NAME{n,1})>=length(str_to_find) && strcmp(Metabolites(i).NAME{n,1}(1:length(str_to_find)),str_to_find)==1
                 Metabolites(i).NAME{end+1,1}=Metabolites(i).NAME{n,1}(numel(str_to_find)+1:numel(Metabolites(i).NAME{n,1}));
             end
         end
          str_to_find='(E)-';
         for n=1:Length_of_name
             if length(Metabolites(i).NAME{n,1})>=length(str_to_find) && strcmp(Metabolites(i).NAME{n,1}(1:length(str_to_find)),str_to_find)==1
                 Metabolites(i).NAME{end+1,1}=Metabolites(i).NAME{n,1}(numel(str_to_find)+1:numel(Metabolites(i).NAME{n,1}));
             end
         end
         str_to_find='R-';
        for n=1:Length_of_name
            if length(Metabolites(i).NAME{n,1})>=length(str_to_find) && strcmp(Metabolites(i).NAME{n,1}(1:length(str_to_find)),str_to_find)==1
                Metabolites(i).NAME{end+1,1}=Metabolites(i).NAME{n,1}(numel(str_to_find)+1:numel(Metabolites(i).NAME{n,1}));
            end
        end
         str_to_find='S-';
        for n=1:Length_of_name
            if length(Metabolites(i).NAME{n,1})>=length(str_to_find) && strcmp(Metabolites(i).NAME{n,1}(1:length(str_to_find)),str_to_find)==1
                Metabolites(i).NAME{end+1,1}=Metabolites(i).NAME{n,1}(numel(str_to_find)+1:numel(Metabolites(i).NAME{n,1}));
            end
        end
    end
        Metabolites(i).NAME=unique(Metabolites(i).NAME);
        Name_for_glycan=Metabolites(i).NAME;
 
        %MASS
        index_position=regexp(raw_data,'MASS');
        if isempty(index_position)==0
            k=find(new_lines==index_position); 
            Metabolites(i).MASS=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'MASS','');
            Metabolites(i).MASS=regexprep(Metabolites(i).MASS,'\s+',' ');
            Metabolites(i).MASS=strtrim(Metabolites(i).MASS);
            %if ~isempty(Metabolites(i).MASS)
             %   Metabolites(i).MASS=strcat(Metabolites(i).MASS,' ');
            position_of_space=regexp(Metabolites(i).MASS,' ');
            if ~isempty(position_of_space)
                Metabolites(i).MASS=Metabolites(i).MASS(1:position_of_space(1));
            end
        end
        if isempty(Metabolites(i).MASS)
            index_position=regexp(raw_data,'MOL_WEIGHT');
            if isempty(index_position)==0
                k=find(new_lines==index_position); 
                Metabolites(i).MASS=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'MOL_WEIGHT','');
                Metabolites(i).MASS=regexprep(Metabolites(i).MASS,'\s+',' ');
                Metabolites(i).MASS=strtrim(Metabolites(i).MASS);
            end
        end
        %reaction
        index_position=regexp(raw_data,'REACTION');
        if isempty(index_position)==0
            k=find(new_lines==index_position); 
            Metabolites(i).REACTION=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'REACTION','');
            Metabolites(i).REACTION=regexprep(Metabolites(i).REACTION,'\s+',' ');
            Metabolites(i).REACTION=strtrim(Metabolites(i).REACTION);
            Metabolites(i).REACTION=split(Metabolites(i).REACTION,{' '});
            Metabolites(i).REACTION=strtrim(Metabolites(i).REACTION);
        end
        %PATHWAY
        index_position=regexp(raw_data,'PATHWAY');
        if isempty(index_position)==0
            index_position=index_position(1);
            k=find(new_lines==index_position(1)); 
            Metabolites(i).PATHWAY=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'PATHWAY','');
            Metabolites(i).PATHWAY=regexprep(Metabolites(i).PATHWAY,'\s+',' ');
            Metabolites(i).PATHWAY=split(Metabolites(i).PATHWAY,{'map'});
            Metabolites(i).PATHWAY(1,:)=[];
            for k=1:length(Metabolites(i).PATHWAY)
                Metabolites(i).PATHWAY{k}=Metabolites(i).PATHWAY{k}(1:5);
            end
            %m(i).PATHWAY{i}= m(i).PATHWAY{i}(1:5)
            %m(i).PATHWAY=insertBefore(m(i).PATHWAY,1,'map');
            Metabolites(i).PATHWAY=strtrim(Metabolites(i).PATHWAY);
        else
            Metabolites(i).PATHWAY=[];
        end
        %ENZYME
        index_position=regexp(raw_data,'ENZYME');
        if isempty(index_position)==0
            k=find(new_lines==index_position(1)); 
            Metabolites(i).ENZYME=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'ENZYME','');
            Metabolites(i).ENZYME=regexprep(Metabolites(i).ENZYME,'\s+',' ');
            Metabolites(i).ENZYME=strtrim(Metabolites(i).ENZYME);
            Metabolites(i).ENZYME=split(Metabolites(i).ENZYME,' ');
            Metabolites(i).ENZYME=strtrim(Metabolites(i).ENZYME);
        end
        if isempty(Metabolites(i).NAME)
            index_position=regexp(raw_data,'COMPOSITION');
            if ~isempty(index_position)
                k=find(new_lines==index_position, 1); 
                if isempty(k)%==1 Since Glycans here are not in a format that can be analysed by this code, the code below enters all the data for the glycans
                else
                    k=find(new_lines==index_position); 
                    Metabolites(i).NAME=regexprep(raw_data(new_lines(k)+11:(new_lines(k+1)-1)),'NAME','');
                    Metabolites(i).NAME=regexprep(Metabolites(i).NAME,'\s+',' ');
                    Metabolites(i).COMPOSITION=strtrim(split(Metabolites(i).NAME,{'; '}));
                    Metabolites(i).NAME=split(Metabolites(i).NAME,{'; '});
                    Metabolites(i).NAME=strtrim(Metabolites(i).NAME); 
                end
            end
        elseif ismember(Metabolites(i).ENTRY(1),'G')
                index_position=regexp(raw_data,'COMPOSITION');
                if ~isempty(index_position)
                    k=find(new_lines==index_position); 
                    Metabolites(i).NAME=regexprep(raw_data(new_lines(k)+11:(new_lines(k+1)-1)),'NAME','');
                    Metabolites(i).NAME=regexprep(Metabolites(i).NAME,'\s+',' ');
                    Metabolites(i).COMPOSITION=strtrim(split(Metabolites(i).NAME,{'; '}));
                    Metabolites(i).NAME=split(Metabolites(i).NAME,{'; '});
                    Metabolites(i).NAME=strtrim(Metabolites(i).NAME);
                    Metabolites(i).NAME=strcat(Name_for_glycan,'@',Metabolites(i).NAME);
                    Metabolites(i).NAME=split(Metabolites(i).NAME,'@');
                    if size(Metabolites(i).NAME,2)>1
                        Metabolites(i).NAME=cat(1,Metabolites(i).NAME(:,1),Metabolites(i).NAME(:,2));
                    end
                end
        end
        index_position=regexp(raw_data,'REMARK');
        if isempty(index_position)==0
            k=find(new_lines==index_position); 
            Metabolites(i).SAME_AS=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'REMARK','');
            Metabolites(i).SAME_AS=regexprep(Metabolites(i).SAME_AS,'\s+',' ');
            Metabolites(i).SAME_AS=strtrim(Metabolites(i).SAME_AS);
            position_of_colon=regexp(Metabolites(i).SAME_AS,':');
            Metabolites(i).SAME_AS=Metabolites(i).SAME_AS(position_of_colon+2:end);
            if strcmp(Metabolites(i).SAME_AS(1),'C')
                Metabolites(i).SAME_AS=Metabolites(i).SAME_AS;
            elseif strcmp(Metabolites(i).SAME_AS(1),'G')
                Metabolites(i).SAME_AS=Metabolites(i).SAME_AS;
            else
                Metabolites(i).SAME_AS=[];
            end
        end
end
for i=numCompounds:-1:1
    if isempty(Metabolites(i).NAME)
        Metabolites(i)=[];
    end
end
close(handleWaitbar)
disp(toc)
end % function

%dont delete weight and exact mass(done)
%compare weight between 2 compounds before merging(done)
%absolute difference between masses <5 (for now) make varibale
%max_diff=x(done)
%should only take C0000x id if it is larger than the first C000x??
%for later
%make a list of deleted(count glycans/compounds)-next--> problem with
%List_of_index--> does not delete/remove glycan duplicates......
%make a function to match names to compound code(Metabolite_expand_name and
%output id from Metabolite_expand_id
%search pubmed(search terms)- metabolites,metabolome associated with
%colorectal cancer(search queries- choose which one to pick,literature review)

%check if reaction pair(replacing name works properly)

%find if chromosome is the same--> and if the snp is larger or smaller than
%the gene(to know if the snp is inside or outside the gene) if outside gene
%diff (done:D)

%make a workflow for all the functions-->in a way that someone can
%replicate the final database
