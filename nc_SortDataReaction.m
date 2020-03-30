%reactionList=rt([1:11296]);
%reactionList=rt([1:100 1000:2000 11200:11296]);
function Reaction=nc_SortDataReaction(reactionList)
    Reaction=struct([]);
    number_of_Reactions=length(reactionList);
    handleWaitbar=waitbar(0,'Please wait...');
    for i=1:number_of_Reactions
        waitbar(i/number_of_Reactions,handleWaitbar,['Sorting metabolite & reaction: ' num2str(i) ' of ' num2str(number_of_Reactions)]);
        raw_data=reactionList(i).txt; % get specific entry based on ind (index)
        new_lines=regexp(raw_data,'[\n]'); % find all new lines
        new_lines_with_space=regexp(raw_data,'[\n] '); % find all new lines where the next line starts with a space
        new_lines(ismember(new_lines,new_lines_with_space))=[]; % remove the ones that start with a space from the list in f
        new_lines=new_lines+1;
        new_lines=[1 new_lines]; %#ok % add the first text character
        Reaction(i).KEGG_ID=reactionList(i).KEGG_ID;
        %find entry
        index_position=regexp(raw_data,'ENTRY'); %find the position of ENTRY in txt
        k=find(new_lines==index_position); %find th position where the position of ENTRY(where ENTRY starts in a new line)
        Reaction(i).ENTRY=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'ENTRY','');%replace(txt(from postion f(k):before the start of new line)-> replace entry with empty space(remove entry from the str)
        Reaction(i).ENTRY=regexprep(Reaction(i).ENTRY,'\s+',' ');
        Reaction(i).ENTRY=strtrim(Reaction(i).ENTRY);
        Reaction(i).ENTRY=Reaction(i).ENTRY(1:6);
        %names
        index_position=regexp(raw_data,'NAME');%find the position of NAME in txt
        if ~isempty(index_position)%if the position of NAME is not empty-->
            k=find(new_lines==index_position(1));%K=position where NAME=new line, find the position where the position of NAME=f(where NAME starts in a new line) *this lines give me error 'Matrix dimensions must agree.'*
            if isempty(k)%==1
                k=find(new_lines<index_position);%K=find psotion where new line is greater than position of NAME in text 
                k=k(end); 
            else
                
            end
            Reaction(i).NAME=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'NAME','');
            Reaction(i).NAME=regexprep(Reaction(i).NAME,'\s+',' ');
            Reaction(i).NAME=split(Reaction(i).NAME,{'; '});
            Reaction(i).NAME=strtrim(Reaction(i).NAME);
        end     
        %DEFINITION
        %if isempty(m(i).DEFINITION)
%        index_position=regexp(raw_data,'DEFINITION');
%        if ~isempty(index_position)
%            k=find(new_lines==index_position); 
%            if isempty(k)%==1
%                k=find(new_lines<index_position);
%                k=k(end); 
%            end
%            Reaction(i).DEFINITION=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'DEFINITION','');
%            Reaction(i).DEFINITION=regexprep(Reaction(i).DEFINITION,'\s+',' ');
%            %m(i).NAME=split(m(i).NAME,{'; '});
%            Reaction(i).DEFINITION=strtrim(Reaction(i).DEFINITION);
%        end
        %end
        %Genes
%        index_position=regexp(raw_data,'EQUATION');
%        if isempty(index_position)==0
%            k=find(new_lines==index_position); 
%            Reaction(i).EQUATION=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'EQUATION','');
%            Reaction(i).EQUATION=regexprep(Reaction(i).EQUATION,'\s+',' ');
            %m(i).GENES=split(m(i).REACTION,{' '});
%        end
        %RCLASS
        index_position=regexp(raw_data,'RCLASS');
        if isempty(index_position)==0
            index_position=index_position(1);
            k=find(new_lines==index_position); 
            Reaction(i).RCLASS=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'RCLASS','');
            Reaction(i).RCLASS=regexprep(Reaction(i).RCLASS,'\s+',' ');
            Reaction(i).RCLASS=strtrim(Reaction(i).RCLASS);
       %     m(i).ORGANISM=split(m(i).GENES,{'map'});
        %    m(i).ORGANISM=insertBefore(m(i).GENES,1,'map');
            find_RC=regexp(Reaction(i).RCLASS,'R');
            for k=numel(find_RC):-1:1
                remove_list=(Reaction(i).RCLASS(find_RC(k):find_RC(k)+7));
                Reaction(i).RCLASS=erase(Reaction(i).RCLASS,remove_list);
            end
            %Reaction(i).RCLASS=Reaction(i).RCLASS(find_RC+8:RCLASS_numel);
            Reaction(i).RCLASS=strtrim(Reaction(i).RCLASS);
            Reaction(i).RCLASS=split(Reaction(i).RCLASS,' ');
            %Reaction(i).RCLASS(1,:)=[];
            Reaction(i).RCLASS=unique(Reaction(i).RCLASS);
        end
        %ENZYME
        index_position=regexp(raw_data,'ENZYME');
        if isempty(index_position)==0
            k=find(new_lines==index_position(1)); 
            Reaction(i).ENZYME=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'ENZYME','');
            Reaction(i).ENZYME=regexprep(Reaction(i).ENZYME,'\s+',' ');
            Reaction(i).ENZYME=strtrim(Reaction(i).ENZYME);
            Reaction(i).ENZYME=split(Reaction(i).ENZYME,' ');
            
        end
        %REMARK
        index_position=regexp(raw_data,'REMARK');
        if isempty(index_position)==0
            k=find(new_lines==index_position); 
            Reaction(i).REMARK=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'REMARK','');
            Reaction(i).REMARK=regexprep(Reaction(i).REMARK,'\s+',' ');
            Reaction(i).REMARK=strtrim(Reaction(i).REMARK);
            %position_R=regexp(Reaction(i).REMARK,'R');
            %length_REMARK=length(Reaction(i).REMARK);
            %Reaction(i).REMARK=Reaction(i).REMARK(position_R+1:length_REMARK);
        %    m(i).ENZYME=split(m(i).ENZYME,' ')
        end
        %PATHWAY
        index_position=regexp(raw_data,'PATHWAY');
        if isempty(index_position)==0
            k=find(new_lines==index_position); 
            Reaction(i).PATHWAY=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'PATHWAY','');
            Reaction(i).PATHWAY=regexprep(Reaction(i).PATHWAY,'\s+',' ');
            Reaction(i).PATHWAY=strtrim(Reaction(i).PATHWAY);
            length_PATHWAY=length(Reaction(i).PATHWAY);
            
            rn_position=regexp(Reaction(i).PATHWAY,' rn');
            Indexing_Pathway=Reaction(i).PATHWAY;
            pat=[];
            for h=1:length(rn_position)
                if isempty(Reaction(i).PATHWAY)
                   Reaction(i).PATHWAY={h};
                   Reaction(i).PATHWAY=Indexing_Pathway(rn_position(h)+3:rn_position(h)+7);
                else
                    %Reaction(i).PATHWAY(h)=cellstr(Indexing_Pathway(rn_position(h)+3:rn_position(h)+7));
                    pat=cat(1,pat,cellstr(Indexing_Pathway((rn_position(h)+3):(rn_position(h)+7))));
                end
            end
            Reaction(i).PATHWAY=pat;
                
            
            %Reaction(i).PATHWAYNAME=Reaction(i).PATHWAY(8:length_PATHWAY);
           % Reaction(i).PATHWAY=Reaction(i).PATHWAY(3:7);
            %Reaction(i).PATHWAYNAME=strtrim(Reaction(i).PATHWAYNAME);

        end
        %Orthology
        index_position=regexp(raw_data,'ORTHOLOGY');
        if isempty(index_position)==0
            k=find(new_lines==index_position); 
            Reaction(i).ORTHOLOGY=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'ORTHOLOGY','');
            Reaction(i).ORTHOLOGY=regexprep(Reaction(i).ORTHOLOGY,'\s+',' ');
        %    m(i).ENZYME=split(m(i).ENZYME,' ')
            Reaction(i).ORTHOLOGY=strtrim(Reaction(i).ORTHOLOGY);
            Reaction(i).ORTHOLOGY=Reaction(i).ORTHOLOGY(1:7);
            Reaction(i).ORTHOLOGY=strtrim(Reaction(i).ORTHOLOGY);
        end
    end
    close(handleWaitbar)
end