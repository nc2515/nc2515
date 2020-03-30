%ClassList=ct([1:3161]);
function Class=nc_SortDataClass(ClassList)
    Class=struct([]);
    numClass=length(ClassList);
    handleWaitbar=waitbar(0,'Please wait...');
    for i=1:numClass
        waitbar(i/numClass,handleWaitbar,['Processing ' num2str(i) ' of ' num2str(numClass)]);
        raw_data=ClassList(i).txt; % get specific entry based on ind (index)
        new_lines=regexp(raw_data,'[\n]'); % find all new lines
        new_lines_with_space=regexp(raw_data,'[\n] '); % find all new lines where the next line starts with a space
        new_lines(ismember(new_lines,new_lines_with_space))=[]; % remove the ones that start with a space from the list in f
        new_lines=new_lines+1;
        new_lines=[1 new_lines]; %#ok add the first text character
        Class(i).KEGG_ID=ClassList(i).KEGG_ID;
        %find entry
        index=regexp(raw_data,'ENTRY'); %find the position of ENTRY in txt
        k=find(new_lines==index); %find th position where the position of ENTRY=f(where ENTRY starts in a new line)
        Class(i).ENTRY=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'ENTRY','');%replace(txt(from postion f(k):before the start of new line)-> replace entry with empty space(remove entry from the str)
        Class(i).ENTRY=regexprep(Class(i).ENTRY,'\s+',' ');
        Class(i).ENTRY=regexprep(Class(i).ENTRY,'EC',' ');
        Class(i).ENTRY=strtrim(Class(i).ENTRY);
        index_of_space=regexp(Class(i).ENTRY,' ');
        Class(i).ENTRY=Class(i).ENTRY(1:index_of_space-1);
        %DEFINITION
        index=regexp(raw_data,'DEFINITION');%find the position of NAME in txt
        if ~isempty(index)%if the position of NAME is not empty-->
            k=find(new_lines==index(1));%K=position where NAME=new line, find the position where the position of NAME=f(where NAME starts in a new line) *this lines give me error 'Matrix dimensions must agree.'*
            if isempty(k)%==1
                k=find(new_lines<index);%K=find psotion where new line is greater than position of NAME in text 
                k=k(end); 
            end
            Class(i).DEFINITION=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'DEFINITION','');
            Class(i).DEFINITION=regexprep(Class(i).DEFINITION,'\s+',' ');
            %m(i).DEFINITION=split(m(i).DEFINITION,{'; '});
            %m(i).DEFINITION=strtrim(m(i).DEFINITION);
        end   
        %RPAIR
        index=regexp(raw_data,'RPAIR');
        if isempty(index)==0
            index=index(1);
            k=find(new_lines==index);
            Class(i).RPAIR=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'RPAIR','');
            Class(i).RPAIR=regexprep(Class(i).RPAIR,'\s+',' ');
            Class(i).RPAIR=strtrim(Class(i).RPAIR);
            Class(i).RPAIR=split(Class(i).RPAIR,{' '});
        %    m(i).ORGANISM=insertBefore(m(i).GENES,1,'map');
            %for k=1:length(m(i).RCLASS)
            %   m(i).RCLASS{k}=m(i).RCLASS{k}(3:length(m(i).RCLASS));
            %end
            %m(i).RCLASS=split(m(i).CLASS,'RC');
        end
        %REACTION
        index=regexp(raw_data,'REACTION');
        if isempty(index)==0
            index=index(1);
            k=find(new_lines==index); 
            Class(i).REACTION=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'REACTION','');
            Class(i).REACTION=regexprep(Class(i).REACTION,'\s+',' ');
            Class(i).REACTION=strtrim(Class(i).REACTION);
            Class(i).REACTION=split(Class(i).REACTION,{' '});
        end
        %ENZYME
        index=regexp(raw_data,'ENZYME');
        if isempty(index)==0
            index=index(1);
            k=find(new_lines==index);
            Class(i).ENZYME=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'ENZYME','');
            Class(i).ENZYME=regexprep(Class(i).ENZYME,'\s+',' ');
            Class(i).ENZYME=strtrim(Class(i).ENZYME);
            Class(i).ENZYME=split(Class(i).ENZYME,{' '});
        end
        %PATHWAY
        index=regexp(raw_data,'PATHWAY');
        if ~isempty(index)
            k=find(new_lines==index);
            if isempty(k)%==1
                k=find(new_lines<index);
                k=k(end);
            end
            Class(i).PATHWAY=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'PATHWAY','');
            Class(i).PATHWAY=regexprep(Class(i).PATHWAY,'\s+',' ');
            Class(i).PATHWAY=strtrim(Class(i).PATHWAY);
            Class(i).PATHWAY=split(Class(i).PATHWAY,{'rn'});
            %for x=1:length(m(i).PATHWAY);
            %     s=regexp(m(i).PATHWAY{x,1},'\s+');
            %    m(i).PATHWAY{x,1}=m(i).PATHWAY{x,1}(1:s(1)-1);
            %end
        end
        %ORTHOLOGY
        index=regexp(raw_data,'ORTHOLOGY');
        if isempty(index)==0
            k=find(new_lines==index);
            Class(i).ORTHOLOGY=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'ORTHOLOGY','');
            %m(i).GENES=regexprep(m(i).GENES,'\s+',' ');
            %m(i).GENES=split(m(i).GENES,{'?'});
            %s=regexp(m(i).GENES,'[\n] ');
             Class(i).ORTHOLOGY=split(Class(i).ORTHOLOGY,{'            '});
             Class(i).ORTHOLOGY=strtrim(Class(i).ORTHOLOGY);
        end
    end
    close(handleWaitbar)
end