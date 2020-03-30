%enzymeList=et([1:7637]);
%enzymeList=et([1:10 100:110 1000:1010]);
%I am up to clearing up NAMES GENES and ORGANISMS, I should be adding
%reaction class into the group too and then, clean up and move on onto
%DNn_rt
function Enzyme=nc_SortDataEnzyme(enzymeList)
    Enzyme=struct([]);
    numEnzymes=length(enzymeList);
    for_gene_match=["(",")","HSA:"];
    handleWaitbar=waitbar(0,'Please wait...');
    for i=1:numEnzymes
        waitbar(i/numEnzymes,handleWaitbar,['Processing ' num2str(i) ' of ' num2str(numEnzymes) 'for' num2str(toc)]);
        raw_data=enzymeList(i).txt; % get specific entry based on ind (index)
        new_lines=regexp(raw_data,'[\n]'); % find all new lines
        new_lines_with_space=regexp(raw_data,'[\n] '); % find all new lines where the next line starts with a space
        new_lines(ismember(new_lines,new_lines_with_space))=[]; % remove the ones that start with a space from the list in f
        new_lines=new_lines+1;
        new_lines=[1 new_lines]; %#ok add the first text character
        Enzyme(i).KEGG_ID=enzymeList(i).KEGG_ID;
        %find entry
        index_position=regexp(raw_data,'ENTRY'); %find the position of ENTRY in txt
        k=find(new_lines==index_position); %find th position where the position of ENTRY=f(where ENTRY starts in a new line)
        Enzyme(i).ENTRY=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'ENTRY','');%replace(txt(from postion f(k):before the start of new line)-> replace entry with empty space(remove entry from the str)
        Enzyme(i).ENTRY=regexprep(Enzyme(i).ENTRY,'\s+',' ');
        Enzyme(i).ENTRY=regexprep(Enzyme(i).ENTRY,'EC',' ');
        Enzyme(i).ENTRY=strtrim(Enzyme(i).ENTRY);
        r=regexp(Enzyme(i).ENTRY,' ');
        Enzyme(i).ENTRY=Enzyme(i).ENTRY(1:r-1);
        %names
        index_position=regexp(raw_data,'NAME');%find the position of NAME in txt
        if ~isempty(index_position)%if the position of NAME is not empty-->
            k=find(new_lines==index_position(1));%K=position where NAME=new line, find the position where the position of NAME=f(where NAME starts in a new line) *this lines give me error 'Matrix dimensions must agree.'*
            if isempty(k)%==1
                k=find(new_lines<index_position);%K=find psotion where new line is greater than position of NAME in text 
                k=k(end);
            else
                
            end
            Enzyme(i).NAME=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'NAME','');
            Enzyme(i).NAME=regexprep(Enzyme(i).NAME,'\s+',' ');
            Enzyme(i).NAME=split(Enzyme(i).NAME,{'; '});
            Enzyme(i).NAME=strtrim(Enzyme(i).NAME);
        end   
        %CLASS
        index_position=regexp(raw_data,'CLASS');
        if isempty(index_position)==0
            index_position=index_position(1);
            k=find(new_lines==index_position);
            Enzyme(i).CLASS=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'CLASS','');
            Enzyme(i).CLASS=regexprep(Enzyme(i).CLASS,'\s+',' ');
            Enzyme(i).CLASS=strtrim(Enzyme(i).CLASS);
       %     m(i).ORGANISM=split(m(i).GENES,{'map'});
        %    m(i).ORGANISM=insertBefore(m(i).GENES,1,'map');
            %for k=1:length(m(i).RCLASS)
            %   m(i).RCLASS{k}=m(i).RCLASS{k}(3:length(m(i).RCLASS));
            %end
            %m(i).RCLASS=split(m(i).CLASS,'RC');
        end
        %REACTION
        index_position=regexp(raw_data,'ALL_REAC');
        if isempty(index_position)==0
            %index_position=index_position;
            k=find(new_lines==index_position); 
            Enzyme(i).REACTION=regexprep(raw_data(new_lines(k)+8:(new_lines(k+1)-1)),'REACTION','');
            Enzyme(i).REACTION=regexprep(Enzyme(i).REACTION,'>',' ');
            Enzyme(i).REACTION=regexprep(Enzyme(i).REACTION,'(other)',' ');
            Enzyme(i).REACTION=regexprep(Enzyme(i).REACTION,')',' ');
            Enzyme(i).REACTION=regexprep(Enzyme(i).REACTION,'(',' ');
            Enzyme(i).REACTION=regexprep(Enzyme(i).REACTION,';',' ');
            Enzyme(i).REACTION=regexprep(Enzyme(i).REACTION,'\s+',' ');
            Enzyme(i).REACTION=strtrim(Enzyme(i).REACTION);
            Enzyme(i).REACTION=split(Enzyme(i).REACTION,' ');
            %b=regexp(Enzyme(i).REACTION,'RN:R');
            %c=regexp(Enzyme(i).REACTION,']');
            %Enzyme(i).REACTION=Enzyme(i).REACTION(b+4:c-1);
            %Enzyme(i).REACTION=split(Enzyme(i).REACTION,'R');
        end
         %Substrate
        index_position=regexp(raw_data,'SUBSTRATE');
        if isempty(index_position)==0
            index_position=index_position(1);
            k=find(new_lines==index_position);
            Enzyme(i).SUBSTRATE=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'SUBSTRATE','');
            Enzyme(i).SUBSTRATE=regexprep(Enzyme(i).SUBSTRATE,'\s+',' ');
            Enzyme(i).SUBSTRATE=strtrim(Enzyme(i).SUBSTRATE);
            Enzyme(i).SUBSTRATE=split(Enzyme(i).SUBSTRATE,'+');
            for k=1:numel(Enzyme(i).SUBSTRATE)
                index_of_CPD=regexp(Enzyme(i).SUBSTRATE{k,1},'CPD:');
                index_bracket=regexp(Enzyme(i).SUBSTRATE{k,1},']');
                Enzyme(i).SUBSTRATE{k,1}=Enzyme(i).SUBSTRATE{k,1}(index_of_CPD+4:index_bracket-1);
            end
        %    index_of_CPD=regexp(Enzyme(i).SUBSTRATE,'CPD:');
        %    number_of_indexes=numel(index_of_CPD);
        %    index_bracket=regexp(m(i).SUBSTRATE,']');
        %    for x=1:number_of_indexes
        %        Enzyme(i).SUBSTRATE=Enzyme(i).SUBSTRATE(index_of_CPD(x)+4:index_of_CPD(x)+7);
        %    end
        end
        %PRODUCT
        index_position=regexp(raw_data,'PRODUCT');
        if isempty(index_position)==0
            index_position=index_position(1);
            k=find(new_lines==index_position);
            Enzyme(i).PRODUCT=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'PRODUCT','');
            Enzyme(i).PRODUCT=regexprep(Enzyme(i).PRODUCT,'\s+',' ');
            Enzyme(i).PRODUCT=strtrim(Enzyme(i).PRODUCT);
            Enzyme(i).PRODUCT=split(Enzyme(i).PRODUCT,';');
            for k=1:numel(Enzyme(i).PRODUCT)
                index_of_CPD=regexp(Enzyme(i).PRODUCT{k,1},'CPD:');
                index_bracket=regexp(Enzyme(i).PRODUCT{k,1},']');
                Enzyme(i).PRODUCT{k,1}=Enzyme(i).PRODUCT{k,1}(index_of_CPD+4:index_bracket-1);
            end
        %    index_of_CPD=regexp(Enzyme(i).SUBSTRATE,'CPD:');
        %    number_of_indexes=numel(index_of_CPD);
        %    index_bracket=regexp(m(i).SUBSTRATE,']');
        %    for x=1:number_of_indexes
        %        Enzyme(i).SUBSTRATE=Enzyme(i).SUBSTRATE(index_of_CPD(x)+4:index_of_CPD(x)+7);
        %    end
        end
        %Genes
        index_position=regexp(raw_data,'GENES');
        if isempty(index_position)==0
            k=find(new_lines==index_position);
            Enzyme(i).GENES=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'GENES','');
            %m(i).GENES=regexprep(m(i).GENES,'\s+',' ');
            %m(i).GENES=split(m(i).GENES,{'?'});
            %s=regexp(m(i).GENES,'[\n] ');
            %Enzyme(i).GENES=split(Enzyme(i).GENES,{'         '});
            Enzyme(i).GENES=strtrim(Enzyme(i).GENES);
            if strcmp(Enzyme(i).GENES(1:3),'HSA')
                find_space=regexp(Enzyme(i).GENES,'  ');
                b=Enzyme(i).GENES(1:find_space(1));
                Enzyme(i).GENES=erase(eraseBetween(b,'(',')'),for_gene_match);
                Enzyme(i).GENES=strtrim(Enzyme(i).GENES);
                Enzyme(i).GENES=split(Enzyme(i).GENES,' ',1);
            else
                Enzyme(i).GENES=[];
            end
                    
        end
        %Organisms
        index_position=regexp(raw_data,'GENES');
        if isempty(index_position)==0
            index_position=index_position(1);
            k=find(new_lines==index_position); 
            Enzyme(i).ORGANISM=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'GENES','');
            %m(i).ORGANISM=regexprep(m(i).GENES,'\s+',' ');
            %m(i).ORGANISM=split(m(i).GENES,{'map'});
            %m(i).ORGANISM=insertBefore(m(i).GENES,1,'map');
            Enzyme(i).ORGANISM=split(Enzyme(i).ORGANISM,{'         '});
            Enzyme(i).ORGANISM=strtrim(Enzyme(i).ORGANISM);
            for x=1:length(Enzyme(i).ORGANISM)
                position_of_semicolon=regexp(Enzyme(i).ORGANISM{x,1},':');
                Enzyme(i).ORGANISM{x,1}=Enzyme(i).ORGANISM{x,1}(1:position_of_semicolon-1);
            end
        end
        %ENZYME
        index_position=regexp(raw_data,'ENZYME');
        if isempty(index_position)==0
            k=find(new_lines==index_position(1));
            Enzyme(i).ENZYME=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'ENZYME','');
            Enzyme(i).ENZYME=regexprep(Enzyme(i).ENZYME,'\s+',' ');
            Enzyme(i).ENZYME=split(Enzyme(i).ENZYME,' ');
        end
    end
    close(handleWaitbar)
end