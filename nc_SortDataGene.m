%human_gene_List=gt([1:3641]);
%human_gene_List=gt([1:10 100:110 1000:1010]);
function Human_genes=nc_SortDataGene(human_gene_List)
    Human_genes=struct([]);
    number_of_Genes=length(human_gene_List);
    handleWaitbar=waitbar(0,'Please wait...');
    for i=1:number_of_Genes
        waitbar(i/number_of_Genes,handleWaitbar,['Processing ' num2str(i) ' of ' num2str(number_of_Genes)]);
        raw_data=human_gene_List(i).txt; % get specific entry based on ind (index)
        enzyme_from_geneList=human_gene_List(i).ec;
        new_lines=regexp(raw_data,'[\n]'); % find all new lines
        new_lines_with_space=regexp(raw_data,'[\n] '); % find all new lines where the next line starts with a space
        new_lines(ismember(new_lines,new_lines_with_space))=[]; % remove the ones that start with a space from the list in f
        new_lines=new_lines+1;
        new_lines=[1 new_lines]; %#ok add the first text character
        Human_genes(i).KEGG_ID=human_gene_List(i).KEGG_ID;
        %find entry
        index=regexp(raw_data,'ENTRY'); %find the position of ENTRY in txt
        k=find(new_lines==index(1)); %find th position where the position of ENTRY=f(where ENTRY starts in a new line)
        Human_genes(i).ENTRY=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'ENTRY','');%replace(txt(from postion f(k):before the start of new line)-> replace entry with empty space(remove entry from the str)
        Human_genes(i).ENTRY=regexprep(Human_genes(i).ENTRY,'\s+',' ');
        Human_genes(i).ENTRY=regexprep(Human_genes(i).ENTRY,'EC',' ');
        Human_genes(i).ENTRY=strtrim(Human_genes(i).ENTRY);
        position_of_space=regexp(Human_genes(i).ENTRY,' ');
        Human_genes(i).ENTRY=Human_genes(i).ENTRY(1:position_of_space-1);
        %names
        index=regexp(raw_data,'NAME');
        if ~isempty(index)
            k=find(new_lines==index(1));
            if isempty(k)%==1
                k=find(new_lines<index);
                k=k(end);
            else
                
            end
            Human_genes(i).NAME=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'NAME','');
            Human_genes(i).NAME=regexprep(Human_genes(i).NAME,'\s+',' ');
            Human_genes(i).NAME=split(Human_genes(i).NAME,{','});
            Human_genes(i).NAME=strtrim(Human_genes(i).NAME);
        end     
        %DEFINITION
        index=regexp(raw_data,'DEFINITION');%find the position of NAME in txt
        if ~isempty(index)%if the position of NAME is not empty-->
            k=find(new_lines==index(1));%K=position where NAME=new line, find the position where the position of NAME=f(where NAME starts in a new line) *this lines give me error 'Matrix dimensions must agree.'*
            if isempty(k)%==1
                k=find(new_lines<index);%K=find psotion where new line is greater than position of NAME in text 
                k=k(end); 
            else
                
            end
            Human_genes(i).DEFINITION=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'DEFINITION','');
            Human_genes(i).DEFINITION=regexprep(Human_genes(i).DEFINITION,'\s+',' ');
            %m(i).DEFINITION=split(m(i).DEFINITION,{'; '});
            %m(i).DEFINITION=strtrim(m(i).DEFINITION);
        end   
        %ENZYME
        Human_genes(i).ENZYME=enzyme_from_geneList(4:numel(enzyme_from_geneList));
        Human_genes(i).ENZYME=strtrim(Human_genes(i).ENZYME);
        %ORTHOLOGY
        index=regexp(raw_data,'ORTHOLOGY');
        if isempty(index)==0
            index=index(1);
            k=find(new_lines==index);
            Human_genes(i).ORTHOLOGY=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'ORTHOLOGY','');
            Human_genes(i).ORTHOLOGY=regexprep(Human_genes(i).ORTHOLOGY,'\s+',' ');
            Human_genes(i).ORTHOLOGY=strtrim(Human_genes(i).ORTHOLOGY);
            %m(i).ORTHOLOGY=split(m(i).ORTHOLOGY,{''});
        %    m(i).ORGANISM=insertBefore(m(i).GENES,1,'map');
            %for k=1:length(m(i).RCLASS)
            %   m(i).RCLASS{k}=m(i).RCLASS{k}(3:length(m(i).RCLASS));
            %end
            %m(i).RCLASS=split(m(i).CLASS,'RC');
        end
        %ORGANISM
        index=regexp(raw_data,'ORGANISM');
        if isempty(index)==0
            index=index(1);
            k=find(new_lines==index);
            Human_genes(i).ORGANISM=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'ORGANISM','');
            Human_genes(i).ORGANISM=regexprep(Human_genes(i).ORGANISM,'\s+',' ');
            Human_genes(i).ORGANISM=strtrim(Human_genes(i).ORGANISM);
            %m(i).ORGANISM=split(m(i).ORGANISM,{' '});
        end
        %PATHWAY
        index=regexp(raw_data,'PATHWAY');
        if ~isempty(index)
            k=find(new_lines==index);
            if isempty(k)%==1
                k=find(new_lines<index);
                k=k(end);
            
            end
            Human_genes(i).PATHWAY=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'PATHWAY','');
            Human_genes(i).PATHWAY=regexprep(Human_genes(i).PATHWAY,'\s+',' ');
            Human_genes(i).PATHWAY=strtrim(Human_genes(i).PATHWAY);
            Human_genes(i).PATHWAY=strtrim(Human_genes(i).PATHWAY);
            Human_genes(i).PATHWAY=Human_genes(i).PATHWAY(4:numel(Human_genes(i).PATHWAY));
            Human_genes(i).PATHWAY=split(Human_genes(i).PATHWAY,'hsa');
            for k=1:numel(Human_genes(i).PATHWAY)
                position_of_space=regexp(Human_genes(i).PATHWAY{k},' ');
                Human_genes(i).PATHWAY{k}=Human_genes(i).PATHWAY{k}(1:position_of_space(1)-1);
            end
        end
        %BRITE
        index=regexp(raw_data,'BRITE');
        if isempty(index)==0
            index=index(1);
            k=find(new_lines==index);
            Human_genes(i).BRITE=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'BRITE','');
            Human_genes(i).BRITE=regexprep(Human_genes(i).BRITE,'\s+',' ');
            Human_genes(i).BRITE=strtrim(Human_genes(i).BRITE);
           % m(i).BRITE=split(m(i).BRITE,{' '});
        end
        %  POSITION
        index=regexp(raw_data,'POSITION');
        if isempty(index)==0
            k=find(new_lines==index);
            Human_genes(i).POSITION=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'POSITION','');
            Human_genes(i).POSITION=regexprep(Human_genes(i).POSITION,'\s+',' ');
            %m(i).GENES=split(m(i).GENES,{'?'});
            %s=regexp(m(i).GENES,'[\n] ');
             %m(i).ORTHOLOGY=split(m(i).ORTHOLOGY,{'            '});
             Human_genes(i).POSITION=strtrim(Human_genes(i).POSITION);
        end
        %MOTIF
        index=regexp(raw_data,'MOTIF');
        if isempty(index)==0
            index=index(1);
            k=find(new_lines==index);
            Human_genes(i).MOTIF=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'MOTIF','');
            Human_genes(i).MOTIF=regexprep(Human_genes(i).MOTIF,'\s+',' ');
            Human_genes(i).MOTIF=strtrim(Human_genes(i).MOTIF);
           % m(i).BRITE=split(m(i).BRITE,{' '});
        end
        %DBLINKS
        index=regexp(raw_data,'DBLINKS');
        if isempty(index)==0
            index=index(1);
            k=find(new_lines==index);
            Human_genes(i).DBLINKS=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'DBLINKS','');
            Human_genes(i).DBLINKS=regexprep(Human_genes(i).DBLINKS,'\s+',' ');
            Human_genes(i).DBLINKS=strtrim(Human_genes(i).DBLINKS);
           % m(i).BRITE=split(m(i).BRITE,{' '});
        end
        %STRUCTURE
        index=regexp(raw_data,'STRUCTURE');
        if isempty(index)==0
            index=index(1);
            k=find(new_lines==index); 
            Human_genes(i).STRUCTURE=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'STRUCTURE','');
            Human_genes(i).STRUCTURE=regexprep(Human_genes(i).STRUCTURE,'\s+',' ');
            Human_genes(i).STRUCTURE=strtrim(Human_genes(i).STRUCTURE);
           % m(i).BRITE=split(m(i).BRITE,{' '});
        end
        %AASEQ
        index=regexp(raw_data,'AASEQ');
        if isempty(index)==0
            index=index(1);
            k=find(new_lines==index); 
            Human_genes(i).AASEQ=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'AASEQ','');
            Human_genes(i).AASEQ=regexprep(Human_genes(i).AASEQ,'\s+',' ');
            Human_genes(i).AASEQ=strtrim(Human_genes(i).AASEQ);
           % m(i).BRITE=split(m(i).BRITE,{' '});
        end
        %NTSEQ
        index=regexp(raw_data,'NTSEQ');
        if isempty(index)==0
            index=index(1);
            k=find(new_lines==index); 
            Human_genes(i).NTSEQ=regexprep(raw_data(new_lines(k):(new_lines(k+1)-1)),'NTSEQ','');
            Human_genes(i).NTSEQ=regexprep(Human_genes(i).NTSEQ,'\s+',' ');
            Human_genes(i).NTSEQ=strtrim(Human_genes(i).NTSEQ);
           % m(i).BRITE=split(m(i).BRITE,{' '});
        end
    end
    close(handleWaitbar)
end