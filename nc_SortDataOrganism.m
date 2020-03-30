%organismList=ot([1:5892]);
function Organism_Taxonomy=nc_SortDataOrganism(organismList)
    Organism_Taxonomy=struct([]);
    numOrganism=length(organismList);
    handleWaitbar=waitbar(0,'Please wait...');
    for i=1:numOrganism
        waitbar(i/numOrganism,handleWaitbar,['Processing ' num2str(i) ' of ' num2str(numOrganism)]);
        raw_data=organismList(i).txt; % get specific entry based on ind (index)
        %new_lines=regexp(raw_data,'[\n]'); % find all new lines
        %new_lines_with_space=regexp(raw_data,'[\n] '); % find all new lines where the next line starts with a space
        %new_lines(ismember(new_lines,new_lines_with_space))=[]; % remove the ones that start with a space from the list in f
        %new_lines=new_lines+1;
        %new_lines=[1 new_lines]; % add the first text character
        Organism_Taxonomy(i).KEGG_ID=organismList(i).KEGG_ID;
        %find entry
        %ind=regexp(txt,'ENTRY'); %find the position of ENTRY in txt
        %k=find(f==ind); %find th position where the position of ENTRY=f(where ENTRY starts in a new line)
        %m(i).ENTRY=regexprep(txt(f(k):(f(k+1)-1)),'ENTRY','');%replace(txt(from postion f(k):before the start of new line)-> replace entry with empty space(remove entry from the str)
        %m(i).ENTRY=regexprep(m(i).ENTRY,'\s+',' ');
        %m(i).ENTRY=strtrim(m(i).ENTRY);
        %m(i).ENTRY=m(i).ENTRY(2:6);
        for k=1:numel(raw_data)
            Organism_Taxonomy(i).ENTRY=strtrim(raw_data(1:6));
            position_of_space=regexp(raw_data,' ');  %doesnt work cuz uses bracket to separate name
            Organism_Taxonomy(i).NAME=strtrim(raw_data(position_of_space(1):position_of_space(end)));
            position_of_spaces=regexp(Organism_Taxonomy(i).NAME,' ');
            if ~isempty(Organism_Taxonomy(i).NAME)
                Organism_Taxonomy(i).NAME_abbreviation=Organism_Taxonomy(i).NAME(1:position_of_spaces(1));
            end
            Organism_Taxonomy(i).TAXONOMY=strtrim(raw_data(position_of_space+1:numel(raw_data)));
            Organism_Taxonomy(i).TAXONOMY=split(Organism_Taxonomy(i).TAXONOMY,';');
        end
    end
    close(handleWaitbar)
end