
function [Metabolites,RPAIR_new_table]=Pair_ID_compound(Metabolites,RPAIR_new_table)        
    tic
    handleWaitbar=waitbar(0,'Pairing pID to Metabolites');
    for i=1:numel(Metabolites)
        Metabolites(i).pID={};
    end
    cpd={Metabolites.KEGG_ID}';
    for k=1:size(RPAIR_new_table,2)
        waitbar(k/size(RPAIR_new_table,2),handleWaitbar,['Pairing ID: ' num2str(k) ' of ' num2str(size(RPAIR_new_table,2))]);
        for j=1:size(RPAIR_new_table(k).COMPOUND,2)
            RPAIR_new_table(k).ADJACENCY_COMPOUND{j}=find(ismember(cpd,RPAIR_new_table(k).COMPOUND(j)),1);
            a=RPAIR_new_table(k).ADJACENCY_COMPOUND{1,j};
            if ~isempty(Metabolites(a).pID)
                Metabolites(a).pID{end+1,1}=RPAIR_new_table(k).ID;
            else
                Metabolites(a).pID{end+1,1}=RPAIR_new_table(k).ID;
            end
                Metabolites(a).pID=unique(Metabolites(a).pID);%do you want this to be unique because there are 2 compounds for metabolites
        end
    end
    delete(handleWaitbar)
    toc
end
    