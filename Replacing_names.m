
function [Metabolites,Deleted_entry]=Replacing_names(Metabolites,List_of_index)
tic;    
handleWaitbar=waitbar(0,'Please wait...');
List_to_delete=[]; 
List_to_watch={};
    for i=size(List_of_index,1):-1:1
        waitbar(i/size(List_of_index,1),handleWaitbar,['Renaming metabolites: ' num2str(i) ' of ' num2str(size(List_of_index,1))]);
        if size(List_of_index{i},1)>1
            first=List_of_index{i}(1);
            
            for k=size(List_of_index{i,1},1):-1:2
                j=List_of_index{i}(k);
                OGN=Metabolites(first).NAME;
                OGR=Metabolites(first).REACTION;
                OGP=Metabolites(first).PATHWAY;
                OGE=Metabolites(first).ENZYME;
                List_to_watch{end+1}=OGN;%#ok
                if size(Metabolites(j).NAME,2)>1
                    Metabolites(first).NAME=cat(1,OGN,Metabolites(j).NAME(1:end,1),Metabolites(j).NAME(1:end,2));
                else
                    Metabolites(first).NAME=cat(1,OGN,Metabolites(j).NAME);
                end 
                Metabolites(first).REACTION=cat(1,OGR,Metabolites(j).REACTION);
                Metabolites(first).PATHWAY=cat(1,OGP,Metabolites(j).PATHWAY);
                Metabolites(first).ENZYME=cat(1,OGE,Metabolites(j).ENZYME);
                Metabolites(first).NAME=unique(Metabolites(first).NAME);
                Metabolites(first).REACTION=unique(Metabolites(first).REACTION);
                Metabolites(first).PATHWAY=unique(Metabolites(first).PATHWAY);
                Metabolites(first).ENZYME=unique(Metabolites(first).ENZYME);
                List_to_delete(end+1,1)=j;%#ok
                
            
            end

                
        end
    end
    Deleted_entry={};
    List_to_delete=unique(List_to_delete);
    close(handleWaitbar)
    handleWaitbar=waitbar(0,'Please wait...');
    for k=size(List_to_delete,1):-1:1
        waitbar(k/size(List_to_delete,1),handleWaitbar,['Removing duplicate ' num2str(k) ' of ' num2str(size(List_to_delete,1))]);
        Deleted_entry{k,1}=Metabolites(List_to_delete(k,1)).KEGG_ID;    
        Metabolites(List_to_delete(k,1))=[];
    end
    Deleted_compounds=(Deleted_entry(cell2mat(Deleted_entry(:,1))=='C'));
    Deleted_glycans=(Deleted_entry(cell2mat(Deleted_entry(:,1))=='G'));
    close(handleWaitbar)
    disp(toc);
end



