% 1) Load all the databases
% 2) Enter buildNo(choose between 37 or 38)
% 3) Enter the diff(numerical) : diff is the distant in bp from the gene
%diff is the distant that the SNP will associated with the gene
%this code should take around 1 hour to run
function [Metabolites,Reaction,Class,Enzyme,Human_genes,Organism_Taxonomy,find_SNP_Gene,RPAIR_new_table,List_of_index]=nc_formDatabase(ct,et,gt,mt,ot,rt,intestinalcancer_rsPosition,Gene_Pos,BuildNo,diff)
t=tic;
handleWaitbar1=waitbar(0,'Initializing...');
ClassList=ct([1:3161]);%#ok
enzymeList=et([1:7637]);%#ok
human_gene_List=gt([1:3641]);%#ok
organismList=ot([1:5892]);%#ok
metabolitesList=mt([1:29649]);%#ok
reactionList=rt([1:11296]); %#ok
t_hms = datevec(toc(t)./(60*60*24));
waitbar(0/6,handleWaitbar1,['Sorting compound data (Part 1 of 6). Time elapsed: ' ...
    num2str(round(t_hms(4))) 'h ' num2str(round(t_hms(5))) 'm ' num2str(round(t_hms(6))) 's']);
Metabolites=nc_SortDataCompounds(metabolitesList);
Class=nc_SortDataClass(ClassList);
Enzyme=nc_SortDataEnzyme(enzymeList);
Human_genes=nc_SortDataGene(human_gene_List);
Reaction=nc_SortDataReaction(reactionList);
Organism_Taxonomy=nc_SortDataOrganism(organismList);
find_SNP_Gene=nc_SortDataGene_SNP(Gene_Pos,intestinalcancer_rsPosition,BuildNo,diff);
t_hms = datevec(toc(t)./(60*60*24));
waitbar(1/6,handleWaitbar1,['Expanding metabolite names (Part 2 of 6). Time elapsed: ' ...
    num2str(round(t_hms(4))) 'h ' num2str(round(t_hms(5))) 'm ' num2str(round(t_hms(6))) 's']);
List_of_index=Expanding_Metabolite_names(Metabolites);
t_hms = datevec(toc(t)./(60*60*24));
waitbar(2/6,handleWaitbar1,['Replacing names (Part 3 of 6). Time elapsed: ' ...
    num2str(round(t_hms(4))) 'h ' num2str(round(t_hms(5))) 'm ' num2str(round(t_hms(6))) 's']);
[Metabolites]=Replacing_names(Metabolites,List_of_index);
t_hms = datevec(toc(t)./(60*60*24));
waitbar(3/6,handleWaitbar1,['Removing duplicate reactions (Part 4 of 6). Time elapsed: ' ...
    num2str(round(t_hms(4))) 'h ' num2str(round(t_hms(5))) 'm ' num2str(round(t_hms(6))) 's']);
Reaction=Removing_duplicates_reaction(Reaction,List_of_index);
t_hms = datevec(toc(t)./(60*60*24));
waitbar(4/6,handleWaitbar1,['Pairing compounds-reaction (Part 5 of 6). Time elapsed: ' ...
    num2str(round(t_hms(4))) 'h ' num2str(round(t_hms(5))) 'm ' num2str(round(t_hms(6))) 's']);
[RPAIR_new_table,Reaction]=Pair_Reaction_compound(Reaction);
t_hms = datevec(toc(t)./(60*60*24));
waitbar(5/6,handleWaitbar1,['Pairing compound ID (Part 6 of 6). Time elapsed: ' ...
    num2str(round(t_hms(4))) 'h ' num2str(round(t_hms(5))) 'm ' num2str(round(t_hms(6))) 's']);
[Metabolites,RPAIR_new_table]=Pair_ID_compound(Metabolites,RPAIR_new_table);
clear ct et gt mt ot rt intestinalcancer_rsPosition Gene_Pos ClassList enzymeList human_gene_List organismList metabolitesList reactionList;
t_hms = datevec(toc(t)./(60*60*24));
waitbar(6/6,handleWaitbar1,['Forming database completed. Time elapsed: ' ...
    num2str(round(t_hms(4))) 'h ' num2str(round(t_hms(5))) 'm ' num2str(round(t_hms(6))) 's']);
end



%what do I do with enzymes and gene data-> this will be connected to
%reactions via enzyme code, and enzyme code will connect to gene by
%geneID-->SNP locations

%pathway can be taken from both metabolites and reactions but since one
%RPAIR represents multiple reactions, and each reaction is associated to
%multiple enzymes-->GENE-->SNP
%how do I link the reaction to enyzme to genes(PID_?)


%put sortenzyme and sortgene in the workflow


%how metabonetwork1.0 is improved-->
%number of compounds(mainly about database)(show host genes compare to
%microbial(which reaction) in related to snps and traits))
%need to find out how metabonetwork 1 and 2 are different -->what cant i do
%with metabonetwork one compare to metabonetwork 2 ????


%literature search and find how many studies for colonic neoplasm,colorectal adenoma, colorectal cancer, colorectal mucinous adenocarcinoma
%metastatic colorectal cancer, rectum cancer 
%synonyms for metabolomics, metabolic profiling/fingerprinting/phenotyping
%and put it all in an excel/ count how many hits i get per seach term and
%then record the search hits for each term
%look and make an excel sheet for all the metabolite for that search
%term/papers

%
%