function [A,K,E,N,C,I]=create_database(Organism_Taxonomy,Enzyme,Class,Human_genes,Metabolites,Reaction,RPAIR_new_table)
Organism_Taxonomy=Organism_Taxonomy';
Enzyme=Enzyme';
Class=Class';
Human_genes=Human_genes';
Metabolites=Metabolites';
Reaction=Reaction';
RPAIR_new_table=RPAIR_new_table';

allPhyla={};
for k=1:size(Organism_Taxonomy,1)
    allPhyla{k,1}=Organism_Taxonomy(k).TAXONOMY{3};
end

microbeList={'Actinobacteria','Bacteroidetes','Cyanobacteria','Firmicutes - Bacilli','Firmicutes - Clostridia','Firmicutes - Others','Fusobacteria','Tenericutes','Verrucomicrobia','Alphaproteobacteria','Betaproteobacteria','Deltaproteobacteria','Gammaproteobacteria - Enterobacteria','Gammaproteobacteria','Gammaproteobacteria - Others','Epsilonproteobacteria'}'; % list of species to include

org_include=zeros(k,1);
org_include(1)=1; % inlcude Homo sapiens
org_include(ismember(allPhyla,microbeList))=1;

Metabolites(24086) = []; % remove glycan without a name - but this could mess up adding the right adjacencies, but chances are small as it's second-to-last
K={Metabolites.KEGG_ID}'; % list of compound/glycan codes (size: [n,1])

A=zeros(numel(K),numel(K)); % adjacency matrix

org={Organism_Taxonomy.NAME_abbreviation}'; % get all organism codes
for k=1:size(org,1)
    org{k}=upper(strtrim(org{k}));
end

E={}; % list of enzymes
HSA=[]; % human enzyme or not
for k=1:size(Enzyme,1)
    HSA(k,1)=0;
    if ~isempty(Enzyme(k).ORGANISM)
        if any(ismember(org,Enzyme(k).ORGANISM))
            E{end+1,1}=Enzyme(k).ENTRY;
            if ismember({'HSA'},Enzyme(k).ORGANISM)
                HSA(k,1)=1;
            end
        end
    end
end

HSA_E={Enzyme.ENTRY}';
HSA_E(HSA==0)=[]; % Homo sapiens enzymes only

include_pairs={};
for k=1:numel(Reaction)
    if ~isempty(Reaction(k).ENZYME)
        if any(ismember(Reaction(k).ENZYME,E))
            include_pairs=[include_pairs;Reaction(k).pID];
        end
    end
end
include_pairs=unique(include_pairs);

for k=1:size(RPAIR_new_table,1)
    if ismember({RPAIR_new_table(k).ID},include_pairs)
        ii=[RPAIR_new_table(k).ADJACENCY_COMPOUND{1} RPAIR_new_table(k).ADJACENCY_COMPOUND{2}];
        A(ii(1),ii(2))=1;
        A(ii(2),ii(1))=1;
    end
end

N={};
C={};
I=[];
for k=1:size(Metabolites,1)
    for j=1:size(Metabolites(k).NAME,1)
        if ~isempty(Metabolites(k).NAME{j,1})
            N=[N;Metabolites(k).NAME{j,1}];
            C=[C;Metabolites(k).KEGG_ID];
            I=[I;k];
        end
    end
end
end