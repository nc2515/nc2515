function [ES,NS,Table]=Annotating_MetaboPathways(RPAIR_new_table,Enzyme,Reaction,find_SNP_Gene,ES,NS,Metabolites)
CP={};
CPi=[];
CPSNP={};
CPGENE={};
CPEnz={};
Table={};
cID_list={Metabolites.KEGG_ID};
cID_list=cID_list';
Index_for_K=[];%#ok
match=["[","]"];%#ok
RPAIR={RPAIR_new_table.CPAIR}';
Enzyme_list={Enzyme.ENTRY}';
GeneID_list=cell2mat({find_SNP_Gene.GeneID})';
tic;
handleWaitbar=waitbar(0,'Annotating...');
size_of_ES=size(ES,1);
%determine max number of possible reactions and SNPs
for i=1:length(RPAIR_new_table)
a(i)=length(RPAIR_new_table(i).REACTION); %#ok<AGROW>
end
max_reaction=max(a);
for i=1:length(find_SNP_Gene)
b(i)=length(find_SNP_Gene(i).SNP); %#ok<AGROW>
end
max_SNP=max(b);
%create a colour map from black to red for number of SNPs
red = [1,0,0];
darkred=[0.5,0,0];
blue=[0,0,1];
darkblue=[0,0,0.5];
black = [0,0,0];
colors_red = [linspace(darkred(1),red(1),max_SNP+1)', linspace(darkred(2),red(2),max_SNP+1)', linspace(darkred(3),red(3),max_SNP+1)'];%#ok
colors_blue= [linspace(darkblue(1),blue(1),max_SNP+1)', linspace(darkblue(2),blue(2),max_SNP+1)', linspace(darkblue(3),blue(3),max_SNP+1)'];%#ok
cols=jet(10);
CN=[cols(1,1),cols(1,2),cols(1,3)];%'colonic neoplasm' dark blue 0,0,0.67   255,255,170
CA=[1,0,0];%'colorectal adenoma' red
CRC=[cols(3,1),cols(3,2),cols(3,3)];%'colorectal cancer'ligth blue  0,0.33,1 0,85,255
CMA=[cols(4,1),cols(4,2),cols(4,3)];%'colorectal mucinous adenocarcinoma'ligther blue 0,0.67,1 0,170,255
GA=[cols(5,1),cols(5,2),cols(5,3)];%'gastric adenocarcinoma' very light blue 0,1,1 0,255,255
GC=[cols(6,1),cols(6,2),cols(6,3)];%'gastric carcinoma' light green .33,1,.67, 85,255,170
GCC=[cols(7,1),cols(7,2),cols(7,3)];%'gastric cardia carcinoma') yellow green .67,1,.33 170,255,85
MCRC=[0,0,1];%'metastatic colorectal cancer') blue 0,0,1  0,0,255
SINT=[cols(9,1),cols(9,2),cols(9,3)];%'small intestine neuroendocrine tumor') orange 1,0.67 255,170,0
mul_traits=[cols(10,1),cols(10,2),cols(10,3)]; %yellow, 1,1,0 255,255,0
index_traits={};
%gradient for line width based on max #reactions
x=[1:max_reaction]; %#ok<NBRAK>
xlog=log10(1+(9*x/max(x)));
x_w=1+2*xlog;
for i=1:size(ES,1)
    number_of_reaction=[];
    number_of_SNP=[];
    human_reaction=[];
    waitbar(i/size_of_ES,handleWaitbar,['Annotating: ' num2str(i) ' of ' num2str(size_of_ES) ' (for' num2str(toc)]);
    Index_for_K=str2double(strsplit(regexprep(ES(i).strUD,'[[]]',''))); 
    %Index_for_K=Index_for_K';
    %Node1=str2double(cell2mat(Index_for_K{i}(1)));
    %Node2=str2double(cell2mat(Index_for_K{i}(2)));
    if Index_for_K(1)~=0 && Index_for_K(2)~=0
        if str2double(NS(Index_for_K(2)).K{1}(2:6))>str2double(NS(Index_for_K(1)).K{1}(2:6))
            cp=[NS(Index_for_K(1)).K{1} '_' NS(Index_for_K(2)).K{1}];
        else
            cp=[NS(Index_for_K(2)).K{1} '_' NS(Index_for_K(1)).K{1}];
        end
        Index_for_PID=find(strcmp(cp,RPAIR));
        if ~isempty(Index_for_PID)
                number_of_reaction=size(RPAIR_new_table(Index_for_PID).REACTION,2);%display this data on MBP somehow(thickerlines)
                Index_for_Reaction=RPAIR_new_table(Index_for_PID).ADJACENCY_REACTION;
                Index_for_Reaction=Index_for_Reaction{1};
                        number_of_enzymes=size(Reaction(Index_for_Reaction).ENZYME,2);%#ok<NASGU> %Number of enzymes that can do this reactions
                Enzyme_ID=Reaction(Index_for_Reaction).ENZYME;
                if ~isempty(Enzyme_ID)
                    Index_for_Enzyme=find(ismember(Enzyme_list,Enzyme_ID));
                    if  ~isempty(Index_for_Enzyme)
                        index_traits={};
                        for k=1:size(Index_for_Enzyme,1)
                            if ~isempty(Enzyme(Index_for_Enzyme(k)).GENES) % GENES is not empty (i.e. human)
                                Gene_IDs=str2double(Enzyme(Index_for_Enzyme(k)).GENES);%what todo if there are more than 1 gene
                                Index_for_Gene=find(ismember(GeneID_list,Gene_IDs));
                                n=0;
                                for j=1:size(Index_for_Gene,1)
                                    n=n+size(find_SNP_Gene(Index_for_Gene(j)).SNP,1);
                                    
                                    if ~isempty(find_SNP_Gene(Index_for_Gene(j)).SNP)
                                        index_traits{end+1}=Index_for_Gene(j); %#ok
                                        CP{end+1,1}=cp;
                                        CPi(end+1,1)=size(CP,1);
                                        CPSNP{end+1,1}=find_SNP_Gene(Index_for_Gene(j)).SNP;
                                        CPGENE{end+1,1}=Enzyme(Index_for_Enzyme(k)).GENES;
                                        CPEnz{end+1,1}=Enzyme(Index_for_Enzyme(k)).KEGG_ID;
                                    end
                                end
                                number_of_SNP(k)=n; %#ok<AGROW>
                                if ~isempty(Enzyme(Index_for_Enzyme(k)).ORGANISM)&& strcmp(Enzyme(Index_for_Enzyme(k)).ORGANISM(1),'HSA')
                                    human_reaction(k)=1; %#ok<AGROW>
                                else %if ORGANISM is not HSA (i.e. not human)
                                    human_reaction(k)=0;%#ok<AGROW>
                                    number_of_SNP(k)=0;%#ok<AGROW>
                                end 
                            else %GENES is empty (i.e. not human)
                                human_reaction(k)=0;%#ok<AGROW>
                                number_of_SNP(k)=0;%#ok<AGROW>
                            end
                         end
                    end
                else %No enzyme
                    number_of_SNP=0;
                end
        end
        %change width based on #reaction
        if ~isempty(number_of_reaction)
            ES(i).LW=x_w(number_of_reaction);
        end
        %change colour based on human vs non-human
        if sum(human_reaction)==length(human_reaction) %only human genes
            ES(i).ES='-';
        elseif sum(human_reaction)==0 %only non-human genes
            ES(i).ES=':';
        else 
            ES(i).ES='--';
        end
        
        %change line style based on with/without SNP
            
            
        if sum(number_of_SNP)==0
            ES(i).C=black;
        elseif numel(index_traits)>1 && sum(number_of_SNP)>0
            ES(i).C=mul_traits;
        elseif strcmp(find_SNP_Gene(Index_for_Gene(j)).trait,'colonic neoplasm')
             ES(i).C=CN;
        elseif strcmp(find_SNP_Gene(Index_for_Gene(j)).trait,'colorectal adenoma')
             ES(i).C=CA;
        elseif strcmp(find_SNP_Gene(Index_for_Gene(j)).trait,'colorectal cancer')
             ES(i).C=CRC;
        elseif strcmp(find_SNP_Gene(Index_for_Gene(j)).trait,'colorectal mucinous adenocarcinoma')
             ES(i).C=CMA;
        elseif strcmp(find_SNP_Gene(Index_for_Gene(j)).trait,'gastric adenocarcinoma')
             ES(i).C=GA;
        elseif strcmp(find_SNP_Gene(Index_for_Gene(j)).trait,'gastric carcinoma')
             ES(i).C=GC;
        elseif strcmp(find_SNP_Gene(Index_for_Gene(j)).trait,'gastric cardia carcinoma')
             ES(i).C=GCC;
        elseif strcmp(find_SNP_Gene(Index_for_Gene(j)).trait,'metastatic colorectal cancer')
             ES(i).C=MCRC;
        elseif strcmp(find_SNP_Gene(Index_for_Gene(j)).trait,'small intestine neuroendocrine tumor')
             ES(i).C=SINT;
            
        end
    end
end
uCP=unique(CP);
for i=1:size(uCP,1)
    a=ismember(CP,uCP(i));
    Table{i,1}=uCP{i};
    Table{i,2}=CPGENE(a);
    Table{i,3}=CPSNP(a);
    Table{i,4}=CPEnz(a);
    Table{i,5}=Metabolites(find(ismember(cID_list,Table{i,1}(1:6)))).NAME;
    Table{i,6}=Metabolites(find(ismember(cID_list,Table{i,1}(8:13)))).NAME;
    
end
    
delete(handleWaitbar)
end