
function [find_SNP_Gene,number_of_snp]=nc_SortDataGene_SNP(Gene_Pos,intestinalcancer_rsPosition,intestinalcancer_GWAS,diff)
%Find all matching Chromosome between Gene_Pos and inntestinalcancer_rsPosition
find_SNP_Gene=struct([]);
BuildNo=38;
gwas_ID=intestinalcancer_GWAS(:,1);
c={};
number_of_snp={};
if BuildNo==37
    ChrCol=2;
    StartCol=3;
    EndCol=4;
elseif BuildNo==38
    ChrCol=5;
    StartCol=6;
    EndCol=7;
else
    disp('Please enter Build 37 or 38');
    return
end
if isnumeric(diff)==1
    if rem(diff,1)==0 || diff==0
        chr=sortrows(Gene_Pos,StartCol);
        chr=sortrows(chr,ChrCol);
            for i=1:size(chr,1)
                g=(intestinalcancer_rsPosition(cell2mat(intestinalcancer_rsPosition(:,3))>=(chr(i,StartCol)-diff)&...
                    cell2mat(intestinalcancer_rsPosition(:,3))<=(chr(i,EndCol)+diff)&...
                    cell2mat(intestinalcancer_rsPosition(:,2))==chr(i,ChrCol),:));
                g=table2cell(unique(cell2table(g),'rows'));
                find_SNP_Gene(i).GeneID=chr(i,1);
                find_SNP_Gene(i).Chr=chr(i,ChrCol);
                find_SNP_Gene(i).GeneStart=chr(i,StartCol);
                find_SNP_Gene(i).GeneEnd=chr(i,EndCol);
                find_SNP_Gene(i).SNP=g;
                find_SNP_Gene(i).trait={};                 
            end
            for k=1:3641
                if ~isempty(find_SNP_Gene(k).SNP)
                    gwas_ID_SNP={find_SNP_Gene(k).SNP{:,1}}';
                    c=find(ismember(gwas_ID,gwas_ID_SNP));
                    for j=1:numel(c)
                        find_SNP_Gene(k).trait{end+1,1}=intestinalcancer_GWAS{c(j),2};
                    end
                    if size(find_SNP_Gene(k).SNP,1)>1
                        for j=1:size(find_SNP_Gene(k).SNP,1)
                            number_of_snp{end+1,1}=find_SNP_Gene(k).SNP(j,1:4);
                        end
                    else
                        number_of_snp{end+1,1}=find_SNP_Gene(k).SNP;
                    end
                end
            end
    elseif rem(diff,1)~=0 || diff<0
        disp('Please enter positive integers ;0');
    end
else
    disp('Please enter positive integers :0');
end
end