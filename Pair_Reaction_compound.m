%RPAIR function
function [RPAIR_new_table,Reaction]=Pair_Reaction_compound(Reaction)
    RPAIR=struct([]);
    RPAIR_new_table=struct([]);
    Reaction_length=numel(Reaction);
    handleWaitbar=waitbar(0,'Pairing reaction...');
    for i=1:Reaction_length
        waitbar(i/Reaction_length,handleWaitbar,['Pairing reaction: ' num2str(i) ' of ' num2str(Reaction_length)]);
        RPAIR(i).ID=Reaction(i).ENTRY;
        RPAIR(i).RCLASS=Reaction(i).RCLASS;
        RPAIR(i).COMPOUND=RPAIR(i).RCLASS;
        %RPAIR(i).number_of_RPAIR=[];
        if ~isempty(RPAIR(i).RCLASS)
            RPAIR(i).COMPOUND=cellstr(split(string(RPAIR(i).RCLASS),'_'));
        end
        if ~isempty(RPAIR(i).COMPOUND)
            if numel(RPAIR(i).COMPOUND)==2
                RPAIR(i).COMPOUND{1,2}=RPAIR(i).COMPOUND{2,1};
                RPAIR(i).COMPOUND(2,:)=[];
            end
        else
        end
        num_row=size(RPAIR(i).COMPOUND);
        for k=1:num_row(1)
            sort(RPAIR(i).COMPOUND(k, :));
        end
        RPAIR(i).number_of_RPAIR=numel(RPAIR(i).RCLASS);
    end
    delete(handleWaitbar)
    RPAIR=RPAIR(find(cellfun(@(v)~isempty(v),{RPAIR.RCLASS}))); %#ok<FNDSB>
    
    pbase='P00000';
    rplist={};
    rcn={Reaction.KEGG_ID}';
    %cpd={Metabolites.KEGG_ID}';
    handleWaitbar=waitbar(0,'Compiling data...');
    for k=1:size(RPAIR,2)
        waitbar(k/size(RPAIR,2),handleWaitbar,['Compiling data: ' num2str(k) ' of ' num2str(size(RPAIR,2))]);
        for j=1:size(RPAIR(k).RCLASS,1)
            rc=RPAIR(k).RCLASS{j};
            if  ~ismember(string(rc),string(rplist))        %isempty(find(ismember(rplist,rc),1))
                pt=pbase;
                ind=size(RPAIR_new_table,2)+1;
                pt((end-length(num2str(ind))+1):end)=num2str(ind);
                RPAIR_new_table(ind).ID=pt;
                RPAIR_new_table(ind).CPAIR=rc;
                RPAIR_new_table(ind).COMPOUND=RPAIR(k).COMPOUND(j,1:2);
                RPAIR_new_table(ind).REACTION{1}=RPAIR(k).ID;
                rplist{ind,1}=rc; %#ok
            else 
                % add reaction to existing entry
                f = cellfun(@(x)ismember(string(x),string(rc)),rplist);%find index of where rc is the same as rplist
                RPAIR_new_table(find(f)).REACTION(end+1)={RPAIR(k).ID}; %#ok<FNDSB>  %may be due to multiple indexing
            end
            %RPAIR_new_table(k).ADJACENCY(j)=find(ismember(rcn,RPAIR_new_table(k).REACTION(j)),1);
        end
    end
    delete(handleWaitbar)
    handleWaitbar=waitbar(0,'Calculating adjacency...');
    for k=1:size(RPAIR_new_table,2)
        waitbar(k/size(RPAIR_new_table,2),handleWaitbar,['Calculating adjacency: ' num2str(k) ' of ' num2str(size(RPAIR_new_table,2))]);
        for j=1:size(RPAIR_new_table(k).REACTION,1)
            RPAIR_new_table(k).ADJACENCY_REACTION{j}=find(ismember(rcn,RPAIR_new_table(k).REACTION(j)),1);
        end
    end
    delete(handleWaitbar)
    [Reaction.pID]=deal([]);
    rct={Reaction.KEGG_ID}';
    handleWaitbar=waitbar(0,'Compiling data...');
    for k=1:size(RPAIR_new_table,2)
        waitbar(k/size(RPAIR_new_table,2),handleWaitbar,['Compiling data: ' num2str(k) ' of ' num2str(size(RPAIR_new_table,2))]);
        rc=RPAIR_new_table(k).REACTION';
        f=find(ismember(rct,rc));
        for j=1:size(f,1)
            Reaction(f(j)).pID{end+1,1}=RPAIR_new_table(k).ID;
            Reaction(f(j)).pID=unique(Reaction(f(j)).pID);
        end
    end
    delete(handleWaitbar)
end        
        
        %for j=1:size(RPAIR_new_table(k).REACTION,2)
        %    a=RPAIR_new_table(k).ADJACENCY_REACTION{1,j};
        %    if ~isempty(Reaction(a).pID)
        %        Reaction(a).pID{end+1,1}=RPAIR_new_table(k).ID;
        %    else
        %        Reaction(a).pID{end+1,1}=RPAIR_new_table(k).ID;
        %    end
        %    Reaction(a).pID=unique(Reaction(a).pID);%do you want this to be unique because there are 2 compounds for metabolites
        %end


%add the adjacency goes into the new RPAIR table(h)
%Metabolite entries add the RPAIR id for entry(h)
%Create an empty cells in REACTION and METABOLITES--> find (h) rcn and
%cpd use this to assign Metabolites(h).RPAIR(end+1,1)={RPAIR_new_table(k).ID};
%done for reaction list by adding the pID into the Reaction


%run the nc_SortDataReaction and SortDataMetabolites
%the run Pair_Reaction_compound and Pair_ID_compound