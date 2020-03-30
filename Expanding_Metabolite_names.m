%Name
function [List_of_index]=Expanding_Metabolite_names(Metabolites)
    tic;
    handleWaitbar=waitbar(0,'Please wait...');
    Metabolite_expanded_name={};
    Metabolite_expanded_id={};
    Metabolite_expanded_index=[];
    List_of_index={};
    Mass_diff=5;
    %making index for when there are duplicates

    for i=1:size(Metabolites,2)
        waitbar(i/size(Metabolites,2),handleWaitbar,['Expanding metabolite names:' num2str(i) ' of ' num2str(size(Metabolites,2)) ' for ' num2str(toc) ' seconds']);
        size_of_Name=size(Metabolites(i).NAME,1);
        for j=1:size_of_Name
            Metabolite_expanded_name{end+1,1}=Metabolites(i).NAME{j,1}; %#ok
            Metabolite_expanded_id{end+1,1}=Metabolites(i).KEGG_ID; %#ok
            Metabolite_expanded_index(end+1,1)=i; %#ok
            
        end
    
    end
    
    close(handleWaitbar)
    handleWaitbar=waitbar(0,'Please wait...');
    for i=1:size(Metabolites,2)
        waitbar(i/size(Metabolites,2),handleWaitbar,['Compiling metabolite data:' num2str(i) ' of ' num2str(size(Metabolites,2)) ' for ' num2str(toc) 'seconds']);
        Same_as={};
        if ~isempty(Metabolites(i).NAME)
            if ~isempty(Metabolites(i).SAME_AS)
                Same_as=ismember(Metabolite_expanded_id,Metabolites(i).SAME_AS);%finding all the same_as from metabolites and adding it to a
            end 
            a=ismember(Metabolite_expanded_name,Metabolites(i).NAME);%a is finding index of which Metabolite names in Metabolite_expanded_names is same as Metabolites(i).NAME
                
            List_of_index{i,1}=unique(Metabolite_expanded_index(a));    %#ok position in Metabolites where there are duplicates
            for p=size(List_of_index{i,1},1):-1:2  
                index_first=List_of_index{i,1}(1);
                correct_mass=str2double(Metabolites(index_first).MASS);
                if abs((correct_mass)-(str2double(Metabolites(List_of_index{i,1}(p)).MASS)))>Mass_diff || isnan(abs((correct_mass)-(str2double(Metabolites(List_of_index{i,1}(p)).MASS))))
                    List_of_index{i,1}(p)=[];
                elseif isempty(Metabolites(List_of_index{i,1}(p)).MASS)%if mass in metabolite is empty(delete)unless
                    if isempty(Same_as)%Same_as is empty then delete--> else use same_as to add names
                         List_of_index{i,1}(p)=[]; 
                    end
                end
            end
            if size(List_of_index{i,1},1)==1 && isempty(Same_as) %if size of List_of_index{i,1}==1, it means that List_of_index{i,2}==1 too, determine by mass
                List_of_index{i,2}=unique(Metabolite_expanded_id(a));%#ok
                List_of_index{i,2}=List_of_index{i,2}(1);%#ok
                index_for_name=ismember(Metabolite_expanded_id(a),List_of_index{i,2}(1));%since we want only 1 entry from the first List_of_index{i,2}--> the indexing of the name must be the same
                List_of_index{i,3}=Metabolite_expanded_name(a); %#ok
                List_of_index{i,3}=unique(List_of_index{i,3}(index_for_name));%#ok
            elseif ~isempty(Same_as)
                List_of_index{i,1}=unique(cat(1,List_of_index{i,1},Metabolite_expanded_index(Same_as)));%#ok
                List_of_index{i,2}=unique(cat(1,Metabolite_expanded_id(a),Metabolite_expanded_id(Same_as)));%#ok
                List_of_index{i,3}=unique(cat(1,Metabolite_expanded_name(a),Metabolite_expanded_name(Same_as)));%#ok
            else
                List_of_index{i,2}=unique(Metabolite_expanded_id(a));%#ok
                List_of_index{i,3}=unique(Metabolite_expanded_name(a));%#ok
            end
        end
    end
    List_of_index_Mass={};
    for k=1:size(List_of_index,1)
        Sorted_list=sort(List_of_index{k,2});
        if size(List_of_index{k,2},1)>1
            List_of_index_Mass{end+1,1}=Sorted_list;%#ok
        end
    end
    close(handleWaitbar)
    disp(toc)
end
        
%replace both Metabolites and Reactions RPAIR with the C0000i
%1 using
%Metabolite_expanded_index(ismember(Metabolite_expanded_name,Metabolites(k).NAME))'

%find all the indexing for each C0000i

%add all the details(pathway,reaction,NAME) of the greater i to the
%smallest C0000i and remove the duplicates

%for REACTIONS use strfind to compare each of the RPAIR compounds to the
%C0000i which has duplicates. Replace the RPAIR with the updated compound
%names(sorted)


