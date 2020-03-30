function [Reaction,Sorted_list]=Removing_duplicates_reaction(Reaction,List_of_index)
    tic;
    handleWaitbar=waitbar(0,'Initializing...');
    List_of_index_Reaction={};
    for k=1:size(List_of_index,1)
        Sorted_list=sort(List_of_index{k,2});
        if size(List_of_index{k,2},1)>1
            List_of_index_Reaction{end+1,1}=Sorted_list;%#ok
        end
    end
    
    for i=1:size(Reaction,2)
        waitbar(i/size(Reaction,2),handleWaitbar,['Removing duplicates: ' num2str(i) ' of ' num2str(size(Reaction,2))]);
        size_of_RCLASS=size(Reaction(i).RCLASS,1);
        for j=1:size_of_RCLASS
            c={};
            d={};
            Reaction_expanded_RCLASS=strtrim(split(Reaction(i).RCLASS{j,1},'_',2));
            for h=1:size(List_of_index_Reaction,1)

                a=List_of_index_Reaction{h,1}(ismember(List_of_index_Reaction{h,1},Reaction_expanded_RCLASS(1)));
                b=List_of_index_Reaction{h,1}(ismember(List_of_index_Reaction{h,1},Reaction_expanded_RCLASS(2)));
                if ~strcmp(a,List_of_index_Reaction{h,1}(1))
                    c=List_of_index_Reaction{h,1}(1);
                elseif isempty(c)
                    c=Reaction_expanded_RCLASS(1);
                end
                if ~strcmp(b,List_of_index_Reaction{h,1}(1))
                    d=List_of_index_Reaction{h,1}(1);
                elseif isempty(d)
                    d=Reaction_expanded_RCLASS(2);
                end
            end
            %if ~isempty(c) && ~isempty(d)
            if ~strcmp(Reaction(i).RCLASS{j,1},strcat(c,'_',d))
                Reaction(i).RCLASS{j,1}=strcat(c,'_',d);
            end
            %elseif isempty(c) && ~isempty(d) 
             %   Reaction(i).RCLASS{j,1}=strcat(Reaction_expanded_RCLASS(1),'_',d);
          %  elseif ~isempty(c) && isempty(d) 
            %    Reaction(i).RCLASS{j,1}=strcat(c,'_',Reaction_expanded_RCLASS(2));
           % end
        end
    end
    disp(toc)
    close(handleWaitbar)
end 

%using List_of_index, create another list where List_of_index{:,2} size>=2
%with all the names and id, then use that for ismember compare to
%Reaction_exapneded_RCLASS(dont forget to sort)
%for a=1:size(List_of_index_Reaction,1)
 %               for b=size(List_of_index_Reaction,2):-1:2
  %                  if strcmp(List_of_index_Reaction(a,b),Reaction_expanded_RCLASS(1))
   %                     c=List_of_index_Reaction(a,1);
    %                end
     %               istrcmp(List_of_index_Reaction(a,b),Reaction_expanded_RCLASS(2))
      %                  d=List_of_index_Reaction(a,1);
       %             end
%                end
 %           end
  %          if ~isempty(c)
   %             if ~isempty(d)
    %                Reaction(i).RCLASS{j,1}=strcat(c(1),'_',d(1));
     %           else
      %              Reaction(i).RCLASS{j,1}=strcat(c(1),'_',Reaction_expanded_RCLASS(2));
       %         end
        %    elseif ~isempty(d)
         %       Reaction(i).RCLASS{j,1}=strcat(Reaction_expanded_RCLASS(1),'_',d(1));
          %  else
           %     Reaction(i).RCLASS{j,1}=strcat(Reaction_expanded_RCLASS(1),'_',Reaction_expanded_RCLASS(2));
            %end


 