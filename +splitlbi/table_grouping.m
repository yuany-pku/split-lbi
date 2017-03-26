function [] = table_grouping(score,name)
% DESCRIPTION
% Return Grouping result and save it as csv
p = length(score);
[score_des,~] = sort(score,'descend');
for i = 2:p
    if abs(score_des(i) - score_des(i-1)) < 1e-4
        score_des(i) = score_des(i-1);
    end
end
score_uniq = sort(unique(score_des),'descend');
group_index = zeros(p,1);
for i=1:p
    group_index(i) = find(abs(score_uniq - score(i)) < 1e-4);
end
T = table(name,group_index,'VariableNames',{'Candidates' 'Group'});
writetable(T,'Grouping_res.csv','Delimiter',',');
end