array = {NaN; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'FALSE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'FALSE'; NaN; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'FALSE'; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'FALSE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'FALSE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'FALSE'; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; 'TRUE'; NaN; NaN; 'FALSE'; NaN; 'TRUE'; NaN; NaN; NaN; 'FALSE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'FALSE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; 'TRUE'; NaN; NaN; NaN; NaN; NaN; 'TRUE'}

count_true = 0;
count_false = 0;

sum(strcmp(array(:),'TRUE'))

for i = 1:size(array,1)
    if strcmp(array{i},'TRUE')
        count_true = count_true + 1;
    elseif strcmp(array{i},'FALSE')
        count_false = count_false + 1;
    end
end

total_count = count_true + count_false;