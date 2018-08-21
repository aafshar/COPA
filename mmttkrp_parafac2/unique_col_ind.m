function [colind] = unique_col_ind(X)

[~, columns] = find(X);

unq = zeros(size(columns));
ctr = 0;
temp = -1;
for j=1: length(columns)
    if (temp ~= columns(j))
        ctr = ctr + 1;
        unq(ctr) = columns(j);
        temp = columns(j);
    end
end
colind = unq(1:ctr);