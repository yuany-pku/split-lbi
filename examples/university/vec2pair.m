function [left,right] = vec2pair(n)
    x = triu(ones(n),1);
    [i,j] = find(x);
    left = i;
    right = j;
end