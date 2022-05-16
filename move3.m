function [ pos ] = move(x,y,di)
% returns next position based on current position and direction
    
if di == 1
    pos = [x-1 y];
elseif di == 2
    pos = [x-1 y+1];
elseif di == 3
    pos = [x y+1];
elseif di == 4
    pos = [x+1 y+1];
elseif di == 5
    pos = [x+1 y];
elseif di == 6
    pos = [x+1 y-1];
elseif di == 7
    pos = [x y-1];
elseif di == 8
    pos = [x-1 y-1];
elseif di == 0
    pos = [x y];
    
end

end

