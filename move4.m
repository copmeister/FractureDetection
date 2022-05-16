function [ pos ] = move(x,y,di,ss)
% returns next position based on current position and direction
% factoring in step size (ss)
    
if di == 1
    pos = [x-ss y];
elseif di == 2
    pos = [x-ss y+ss];
elseif di == 3
    pos = [x y+ss];
elseif di == 4
    pos = [x+ss y+ss];
elseif di == 5
    pos = [x+ss y];
elseif di == 6
    pos = [x+ss y-ss];
elseif di == 7
    pos = [x y-ss];
elseif di == 8
    pos = [x-ss y-ss];
elseif di == 0
    pos = [x y];
    
end

end

