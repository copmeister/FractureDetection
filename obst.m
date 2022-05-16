function [ dir ] = obst( v )
% This function returns the order of preference for directions for the
% centre point of any given 3x3 matrix

pos=(1:8); % set 8 possible directions

% prevent movement through corners by sealing them shut

% example   010
%           100
%           000
% cannot pass through top left corner

if sum([v(2,1),v(1,2)]) > 1
    pos(pos==8)=[];
end

if sum([v(2,1),v(3,2)]) > 1
    pos(pos==6)=[];
end

if sum([v(3,2),v(2,3)]) > 1
    pos(pos==4)=[];
end

if sum([v(2,3),v(1,2)]) > 1
    pos(pos==2)=[];
end

% prevent movement into obstacle

if v(1,1) == 1
    pos(pos==8)=[];
end

if v(1,2) == 1
    pos(pos==1)=[];
end

if v(1,3) == 1
    pos(pos==2)=[];
end

if v(2,1) == 1
    pos(pos==7)=[];
end

if v(2,3) == 1
    pos(pos==3)=[];
end

if v(3,1) == 1
    pos(pos==6)=[];
end

if v(3,2) == 1
    pos(pos==5)=[];
end

if v(3,3) == 1
    pos(pos==4)=[];
end

% best option

% this ensures that best options for direction are those closest to any
% detected edges without passing through obstacle

op(1) = mean([v(1,1) v(1,2) v(1,3)]);
op(2) = mean([v(1,2) v(1,3) v(2,3)]);
op(3) = mean([v(1,3) v(2,3) v(3,3)]);
op(4) = mean([v(2,3) v(3,3) v(3,2)]);
op(5) = mean([v(3,1) v(3,2) v(3,3)]);
op(6) = mean([v(2,1) v(3,1) v(3,2)]);
op(7) = mean([v(1,1) v(2,1) v(3,1)]);
op(8) = mean([v(2,1) v(1,1) v(1,2)]);

lv = length(pos);

for i = 1:lv
    pos(2,i) = op(pos(1,i))+rand/100; % add small random value to prevent draws
end

% order possibe directions highest mean to lowest
pos = sortrows(pos',-2)';

dir = pos(1,:);
dir(end+1:8) = 0;

end

