function [ spawn ] = spawner( area,x,y,di )

% check by direction whether there is a edge directly inline AND a space
% direcly inline after the edge ([1 0] condition), if this is found that a
% spawn is suggested inline with the parents position and direction in the
% space (0) directly beyone the edge (1), the parents direction is also
% inherited

% if [1 0] is not met, no spawn is returned

if di == 1 
    if sum([area(2,3),area(1,3)]==[1 0])==2
        spawn = [1 x-2 y di];
    else
        spawn = [0 x y di];
    end
end

if di == 2 
    if sum([area(2,4),area(1,5)]==[1 0])==2
        spawn = [1 x-2 y+2 di];
    else
        spawn = [0 x y di];
    end
end

if di == 3 
    if sum([area(3,4),area(3,5)]==[1 0])==2
        spawn = [1 x y+2 di];
    else
        spawn = [0 x y di];
    end
end

if di == 4
    if sum([area(4,4),area(5,5)]==[1 0])==2
        spawn = [1 x+2 y+2 di];
    else
        spawn = [0 x y di];
    end
end

if di == 5
    if sum([area(4,3),area(5,3)]==[1 0])==2
        spawn = [1 x+2 y di];
    else
        spawn = [0 x y di];
    end
end

if di == 6
    if sum([area(4,2),area(5,1)]==[1 0])==2
        spawn = [1 x+2 y-2 di];
    else
        spawn = [0 x y di];
    end
end

if di == 7
    if sum([area(3,2),area(3,1)]==[1 0])==2
        spawn = [1 x y-2 di];
    else
        spawn = [0 x y di];
    end
end

if di == 8
    if sum([area(2,2),area(1,1)]==[1 0])==2
        spawn = [1 x-2 y-2 di];
    else
        spawn = [0 x y di];
    end
end

