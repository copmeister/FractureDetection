function [ stepsize ] = stepcheck( area,dir,pre,prepre )
% checks proposed jump to see if any obstacles will be encountered.

% CAN'T BE AN EVEN AREA!!!

l = size(area,1); % length of sides
vert = area(:,(l+1)/2); % center vertical column
horz = area((l+1)/2,:); % center horizontal row
tbdiag = diag(area); % top(l) to bottom(r) diagonal

for i = 1:l-1
    avtbdiag(i) = sum(sum(area(i:i+1,i:i+1)));
end

btdiag = diag(flip(area)); % bottom(l) to top(r) diagonal

for i = 1:l-1
    avbtdiag(i) = sum(sum(area(l-i:l+1-i,i:i+1)));
end

if dir == 1
    if sum(vert(1:(l-1)/2)) == 0 
        % from bot to area edge in direction 1
        stepsize = pre+prepre;
    else
        stepsize = 1;
    end
end
  
if dir == 2
    if sum(btdiag(((l+1)/2)+1:l)) == 0 && any(avbtdiag(ceil(l/2):l-1)==2) == 0
        % from bot to area edge in direction 2 & check for gaps through 2
        % diagonal pixels
        stepsize = pre+prepre;
    else
        stepsize = 1;
    end
end

if dir == 3
    if sum(horz(((l+1)/2)+1:l)) == 0
        stepsize = pre+prepre;
    else
        stepsize = 1;
    end
end

if dir == 4
    if sum(tbdiag(((l+1)/2)+1:l)) == 0 && any(avtbdiag(ceil(l/2):l-1)==2) == 0
        stepsize = pre+prepre;
    else
        stepsize = 1;
    end
end

if dir == 5
    if sum(vert(((l+1)/2)+1:l)) == 0
        stepsize = pre+prepre;
    else
        stepsize = 1;
    end
end

if dir == 6
    if sum(btdiag(1:(l-1)/2)) == 0 && any(avbtdiag(1:floor(l/2))==2) == 0
        stepsize = pre+prepre;
    else
        stepsize = 1;
    end
end

if dir == 7
    if sum(horz(1:(l-1)/2)) == 0
        stepsize = pre+prepre;
    else
        stepsize = 1;
    end
end

if dir == 8
    if sum(tbdiag(1:(l-1)/2)) == 0 && any(avtbdiag(1:floor(l/2))==2) == 0
        stepsize = pre+prepre;
    else
        stepsize = 1;
    end
end

if dir == 0
    stepsize = 1;
end

