function [ direc ] = dircha( x,y,di,m,n )
% image boundary behaviour, bounces from edge, stuck in corner  

% hits top edge
if x==1 && di==2
    if rand>0.5
    direc = 4;
    else direc = 5;
    end
elseif x==1 && di==8
    if rand>0.5
    direc = 6;
    else direc = 5;
    end
elseif x==1 && di==1
    if rand>0.5
    direc = 4;
    else direc = 6;
    end
elseif x==1 && ismember(di,[3 4 5 6 7]) == 1
    direc = 0;
end
   
% hits bottom edge
if x==m && di==6
    if rand>0.5
    direc = 8;
    else direc = 1;
    end
elseif x==m && di==4
    if rand>0.5
    direc = 2;
    else direc = 1;
    end
elseif x==m && di==5
    if rand>0.5
    direc = 8;
    else direc = 2;
    end
elseif x==m && ismember(di,[1 2 3 7 8]) == 1
    direc = 0;
end

% hits left edge
if y==1 && di==8
    if rand>0.5
    direc = 2;
    else direc = 3;
    end
elseif y==1 && di==6
    if rand>0.5
    direc = 4;
    else direc = 3;
    end
elseif y==1 && di==7
    if rand>0.5
        direc = 2;
    else direc = 4;
    end
elseif y==1 && ismember(di,[1 2 3 4 5]) == 1
    direc = 0;
end

% hits right edge
if y==n && di==2
    if rand>0.5
    direc = 8;
    else direc = 7;
    end
elseif y==n && di==4
    if rand>0.5
    direc = 6;
    else direc = 7;
    end
elseif y==n && di==3
    if rand>0.5
        direc = 8;
    else direc = 6;
    end
elseif y==n && ismember(di,[1 5 6 7 8]) == 1
    direc = 0;
end

% direction unaltered when not at an edge
if x>1 && x<m && y>1 && y<n
    direc = di;
end

% movement stopped if top left corner reached
if x==1 && y==1
    direc = 0;
else
end

% movement stopped if top right corner reached
if x==1 && y==n
    direc = 0;
else
end

% movement stopped if bottom left corner reached
if x==m && y==1
    direc = 0;
else
end

% movement stopped if bottom right corner reached
if x==m && y==n
    direc = 0;
else 
end

% if no direction then confirm stopped
if di == 0
    direc = 0;
end
end

