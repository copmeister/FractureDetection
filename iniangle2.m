function [ direction ] = iniangle( space,x,y )
% sets initial angle (direction) for bot, which remains it prefered
% direction until boundary bounce

for i = 1:8
   for j = 1:8
      
       if atan2d(x-i,j-y) > 112.5 && atan2d(x-i,j-y) <= 157.5
           dir(i,j) = 8;  
       elseif atan2d(x-i,j-y) > 67.5 && atan2d(x-i,j-y) <= 112.5
           dir(i,j) = 1;  
       elseif atan2d(x-i,j-y) > 22.5 && atan2d(x-i,j-y) <= 67.5
           dir(i,j) = 2;  
       elseif atan2d(x-i,j-y) > -22.5 && atan2d(x-i,j-y) <= 22.5
           dir(i,j) = 3;  
       elseif atan2d(x-i,j-y) > -67.5 && atan2d(x-i,j-y) <= -22.5
           dir(i,j) = 4;  
       elseif atan2d(x-i,j-y) > -112.5 && atan2d(x-i,j-y) <= -67.5
           dir(i,j) = 5;  
       elseif atan2d(x-i,j-y) > -157.5 && atan2d(x-i,j-y) <= -112.5
           dir(i,j) = 6;  
       elseif atan2d(x-i,j-y) > 157.5 || atan2d(x-i,j-y) <= -157.5
           dir(i,j) = 7;  
       end
       
   end
end

dir2 = dir';
pos = space';
dir3(1,:) = dir2(:)';
dir3(2,:) = pos(:)';
dir4 = dir3';

[a,~,c] = unique(dir4(:,1));
out = [a, accumarray(c,dir4(:,2))];
out(:,3) = cumsum(out(:,2));

rn = rand;
it = 1;

while it < length(out)+1

if rn <= out(it,3)
    direction = out(it,1);
    it = length(out)+1;
else
end
    
it = it+1;

end 


