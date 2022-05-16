close all;
clear;

% =========================================================================
% --- IMPORT IMAGE --------------------------------------------------------
% =========================================================================

s1 = tic;

A=imread('Picture1.jpg'); % read source image
a2 = imgaussfilt(A, 0.5); % 0.5 low, 1.0 med, 2.0 high blur
B=rgb2gray(a2); % convert from colour to greyscale
C=edge(B,'canny',0.045); % canny edge detection < more detail > less detail
[m,n] = size(C);
thres = 0.20; % 0 = no effect, < less effect, > more effect

Gx = zeros(m,n); % prelocate horizontal grid 
Gy = zeros(m,n); % prelocate vertical grid
direc = zeros(m,n); % prelocate gradient grid - direction
D = zeros(m,n); % prelocate edge intesity grid

for i=2:m-1
    for j=2:n-1
        %Sobel mask for x-direction:
        Gx(i,j)=((2*C(i+1,j)+C(i+1,j-1)+C(i+1,j+1))-(2*C(i-1,j)+C(i-1,j-1)+C(i-1,j+1)));
        %Sobel mask for y-direction:
        Gy(i,j)=((2*C(i,j+1)+C(i+1,j+1)+C(i-1,j+1))-(2*C(i,j-1)+C(i+1,j-1)+C(i-1,j-1)));
      
        %The gradient of the image
        %B(i,j)=abs(Gx)+abs(Gy);
        D(i,j)=sqrt(Gx(i,j).^2+Gy(i,j).^2); % edge intensity
        direc(i,j) = atand(Gy(i,j)/Gx(i,j)); % edge direction
        %D(i,j)=sqrt(Gx.^2+Gy.^2);
      
    end
end

ms = 11; % mask size for variance of pixel gradient
si = (ms-1)/2; % radius of mask
gravar = zeros(m,n); % prelocate gradient variance map
direc(isnan(direc))=0; % convert any NaN values to zero


for i = (ms+1)/2:m-((ms-1)/2)
   for j = (ms+1)/2:n-((ms-1)/2)
   
       gravar(i,j) = var(reshape(direc(i-si:i+si,j-si:j+si),[],1));
       % build gradient variance map, pixel by pixel
       
   end
end

z = gravar/max(max(gravar)); % standardize range to 0:1

z(z<thres)=0; % replace any values under threshold to zero

E = C.*z;   % process to remove any pixels from canny edge detection
            % eliminated during gradient variance
zz = E>0;   % any pixel with values >0 remaining transfered to
            % new map with logical value 1
zzz = zz*255; % turn map into 0's and 255's
zzzz = uint8(zzz); % convert into greyscale image

pre_process = toc(s1)

% =========================================================================
% -------------------------------------------------------------------------
% =========================================================================

s2 = tic;

it = 1; % set start at first iteration
%I1 = imread('pp2.jpg'); 
%I2 = rgb2gray(I1);
I2 = B; % greyscale of source image
%I3 = double(I2);
I3 = zzz; % 0 and 255 map of processed images
[m,n] = size(I3);
ms = round(m/8);
ns = round(n/8);
martyr = [0 0]; % intial death to prevent crash when searching it-1 later


% =========================================================================
% --- INITIAL CONDITIONS --------------------------------------------------
% =========================================================================

nb = 25; % number of bots
maxbp = 150; % max bot population
maxb = 500; % max number of bots allowed to exist
minb = 1; % min number of bots
iterations = 75; % number of iterations
spawnrate = 0.3; % spawn rate (< high chance, > low chance)
mo = 20; % movement to monitor
th = 15; % threshold for average movement (< forgiving > unforgiving) d:15
hm = 7; % diameter for heat detection
exht = 0.8; % extreme heat causes death
ll = 20; % area for final martyr (distance from martyr to edge of zone) 

% =========================================================================
% -------------------------------------------------------------------------
% =========================================================================

activity = zeros(8); % divide image into 64 regions

for i = 1:7
   for j = 1:7
   activity(i,j) = var(reshape(I3(i*ms-ms+1:i*ms,j*ns-ns+1:j*ns),[],1));
   activity(8,j) = var(reshape(I3(8*ms-ms+1:m,j*ns-ns+1:j*ns),[],1));
   activity(i,8) = var(reshape(I3(i*ms-ms+1:i*ms,8*ns-ns+1:n),[],1));
   end
end

prob = activity./sum(sum(activity)); % activity of each region 0 to 1


I3(I3<=100)=0; I3(I3>=101)=255; % final clean of image
I3=I3/255; % turn full image to 0 or 1


% =========================================================================
% -- HEAT MAP -------------------------------------------------------------
% =========================================================================

rad = (hm-1)/2; % radius of heat detection
hemap = zeros(m,n); % prelocate heatmap

ang = 0:1:360; % specify 360 degree heat analysis
pos = round([rad*cosd(ang);rad*sind(ang)])'; % position of each degree
        % in relation to centre accoring to heat radius, rounded to give
        % specific vector location
un = unique(pos,'rows'); % eliminate duplicate values
val = length(un); % number of vector locations on circumfrence of heat 
        % detection according to radius
        
for i = 1+rad:m-rad
   for j = 1+rad:n-rad
      
       for k = 1:val
          
           un(k,3) = I3(un(k,1)+i,un(k,2)+j);
           % determine if heat (pixel) is located at each point on
           % the circumfrence of heat analysis circle
           
       end
 
       hemap(i,j) = sum(un(:,3)).^2;
       % sum all heat detected on the circumfrence and recorded as
       % temperature at each pixel point   
       
   end    
end

% show heat map with temperature bar
colormap(hot); 
imagesc(hemap);
colorbar;
title('Heat Map')

maxt = max(reshape(hemap,[],1)); % maximum temperature on heatmap
mint = min(reshape(hemap,[],1)); % minimum temperature on heatmap
baset = mint+round((maxt-mint)*0.1); % set base temperature accorinding to heat 
                                % range (default is 10%)
                                
extr = maxt*exht; % extreme heat causes death <<<<< trial measure >>>>>

% =========================================================================
% -------------------------------------------------------------------------
% =========================================================================

bl = zeros(iterations,nb*10); % prelocate the botloc results

for i = 1:nb
    %bl(it,1:2) = [16 18]; % first bot
    %bl(it,1:3) = [500 400 1]; % first bot with direction
    %bl(it,11:12) = [5 5]; % second bot
    %bl(it,11:13) = [27 3 2]; % second bot with direction
    
    bl(it,(i*10)-9:(i*10)-8) = [randi(m-2)+1 randi(n-2)+1]; % starting location
    
    %bl(it,(i*10)-7) = randi(8); % prefered direction %%FINISHED
    
    bl(it,(i*10)-7) = iniangle2(prob,round(bl(it,(i*10)-9)/m*8),round(bl(it,(i*10)-8)/n*8));
    % initial angle from iniangle2.m function, each bot aims for active
    % region
    
    bl(it,(i*10)-6) = bl(it,(i*10)-7); % compromised direction
    bl(it,(i*10)-5) = 0; % 0 alive or 1 dead
    bl(it,(i*10)-4) = 1; % step size, initial 1
    bl(it,(i*10)-3) = 0; % distance from mo'th previous location
    bl(it,(i*10)-2) = baset; % bot heat normal value
    bl(it,(i*10)-1) = 0;
    bl(it,(i*10)) = 0;
end

figure
imshow(zzz) % show image

for i = 1:it*nb
    history(i,:) = bl(1,(i*10)-9:(i*10)-7); % enter initial bot details in history
end

heat_map = toc(s2)

s3 = tic;

%%
% =========================================================================
% =========================================================================
% ------ START MAIN BOT LOOP ----------------------------------------------
% =========================================================================
% =========================================================================

while it < iterations
    
for i = 1:nb % BOT LOOP!!

if bl(it+1,(i*10)-5) == 0 % CHECK IF ALIVE FOR THIS ITERATION
   
if sum(bl(it,(i*10)-9:(i*10)-7)) > 0 % EXTRA CHECK IF ALIVE LAST ITERATION

    % calculate distance from previous movement indicator
    if it > mo
        bl(it,(i*10)-3) = sqrt(sum((bl(it,(i*10)-9:(i*10)-8)-bl(it-mo,(i*10)-9:(i*10)-8)).^2));
        % calculate distance from mo'th previous location and record in bot
        % history
    else
    end
    
% CHECK DIRECTION (RETURN AT EDGE OF IMAGE - CORNER CAUSES HALT)
  
    bl(it+1,(i*10)-7) = dircha(bl(it,(i*10)-9),bl(it,(i*10)-8),bl(it,(i*10)-7),m,n);      
    %Sets preferred direction for next move
    
    bl(it+1,(i*10)-6) = bl(it+1,(i*10)-7);
    %Sets compromised default to match preferred
    
    if bl(it+1,(i*10)-7) == 0
    % if preferred direction is zero
    bl(it+1,(i*10)-5) = 1;
    % kill the bot
    else
    end
    

%% CHECK FOR OBSTACLE
    
if bl(it,(i*10)-9) >= 2 && bl(it,(i*10)-9) <= m-1 && bl(it,(i*10)-8) >= 2 && bl(it,(i*10)-8) <= n-1
%Check that bot is not on an edge
    
    area = I3(bl(it,(i*10)-9)-1:bl(it,(i*10)-9)+1,bl(it,(i*10)-8)-1:bl(it,(i*10)-8)+1);
    % capture the 3x3 area around the bot
    
    cand = obst(area);
    % order the directional candidates in order of preference 
    
    prefapprov = 0;
    % preferred direction requires approval
    
    prefloc = move3(bl(it,(i*10)-9),bl(it,(i*10)-8),bl(it,(i*10)-7));
    % find next location if move in preferred direction
    prefcase = [prefloc bl(it,(i*10)-7)];
    % position and direction if moved in preferred direction
    
    if sum(ismember(history,prefcase,'rows')) == 0
           prefapprov = 1;
           % if the position and direction haven't occured, approve the
           % move in principle
    else
    end

    
    if any(cand==bl(it+1,(i*10)-7)) == 0 || prefapprov == 0
    % is the prefered missing from the directional candidate list
        
        %ALL THIS SECTION SKIPPED IF PREFERED DIRECTION IS ALLOWED
        
        pref = 1;
        approv = 0;
        
        
        while approv == 0 && bl(it+1,(i*10)-5) == 0
        
        
        bl(it+1,(i*10)-6) = cand(pref); %Make compromised direction the prefered direction
        
        % if botloc(it+1,(i*10)-6) == 0 then kill bot
            
              
        if bl(it+1,(i*10)-7) == bl(it+1,(i*10)-6)
            tempdir = bl(it+1,(i*10)-7);
        else
            tempdir = bl(it+1,(i*10)-6);
        end
       
        % build test case, using next location based on calculated
        % direction
        testloc = move3(bl(it,(i*10)-9),bl(it,(i*10)-8),tempdir);
        
        % test case to include direction for the next location
        testcase = [testloc tempdir];
        
        % check if testcase has already occured
        if sum(ismember(history,testcase,'rows')) == 0
           approv = 1;
        else
           pref = pref+1;
        end
        
        if testcase(:,3) == 0 || pref >= 9
            bl(it+1,(i*10)-5) = 1; % MAYBE SPAWN ELSEWHERE??????? %
        else
        end
        
        end
        
            
    else
    end
    
else
  
end

%% Step size

    if bl(it+1,(i*10)-6) ~= bl(it,(i*10)-6)
        bl(it+1,(i*10)-4) = 1;
        % if next direction different to previous, step size is 1
    else
    
        if it < 2
            bl(it+1,(i*10)-4) = 1;
            % if less that 2 iterations occured, step size is 1
        else
        
            if bl(it,(i*10)-4) == 1 && bl(it-1,(i*10)-4) ~= 1
                bl(it+1,(i*10)-4) = 1;
                % if current step size is 1 AND previous step size does NOT
                % equal 1, then next step size must equal one (to ensure
                % sequence always begins with 1,1,...)
            else
        
                bl(it+1,(i*10)-4) = bl(it,(i*10)-4) + bl(it-1,(i*10)-4);
                % next step size equals sum of previous 2 step sizes
                % (Fibonacci) IF NO OTHER PREVIOUS CONDITIONS APPLY
        
                    if bl(it,(i*10)-9)-bl(it+1,(i*10)-4) < 1 || bl(it,(i*10)-9)+bl(it+1,(i*10)-4) > m || bl(it,(i*10)-8)-bl(it+1,(i*10)-4) < 1 || bl(it,(i*10)-8)+bl(it+1,(i*10)-4) > n
                        bl(it+1,(i*10)-4) = 1;
                        % if new step size causes bot to move beyone edge
                        % of image, reset step size to 1
                    else
    
                steparea = I3(bl(it,(i*10)-9)-bl(it+1,(i*10)-4):bl(it,(i*10)-9)+bl(it+1,(i*10)-4),bl(it,(i*10)-8)-bl(it+1,(i*10)-4):bl(it,(i*10)-8)+bl(it+1,(i*10)-4));
                % build area to check for edges based on suggested step
                % size
            
                bl(it+1,(i*10)-4) = stepcheck(steparea,bl(it+1,(i*10)-6),bl(it,(i*10)-4),bl(it-1,(i*10)-4));
                % retrieve step size from stepcheck function (1 if edge
                % encountered over step size, unchanged if no obstacle
                % (edge) found

                    end   
    
            end
    
        end

    end

%% MOVE BOT    

    
    if bl(it+1,(i*10)-7) == bl(it+1,(i*10)-6)
    
    bl(it+1,(i*10)-9:(i*10)-8) = move4(bl(it,(i*10)-9),bl(it,(i*10)-8),bl(it+1,(i*10)-7),bl(it+1,(i*10)-4));
    % if compromised direction is SAME as preferred, move accoring to
    % preferred and step size
    
    history(end+1,:) = [bl(it+1,(i*10)-9:(i*10)-8),bl(it+1,(i*10)-7)];
    % update history with new position and direction (preferred)
    
    else
    
    bl(it+1,(i*10)-9:(i*10)-8) = move4(bl(it,(i*10)-9),bl(it,(i*10)-8),bl(it+1,(i*10)-6),bl(it+1,(i*10)-4));
    % if compromised direction is DIFFERENT to preferred, move according to
    % compromised and step size
    
    history(end+1,:) = [bl(it+1,(i*10)-9:(i*10)-8),bl(it+1,(i*10)-6)];
    % update history with new position and direction (compromised)
        
    end

    
%% HEAT CHANGE

bl(it+1,(i*10)-2) = bl(it,(i*10)-2)+hemap(bl(it+1,(i*10)-9),bl(it+1,(i*10)-8))-baset;
% new heat value is equal to: current bot heat PLUS heat value of next
% position MINUS base temperature, e.g. current = 30
                                         % next = 15 (+)
                                         % base = 20 (-)
                                     % new heat = 25 (=)  

if bl(it+1,(i*10)-2) < baset
    bl(it+1,(i*10)-2) = baset;
    % if bot temperature drops below base, then set to base temperature
    % (bots cannot drop below base temperature)
else
end

if bl(it+1,(i*10)-2) > extr
    if bl(it+1,(i*10)-9)-ll >= 1 && bl(it+1,(i*10)-9)+ll <= m && bl(it+1,(i*10)-8)-ll >= 1 && bl(it+1,(i*10)-8)+ll <= n
        martyr(end+1,:) = bl(it+1,(i*10)-9:(i*10)-8);
        bl(it+1,(i*10)-5) = 1;
        bl(it+1,(i*10)-1) = 1;
        % if bot overheats according to the maximum temperature for survival
        % set by 'extr' variable, add to death list (potential fracture
        % location)
    else
    end
else
end
    
         
%% CHECK IF SPAWN NEEDED FOR PREVIOUS MOVE

if bl(it,(i*10)-9) > 2 && bl(it,(i*10)-9) < m-2 && bl(it,(i*10)-8) > 2 && bl(it,(i*10)-8) < n-2
    % can only check for spawn if bot is away from edge of image
    area = I3(bl(it,(i*10)-9)-2:bl(it,(i*10)-9)+2,bl(it,(i*10)-8)-2:bl(it,(i*10)-8)+2);
    % record 5x5 matrix about bot to check for spawn conditions
    
    needbot(:) = spawner(area,bl(it,(i*10)-9),bl(it,(i*10)-8),bl(it,(i*10)-6));
    % spawner function will decide if a spawn is recommended or not
    
    if needbot(1) == 1 && it > 2 && alive(it-1) < maxbp && nb < maxb
        % IF spawn is recommended
        % AND more that 2 iterations have passed
        % AND the bot population from the previous iteration was below the
            % maximum allowed bot population specified by maxbp variable
        % AND the number of bot to have ever lived is less that the total
        % allowed bots allowed to have ever lived specified by the maxb
            % variable
        
        % THEN respawn is decided by chance
        
        randspawn = rand; % random number to check spawn probability against
        
        if randspawn > spawnrate
            % if random number is greater than spawnrate variable then new
            % bot will be spawned
            
            nb = nb+1; % increase number of bots to have ever lived by 1
            bl(:,nb*10-9:nb*10) = zeros;
            % populate bl matrix with 10 new columns initially with zeros            
            bl(it+1,(nb*10)-9:(nb*10)-7) = needbot(2:4);
            % initial bot conditions set in correct iteration row according
            % to initial bot rules, here [x y] location and preferred
            % direction
            bl(it+1,(nb*10)-6) = bl(it+1,(nb*10)-7);
            % set compromised direction equal to preferred direction
            bl(it+1,(nb*10)-4) = 1;
            % set initial step size to 1
        else
        end
        
    else
    end
    
else

end

else
end % EXTRA CHECK IF ALIVE LAST ITERATION

else
end % CHECK IF ALIVE FOR THIS ITERATION

% kill bot if below threshold movement average over movement period

if it > 2*mo
    % only check if the number of iterations completed is more than twice
    % the number of moves to monitor specified by the 'mo' variable
             
    if sum(bl(it+1,(i*10)-9:(i*10)-6))>0
        % if bot had no activity in previous iteration, bot was dead so
        % do not check movement
                
        if mean(bl(it-mo+2:it+1,(i*10)-3)) < th                     
            bl(it+1,(i*10)-5) = 1;
            % if the mean value over the previous 'mo' iterations has fallen
            % below the threshold specified by 'th' then kill the bot
            
        else
        end
    else
    end
    
end

end % BOT LOOP!!

alive(it) = nb-sum(bl(it,5:10:nb*10));
% number of bots minus any bot bots that died in this iteration

it = it+1;
% set new iteration




for i = 1:nb

if any(bl(:,(i*10)-5)) == 1
    bl(it+1,(i*10)-5) = 1;
    % ensure any bots killed this iteration are marked dead for new
    % iteration
else
end

end

if it == 500000 %***
    %bl(it+1,(i*10)-5) = 1;
    bl(it+1,5:10:nb*10) = 1;
    % if iteration exceed this limit kill all bots
else
end %***


hold on;

plot(bl(it-1,2),bl(it-1,1),'g.','MarkerSize',20)

end

bot_loop = toc(s3)

% =========================================================================
% =========================================================================
% ------ END MAIN BOT LOOP ------------------------------------------------
% =========================================================================
% =========================================================================

martyr(all(~martyr,2),:)=[]; % remove none existent martyrs
numart = length(martyr); % number of martys at end
centre = [round(m/2) round(n/2)]; % centre point of image



for i = 1:numart
    
   %if martyr(i,1)-ll >= 1 || martyr(i,1)+ll <= m || martyr(i,2)-ll >= 1 && martyr(i,2)+ll <= n
    
   martyr(i,3) = sqrt(sum((martyr(i,1:2)-centre).^2));
   % distance from centre of image to martyr
   
   martyr(i,4) = var(reshape(I3(martyr(i,1)-ll:martyr(i,1)+ll,martyr(i,2)-ll:martyr(i,2)+ll),[],1));
   % amount of activity in zone around matryr
   
   tdirec = reshape(direc(martyr(i,1)-ll:martyr(i,1)+ll,martyr(i,2)-ll:martyr(i,2)+ll),[],1);
   % list of all gradients contained in zone
   tdirec(all(~tdirec,2),:)=[];
   % remove all zero gradients
   martyr(i,5) = var(tdirec);
   % amount of variance in gradient across the zone around martyr
      
end

dis = 1;
act = 2;
gra = 2;

martyr(:,6)=martyr(:,3)-min(martyr(:,3));
martyr(:,7)=(martyr(:,6)/max(martyr(:,6)));
% alter distance to centre to a 0-1 rating

martyr(:,8)=martyr(:,4)-min(martyr(:,4));
martyr(:,9)=(martyr(:,8)/max(martyr(:,8)));
martyr(:,13)=2*abs(0.5-martyr(:,9)); % place 13,14,15 to work favour  
martyr(:,14)=martyr(:,13)-min(martyr(:,13)); % the bot suurounded by not
martyr(:,15)=(martyr(:,14)/max(martyr(:,14))); % least or most activity
% alter activity score to a 0-1 rating

martyr(:,10)=martyr(:,5)-min(martyr(:,5));
martyr(:,11)=(martyr(:,10)/max(martyr(:,10)));
% alter gradient variance to a 0-1 rating

martyr(:,12) = dis*martyr(:,7)+act*martyr(:,9)+gra*martyr(:,11);
martyr(:,16) = dis*martyr(:,7)+act*martyr(:,15)+gra*martyr(:,11);
% total score for each martyr for distance, activity and gradient variance

morder = sortrows(martyr,12);
morder2 = sortrows(martyr,16);
% reorder list according to total score



%%
% =========================================================================
% ------ PLOT BOT MOVEMENTS -----------------------------------------------
% =========================================================================

hold on;

plot(history(:,2),history(:,1),'r.')

hold on
martyr(all(~martyr,2),:)=[];

figure
plot(alive)
title('Bot population')
xlabel('Iterations')
ylabel('Bot population')

figure
plot(bl(:,8:10:nb*10))
hold on
plot([1,it],[baset baset], 'g--', 'linewidth',2)
text(1,baset+10,'Base temperature')
plot([1,it],[extr extr], 'r--', 'linewidth',2)
text(1,extr+10,'Extreme heat')
title('Bot temperature records')
xlabel('Iterations')
ylabel('Individual bot temperature')

figure
imshow(A)
hold on
plot(martyr(:,2),martyr(:,1),'b.','MarkerSize',20)
%plot(morder(1,2),morder(1,1),'r.','MarkerSize',30)
plot(morder2(1,2),morder2(1,1),'y.','MarkerSize',30)