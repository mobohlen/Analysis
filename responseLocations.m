clear all; close all;


%% Data files to analyze
 % specifiy dates as Mmm dd yyyy'
startDate = 'Feb 1 2016';
stopDate = 'Mar 28 2016';
dataFolder = 'Arya';
%dataFolder = 'Lana';

%% Load all data files
filePath = ['../Data/',dataFolder];
fileDir = dir(fullfile(filePath,'*.mat'));
d = 0;
for i =  1:length(fileDir)
    tempFile = fileDir(i).name;
    tempRead = load(fullfile(filePath,tempFile));
    if(datenum(tempRead.header.date)>=datenum(startDate) &&... 
        datenum(tempRead.header.date)<=datenum(stopDate))
        d = d+1;
        files{d} = fileDir(i).name;
        dateSort{d} = tempRead.header.date; %location of date in file name
    end
end
[~,sortIndex] = sort(datenum(dateSort)); 
files = files(sortIndex); %reorder summary files by date


 % cue dimension properties 
% col 1: Arya: blue = 1, red = 0 / Lana: horiz = 0, verti = 1 
% col 2: Arya: circle = 1, square = 0 / Lana: white = 0, black = 1
% col 3: Arya: white = 1, black = 0 / Lana: circle = 0, square = 1
% col 4: Arya: horiz = 1, verti = 0 / Lana: blue = 0, red = 1
stimprop(1,:) = [1 1 1 1];
stimprop(2,:) = [1 1 1 0];
stimprop(3,:) = [1 1 0 1];
stimprop(4,:) = [1 1 0 0];
stimprop(5,:) = [1 0 1 1];
stimprop(6,:) = [1 0 1 0];
stimprop(7,:) = [1 0 0 1];
stimprop(8,:) = [1 0 0 0];
stimprop(9,:) = [0 1 1 1];
stimprop(10,:) = [0 1 1 0];
stimprop(11,:) = [0 1 0 1];
stimprop(12,:) = [0 1 0 0];
stimprop(13,:) = [0 0 1 1];
stimprop(14,:) = [0 0 1 0];
stimprop(15,:) = [0 0 0 1];
stimprop(16,:) = [0 0 0 0];

stimWeights = [0.9 0.8 0.7 0.6];

weightedStims(:,1) = stimWeights(1).*stimprop(:,1) + (1-stimWeights(1)).*(1-stimprop(:,1));
weightedStims(:,2) = stimWeights(2).*stimprop(:,2) + (1-stimWeights(2)).*(1-stimprop(:,2));
weightedStims(:,3) = stimWeights(3).*stimprop(:,3) + (1-stimWeights(3)).*(1-stimprop(:,3));
weightedStims(:,4) = stimWeights(4).*stimprop(:,4) + (1-stimWeights(4)).*(1-stimprop(:,4));

N = length(files);
n=0;
% stimulus positions 1-4; 1 = leftUp, 2 = leftDown, 3 = rightUp, 4 = rightDown

for f =1:N
%% Read each file and index according to data
    currentFile = load(fullfile(filePath,files{f}));  
    h = currentFile.header;
    s = currentFile.summary;
   
    indResponse = find(s.respLocation ~= -1 &...
                       s.respLocation ~= 0)'; % indeces of all trials where touchscreen was touched
    indMiss = find(s.respLocation == -2)'; % indeces of trials where touchscreen was hit, stimuli were missed
    indValid= find(s.respLocation ~= -2 &...
                   s.respLocation ~= -1 &...
                   s.respLocation ~= 0)'; % indeces of valid trials, where either cue was hit
    indCorrect = find(s.correct == 1);
    %indTime = find(s.responseTime<=h.trialTime);
    indTime = find(s.responseTime<=(h.trialTime+h.gracePeriod));

    indValid(~ismember(indValid,indTime)) = [];
    indIgnore = find(s.correctLoc == -3);
    indValid(ismember(indValid,indIgnore)) = [];
    
%% Total counts of various trial types
    nResponse = length(indResponse);
    nValid = length(indValid);
    nMiss = length(indMiss);
    nCorrect = length(indCorrect);
    
%% Ignore days with not enough valid responses
if(nValid<30)
    disp(['Warning: n<30, ignoring ' files{f}]);
else
    n=n+1;
    dateI = regexprep(h.date,'_',' ');
    allDates{n} = dateI;
    titleStr = [h.monkey '  ' dateI, sprintf('  (T = %.1f s)',h.trialTime)];
%% Extract data from each specified trial
    if(strcmp(h.monkey,'Arya'))
        cueProbInfo = load('cueProbInfo.mat');
    elseif(strcmp(h.monkey,'Lana'))
        cueProbInfo = load('cueProbInfo_Lana.mat');
    else
        warning('Missing: "cueProbInfo.mat" for specified subject name');
    end
    ind = indValid;
    cueA = s.cueA(ind); cueB = s.cueB(ind);
    respStimulus = s.respStimulus(ind); 
    correct = s.correct(ind);
    responseTime = s.responseTime(ind);
    location = s.respLocation(ind);
    trialTime = s.trialTime(ind);
    locs(n,:) = [sum(location==1),sum(location==2),sum(location==3),sum(location==4)];
    locsprop(n,:) = locs(n,:)/length(ind);
    loc1vsloc4(n,:) = [sum((correct==1)&(location==1)),sum((correct==1)&(location==4))];

    trialTimeAll(n) = trialTime(1);
    
    Aloc = s.cueALoc(ind);
    Bloc = s.cueBLoc(ind);

    cueWDiffs = sum(weightedStims(cueA,:),2)-sum(weightedStims(cueB,:),2);
    cueWDiffs = cueWDiffs(cueWDiffs ~= 0);
    idxHigh = find(abs(cueWDiffs)>=1);
    idxLow = find(abs(cueWDiffs)<1);
    
    pairs = [1 2;1 3;1 4;2 3;2 4;3 4];
    for i = 1:6
       %low
       idxPair = union(find(Aloc==pairs(i,1) & Bloc==pairs(i,2)),...
                       find(Aloc==pairs(i,2) & Bloc==pairs(i,1)));
       
       DecsLow(n,i) = sum(location(intersect(idxPair,idxLow)) == pairs(i,1));
       DecsHigh(n,i) = sum(location(intersect(idxPair,idxHigh)) == pairs(i,1));
       
    end
    
    locsLow(n,:) = [sum(location(idxLow)==1),sum(location(idxLow)==2),...
                    sum(location(idxLow)==3),sum(location(idxLow)==4)];    
    locsHigh(n,:) = [sum(location(idxHigh)==1),sum(location(idxHigh)==2),...
                    sum(location(idxHigh)==3),sum(location(idxHigh)==4)];

    
%% Subjective weight & Model performance




end
end
%%
figure(1); clf;
bar(locsprop);
ylabel('Proportion'); xlabel('Day');
title('Proportion of times selected stimulus locations (Arya)');
l1 = line([0,4.7],[0.58,0.58]);
set(l1,'linewidth',3,'color','r');
l2 = line([18,24],[0.45,0.45]);
set(l2,'linewidth',3,'color','r');
% stimulus positions 1-4; 1 = leftUp, 2 = leftDown, 3 = rightUp, 4 = rightDown
legend('Left Up','Left Down','Right Up','Right Down','Time Pressure (375 ms)');


%%
TPloc = sum(locs([1,2,3,4,19,20,21,22,23],:));
TPstd = std(locs([1,2,3,4,19,20,21,22,23],:));
NPloc = sum(locs(5:18,:));
NPstd = std(locs(5:18,:));
figure(2); clf; bar([NPloc;TPloc]);
hold on;
errorbar([0.73,0.915,1.09,1.28;1.73,1.915,2.09,2.28],[NPloc;TPloc],[NPstd;TPstd],'k.');
set(gca,'XtickLabel',{'No Pressure','Time Pressure'});
legend('Left Up','Left Down','Right Up','Right Down','Location','North');
ylabel('Reponses');
title('Number of times selected stimulus locations (Arya)');
%%
TPlocLow = sum(DecsLow([1,2,3,4,19,20,21,22,23],:))/sum(sum(DecsLow([1,2,3,4,19,20,21,22,23],:)));
TPlocHigh = sum(DecsHigh([1,2,3,4,19,20,21,22,23],:))/sum(sum(DecsHigh([1,2,3,4,19,20,21,22,23],:)));

NPlocLow = sum(DecsLow(5:18,:))/sum(sum(DecsLow(5:18,:)));
NPlocHigh = sum(DecsHigh(5:18,:))/sum(sum(DecsHigh(5:18,:)));
figure(3); clf; bar([NPlocLow;NPlocHigh;TPlocLow;TPlocHigh]);
hold on;
set(gca,'XtickLabel',{'No Pressure','Time Pressure'});
legend('Left Up','Left Down','Right Up','Right Down','Location','North');
ylabel('Reponses');
title('Number of times selected stimulus locations (Arya)');


%%
% TPlocHigh = sum(locsHigh([1,2,3,4,19,20,21,22,23],:))/sum(sum(locsHigh([1,2,3,4,19,20,21,22,23],:)));
% TPstdHigh = std(locsHigh([1,2,3,4,19,20,21,22,23],:));
% TPlocLow = sum(locsLow([1,2,3,4,19,20,21,22,23],:))/sum(sum(locsLow([1,2,3,4,19,20,21,22,23],:)));
% TPstdLow = std(locsLow([1,2,3,4,19,20,21,22,23],:));
% NPlocHigh = sum(locsLow(5:18,:))/sum(sum(locsLow(5:18,:)));
% NPstdHigh = std(locsLow(5:18,:));
% NPlocLow = sum(locsLow(5:18,:))/sum(sum(locsLow(5:18,:)));
% NPstdLow = std(locsLow(5:18,:));
% figure(3); clf; bar([NPlocLow;NPlocHigh;TPlocLow;TPlocHigh]);
% hold on;
% %errorbar([0.73,0.915,1.09,1.28;1.73,1.915,2.09,2.28],[NPloc;TPloc],[NPstd;TPstd],'k.');
% set(gca,'XtickLabel',{'NP Low \Delta','NP High \Delta','TP Low \Delta','TP High \Delta'});
% legend('Left Up','Left Down','Right Up','Right Down','Location','North');
% ylabel('Reponses');
% title('Number of times selected stimulus locations (Arya)');
