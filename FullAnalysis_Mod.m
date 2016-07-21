clear all; close all;

%% Edit the following flags to toggle which analysis to perform/figures to plot
plotAll = 0; 
plotIndiv = 1;
choice_flag = 0;
logReg_flag = 0;
subWeight_flag = 0;
model_flag = 0;
responseTime_flag = 0; %No
correctFlag = 0;
modelCorrFlag = 1;

%% Data files to analyze
 % specifiy dates as Mmm dd yyyy'
startDate = 'Mar 1 2016';
stopDate = 'Apr 4 2016';
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

N = length(files);
n=0;
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
    indTime = find(s.responseTime<=h.trialTime);
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
    trialTime = s.trialTime(ind);
    locs(n,:) = [sum(s.respLocation(ind)==1),sum(s.respLocation(ind)==2),sum(s.respLocation(ind)==3),sum(s.respLocation(ind)==4)]/length(ind);
    if(responseTime_flag)
        figure(); hist(responseTime*1000,100);title(titleStr);
        xlabel('days in order');
        line([trialTime(1),trialTime(1)],[0,10]);
    end
    allMedResponseTimes(n) = median(responseTime);
    trialTimeAll(n) = trialTime(1);
    
    allCorrect(n) = sum(correct)/nValid;
    allHit(n) = (nValid/nResponse);   

    
%% Calculate cue choice proportions


if(choice_flag)
    [choiceCount,nChoice] = cueChoice(stimprop,[cueA,cueB],respStimulus);
    choiceProb = choiceCount./nChoice;
    allChoiceProb(n,:) = choiceProb;
    if(strcmp(h.monkey,'Lana'))
        allChoiceProb(n,:) = 1-choiceProb;
    end
end

%% Logistic regresison, Subjective weight & Model performance
if(logReg_flag || subWeight_flag || model_flag)
   [sumW,zSumW,threshold,bLog10,BayesFactor] = logRegPerformance([cueA,cueB],...
            respStimulus,correct,cueProbInfo);
    allThresholds(n) = threshold;
    subweightAll(n,:) = bLog10([5 4 3 2]);
    if(strcmp(h.monkey,'Lana'))
        % Change model order for Lana. Map Lana model order -> Arya model
        % order
        BayesFactor = BayesFactor([4;3;2;1;10;9;7;8;6;5;14;13;12;11;15]);
    end
    bf(n,:) = BayesFactor;

end

if(modelCorrFlag)
        [rho,p] = getModelCorr2(h.expCond.probDist,cueA,cueB,respStimulus);
        if(strcmp(h.monkey,'Lana'))
        % Change model order for Lana. Map Lana model order -> Arya model
        % order
            rho = rho([4;3;2;1;10;9;7;8;6;5;14;13;12;11;15]);
        end
        rhoAll(n,:) = rho;

%             [maxRho(n),maxRhoI(n)] = max(rho);
end

if(plotIndiv)
    % Plot cue choice proportions
    if(choice_flag) 
        figure();
        barLabels = {'Blue','Red';'Circle','Square';'White','Black';'Horizontal','Vertical'};
        barTitles = {'Color','Shape','Border','Orienation'};
        for i=1:4 % index through cues
            subplot(2,2,i), hold on
            bar([choiceProb(i);1-choiceProb(i)])
            ylim([0 1]);
            set(gca, 'XTick', [1 2], 'XTickLabel', barLabels(i,:))
            title(barTitles{i});
            ylabel('Proportion of Times Picked')
            text(2,0.95,sprintf('n = %d',nChoice(i)));
            lineH = h.expCond.probDist(i);
            % Place line up to heighest weight for each dimension
            lineH = (lineH<0.5)*(1-lineH) + (lineH>=0.5)*lineH;
            plot([0.5 2.5], [lineH lineH], 'k', 'LineWidth', 3)
        end
        suptitle(titleStr);
    end
    
    rho = getModelCorr(h.expCond.probDist,cueA,cueB,respStimulus);
    [maxRho(n),maxRhoI(n)] = max(rho);
    if(n == 2)
       figure; bar(rho); drawnow;
       title(titleStr)
    end
    
%% Logistic regresison, Subjective weight & Model performance
   [sumW,zSumW,threshold,bLog10,BayesFactor] = logRegPerformance([cueA,cueB],...
            respStimulus,correct,cueProbInfo);
        
   % Log regression 
   if(logReg_flag)
        figure();
        plot(sumW(:,3),sumW(:,1)./sumW(:,2),'r.','MarkerSize',20); hold on;
        plot(sumW(:,3),zSumW);
        xlabel('evidence for cue A'); ylabel('ratio of cue A choices');
        title(titleStr);
   end
    
   % Psychometric threshold
    allThresholds(n) = threshold;

    % Subjective Weight
    if(subWeight_flag)
        figure();
        orderedWeights = abs(1-2*sort(h.expCond.probDist)); % puts weights in appropraite order
        %(e.g. 0.2 -> 0.8 for Arya, 0.8 -> 0.2 for Lana)
        plot(orderedWeights,bLog10([5 4 3 2]),'o');
        xlim([0 1]); set(gca,'XTick',0.2:0.2:0.8);
        xlabel('assigned weights'); ylabel('subjective weights');
        b = glmfit(0.2:0.2:0.8, bLog10([5 4 3 2]));
        hold on;
        plot([0,1], b(1) + b(2).*[0,1])

        %set(gca,'ylim',[0 5],'YTick',0:0.2:1.6);
        title(titleStr);
        set(gca,'fontsize',15,'linewidth',2);
    end
    subweightAll(n,:) = bLog10([5 4 3 2]);
    
    % Model performance
    if(model_flag)
        figure(); 
        line([0 16.5],[3 3]); hold on;
        bar(BayesFactor);
        set(gca,'XTick',1:16,'ylim',[0 10]);
        title(titleStr);
        set(gca,'fontsize',15,'linewidth',2);
        xlabel('Model number'); ylabel('Bayes factor');
        xlim([0,16.5]);
    end
    bf(n,:) = BayesFactor;
    
        % Model correlation
    if(modelCorrFlag)
        figure(); 
        line([0 16.5],[3 3],'linewidth',2); hold on;
        bar(rho);
        set(gca,'XTick',1:16,'ylim',[0 1]);
        title(titleStr);
        xlabel('Model number'); ylabel('correlation');
        xlim([0,16.5]);
        set(gca,'fontsize',15,'linewidth',2);

    end
    
end

end

end

%% Cumulative Plots
if(plotAll)
    TPcutoff = 2;
    NPidx = trialTime>=TPcutoff;
    TPidx = trialTime<TPcutoff;
    
        for i = 1:n
            if (trialTimeAll(i)<2)
               Mat(i)=1;
            else
                Mat(i)=0;
            end
        end
        
        p=find(abs(diff(Mat))==1);
        Matrix = p+0.5;

        if(mod(length(Matrix),2))
            Matrix = [0, Matrix];
        else  
             Matrix = [0, Matrix, find(Mat,1,'last')+0.5];
        end    
    
    if(correctFlag)
        f1 = figure('Position', [200, 500, 1300, 500]);
        hold on;
        for i = 1:2:length(Matrix)-1
            area([Matrix(i) Matrix(i+1)], [1 1], 'FaceColor', [.9 .9 .9]);           
        end
        for i = 1:n
            text(i,1.05*allCorrect(i),sprintf('%i',1000*trialTimeAll(i)));
        end
        
        h1 = plot(1:n,allCorrect,'k-o','LineWidth',2);
        set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),...
            'ylim',[0.5,1],'yTick',[0.5:0.1:1]);
        legend(h1,'Correct');
        xlabel('days ran in order');
        T=sprintf('Proportion "Correct" from %s to %s',startDate,stopDate);
        title(T);         
        set(gca,'fontsize',15,'linewidth',2);
    end

    if(choice_flag) 
        f2 = figure('Position', [200, 500, 1300, 500]);
        hold on;
        for i = 1:2:length(Matrix)-1
            area([Matrix(i) Matrix(i+1)], [1 1], 'FaceColor', [.9 .9 .9]);           
        end
        for i = 1:n
            text(i,0.95,sprintf('%i',1000*trialTimeAll(i)));
        end
       
        h1=plot(1:n,allChoiceProb(:,1),'r-o',1:n,allChoiceProb(:,2),'b-o',...
            1:n,allChoiceProb(:,3),'g-o',1:n,allChoiceProb(:,4),'m-o');
        
        plot([0,n],[0.5,0.5],'k');
        set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates), ...
            'ylim',[0.4,1])
        legend(h1,{'Color','Shape','Border','Direction'});
        xlabel('days ran in order');
        ylabel('probability')
        T=sprintf('Cue Choice Preference from %s to %s',startDate,stopDate);
        title(T); 
        set(gca,'fontsize',15,'linewidth',2); 

    end

    if(logReg_flag)
        f3 = figure('Position', [200, 500, 1300, 500]);
        hold on;
        for i = 1:2:length(Matrix)-1
            area([Matrix(i) Matrix(i+1)], [1.1*max(allThresholds) 1.1*max(allThresholds)], 'FaceColor', [.9 .9 .9]);           
        end
        for i = 1:n
            text(i,1.05*allThresholds(i),sprintf('%i',1000*trialTimeAll(i)));
        end
        
        plot(1:n,allThresholds,'k-o');
        set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),...
            'ylim',[0,1.1*max(allThresholds)]);
        xlabel('days ran in order');
        T=sprintf('Thresholds from %s to %s',startDate,stopDate);
        title(T);        
        set(gca,'fontsize',15,'linewidth',2);
    end

    if(subWeight_flag)
        f4 = figure('Position', [200, 500, 1300, 500]);
        hold on;
        for i = 1:2:length(Matrix)-1
            area([Matrix(i) Matrix(i+1)], [1.1*max(subweightAll(:)) 1.1*max(subweightAll(:))], 'FaceColor', [.9 .9 .9]);           
        end
        for i = 1:n
            text(i,1.05*max(subweightAll(:)),sprintf('%i',1000*trialTimeAll(i)));
        end
     
        h1 = plot(1:n,subweightAll(:,4),'r-o',1:n,subweightAll(:,3),'b-o',...
            1:n,subweightAll(:,2),'g-o',1:n,subweightAll(:,1),'m-o');
        
        %plot([0,n],[0.8,0.8],'b:',[0,n],[0.6,0.6],'r:',[0,n],[0.4,0.4],'g:',...
         %   [0,n],[0.2,0.2],'m:','LineWidth',2);
        set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),...%num2cell(1000*trialTimeAll),...
            'ylim',[-0.3,1.1*max(subweightAll(:))]);
        legend(h1,{'Color', 'Shape', 'Border','Direction'});
        xlabel('days ran in order');
        T=sprintf('Subjective Cue Weight from %s to %s',startDate,stopDate);
        title(T);
        set(gca,'fontsize',15,'linewidth',2);
    end

%     if(responseTimeFlag)
%         f5 = figure('Position', [200, 500, 1800, 500]);
%         plot(1:n,allMedResponseTimes,'k','LineWidth',2); hold on;
%         plot(1:n,trialTimeAll,'k:');
%         set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),'XtickLabel',allDates);
%         legend('Median Response Time','Trial Time');
%         xlabel('Date');
%         title('Time (s)');
%     end
    if(model_flag)
        f5 = figure('Position', [200, 500, 1800, 800]);
        %subplot(2,1,1);
        hold on
        for i = 1:2:length(Matrix)-1
            area([Matrix(i) Matrix(i+1)], [15 15], 'FaceColor', [.9 .9 .9]);           
        end
        
        for i = 1:n
            text(i,3.1,sprintf('%i',1000*trialTimeAll(i)));
            [~,maxbI(i)] = max(bf(i,:));
        end
        
        h1 = plot(1:n,bf(:,5),'k-*',1:n,bf(:,11),'k-o');
        
        line([0 n+1],[3 3],'linewidth',2); 
         set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates), ...%num2cell(1000*trialTimeAll),...
        'ylim',[-0.3,1.1*max(subweightAll(:))]);
        legend(h1,{'Model 5','Model 11'});
        T=sprintf('Bayes Model Comparison from %s to %s',startDate,stopDate);
        title(T);
        set(gca,'fontsize',15,'linewidth',2);
        xlabel('days ran in order');
        ylabel('Bayes Factor (normalized to 15)');
        ylim([0,15]);
        %subplot(2,1,2);
        MaxBayesMod=hist(maxbI,1:15);
        %set(gca,'xlim',[0,15],'Xtick',1:15);
        %title('Max Bayes Model');
    end 

    if(modelCorrFlag)
        f6 = figure('Position', [200, 500, 1800, 800]);
        subplot(2,1,1); hold on;
        
        for i = 1:2:length(Matrix)-1
            area([Matrix(i) Matrix(i+1)], [1 1], 'FaceColor', [.9 .9 .9]);           
        end
        
        for i = 1:n
            text(i,0.95,sprintf('%i',1000*trialTimeAll(i)));
            [~,maxrI(i)] = max(rhoAll(i,:));
        end
        
        h1=plot(1:n,rhoAll(:,5),'k-*',1:n,rhoAll(:,11),'k-o',1:n,rhoAll(:,15),'r-o');
        
         set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates), ...%num2cell(1000*trialTimeAll),...
        'ylim',[0,1]);
        legend(h1,{'Model 5','Model 11','Model 15'});
        T=sprintf('Model Correlation from %s to %s',startDate,stopDate);
        title(T);
        set(gca,'fontsize',15,'linewidth',2);
        xlabel('days ran in order');
        ylabel('Correlation coefficient');
        ylim([0,1]);
        subplot(2,1,2);
        hist(maxrI,1:15);
        set(gca,'xlim',[0,15],'Xtick',1:15);
        T=sprintf('Max correlation from %s to %s',startDate,stopDate);
        xlabel('Model number')
        ylabel('count')
        title(T);
        set(gca,'fontsize',15,'linewidth',2);
    end


end



% bar(locs);
% legend



