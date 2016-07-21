clear all; close all;

%% Edit the following flags to toggle which analysis to perform/figures to plot
plotAll = 1; 
plotIndiv = 0;
CumComp = 0;
choice_flag = 0;
logReg_flag = 0;
subWeight_flag = 0;
model_flag = 0;
responseTime_flag = 0; %No
correctFlag = 1; %Not really
modelCorrFlag = 0;



%% Data files to analyze
 % specifiy dates as 'Mmm dd yyyy'
startDate = 'Jun 11 2015';
stopDate = 'Apr 25 2016';
%dataFolder = 'Arya';
dataFolder = 'Lana';

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
TP=0;
P=0;
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
    responseTimeAll(1:length(responseTime),n) = responseTime;
    trialTimeAll(n) = trialTime(1);
    
    allCorrect(n) = sum(correct)/nValid;
    allHit(n) = (nValid/nResponse);
    allindValid(n)=length(indValid);

    if(trialTime(1)<1)
        TP=TP+1;
    CumCorrTP=ones(size(correct,2));
    CumCueATP=ones(size(correct,2));
    CumCueBTP=ones(size(correct,2));
    CumrespStimTP=ones(size(correct,2));
   % CumcueProbInfo=ones(size(correct,2));
    for i=1:TP
    CumCorrTP = vertcat(CumCorrTP,correct);
    CumCueATP = vertcat(CumCueATP,cueA);
    CumCueBTP = vertcat(CumCueBTP,cueB);
    CumrespStimTP = vertcat(CumrespStimTP,respStimulus);
   % CumcueProbInfo = vertcat(CumcueProbInfo,
    end
    CumCorrTP=CumCorrTP(2:end,:);
    CumCueATP = CumCueATP(2:end,:);
    CumCueBTP = CumCueBTP(2:end,:);
    CumrespStimTP = CumrespStimTP(2:end,:);
    else
    P=P+1;
    CumCorr=ones(size(correct,2));
    CumCueA=ones(size(correct,2));
    CumCueB=ones(size(correct,2));
    CumrespStim=ones(size(correct,2));
   % CumcueProbInfo=ones(size(correct,2));
    for i=1:P
    CumCorr = vertcat(CumCorr,correct);
    CumCueA = vertcat(CumCueA,cueA);
    CumCueB = vertcat(CumCueB,cueB);
    CumrespStim = vertcat(CumrespStim,respStimulus);
   % CumcueProbInfo = vertcat(CumcueProbInfo,
    end
    CumCorr=CumCorr(2:end,:);
    CumCueA = CumCueA(2:end,:);
    CumCueB = CumCueB(2:end,:);
    CumrespStim = CumrespStim(2:end,:);
    end
    
    
    
    
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
%    [sumW,zSumW,threshold,bLog10,BayesFactor] = logRegPerformance([CumCueATP,CumCueBTP],...
%             CumrespStimTP,CumCorrTP,cueProbInfo);        
        
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
            set(gca,'fontsize',35);
            %title(barTitles{i},'fontsize',45);
            text(2,0.95,sprintf('n = %d',nChoice(i)));
            lineH = h.expCond.probDist(i);
            % Place line up to heighest weight for each dimension
            lineH = (lineH<0.5)*(1-lineH) + (lineH>=0.5)*lineH;
            plot([0.5 2.5], [lineH lineH], 'k', 'LineWidth', 3)
        end
        suptitle(titleStr);
        %ylabel('Proportion of Times Picked','fontsize',40)
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
        plot(sumW(:,3),sumW(:,1)./sumW(:,2),'r.','MarkerSize',25); hold on;
        plot(sumW(:,3),zSumW);
        set(gca,'fontsize',15,'linewidth',3);
        xlabel('evidence for cue A','fontsize',25); ylabel('ratio of cue A choices','fontsize',25);
        title(titleStr,'fontsize',30);
        
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
        b = glmfit(orderedWeights, bLog10([5 4 3 2]));
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
   % bf(n,:) = BayesFactor;
    
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

if(CumComp)

if(logReg_flag || subWeight_flag || model_flag)
   [sumW,zSumW,threshold,bLog10,BayesFactor] = logRegPerformance([CumCueA,CumCueB],...
            CumrespStim,CumCorr,cueProbInfo);
   if(strcmp(h.monkey,'Arya'))
   [sumWTP,zSumWTP,thresholdTP,bLog10TP,BayesFactorTP] = logRegPerformance([CumCueATP,CumCueBTP],...
            CumrespStimTP,CumCorrTP,cueProbInfo); 
   end
   
    if(strcmp(h.monkey,'Lana'))
        % Change model order for Lana. Map Lana model order -> Arya model
        % order
        BayesFactor = BayesFactor([4;3;2;1;10;9;7;8;6;5;14;13;12;11;15]);
    end
   % allThresholds(n) = threshold;
    subweight = bLog10([5 4 3 2]);

end

   if(logReg_flag)
        figure(); hold on;
        plot(sumW(:,3),sumW(:,1)./sumW(:,2),'k.','MarkerSize',20); 
        h1=plot(sumW(:,3),zSumW,'k-');
        
        plot(sumWTP(:,3),sumWTP(:,1)./sumWTP(:,2),'rx','MarkerSize',10); 
        h2=plot(sumWTP(:,3),zSumWTP,'r--');
        
        legend([h1 h2],{'no time pressure','time pressure'})
        xlabel('evidence for cue A'); ylabel('ratio of cue A choices');
        T=sprintf('Individual Days of Logistic Regression Model Cumulated from %s to %s',startDate,stopDate);
        title(T); 
        set(gca,'fontsize',15,'linewidth',2);
   end
   
   if(subWeight_flag)
        figure(); hold on;
        orderedWeights = abs(1-2*sort(h.expCond.probDist)); % puts weights in appropraite order
        %(e.g. 0.2 -> 0.8 for Arya, 0.8 -> 0.2 for Lana)
        sde = std(subweightAll)./sqrt(size(subweightAll,2));
        errorbar(orderedWeights,bLog10([5 4 3 2]),sde,'k.','MarkerSize',30);
        if(strcmp(h.monkey,'Arya'))
        plot(orderedWeights,bLog10TP([5 4 3 2]),'r.','MarkerSize',30);
        end
        xlim([0 1]); set(gca,'XTick',0.2:0.2:0.8);
        xlabel('assigned weights'); ylabel('subjective weights');
        b = glmfit(0.2:0.2:0.8, bLog10([5 4 3 2]));
        if(strcmp(h.monkey,'Arya'))
        bTP=glmfit(0.2:0.2:0.8, bLog10TP([5 4 3 2]));
        h2=plot([0,1], bTP(1) + bTP(2).*[0,1],'k:','linewidth',3);
        end
        h1=plot([0,1], b(1) + b(2).*[1,0],'k-','linewidth',3);
        h3=plot([0,1],[0,1],'r--','linewidth',3);
        axis([0 1 -inf inf])
        if(strcmp(h.monkey,'Arya'))
        legend([h1, h2, h3],{'no TP','TP','expected'});
        elseif(strcmp(h.monkey,'Lana'))
        legend([h1 h3],{'monkey performance','ideal'}); 
        else
           warning('Missing: no monkey?');
        end
        %set(gca,'ylim',[0 5],'YTick',0:0.2:1.6);
        set(gca,'fontsize',25,'linewidth',3);
        xlabel('assigned weights','fontsize',35); ylabel('subjective weights','fontsize',35);
        T=sprintf('Subjective Weight from %s to %s',startDate,stopDate);
        title(T,'fontsize',45); 
        hold off
        
        
        %% testing
        
        figure();
        orderedWeights = abs(1-2*sort(h.expCond.probDist)); % puts weights in appropraite order
        %(e.g. 0.2 -> 0.8 for Arya, 0.8 -> 0.2 for Lana)
        sde = std(subweightAll)./sqrt(size(subweightAll,2));
        errorbar(orderedWeights,mean(subweightAll),sde,'k.','MarkerSize',25);
        xlim([0 1]); set(gca,'XTick',0.2:0.2:0.8);
        b = glmfit(orderedWeights, mean(subweightAll));
        hold on;
        h1 = plot([0,1], b(1) + b(2).*[0,1],'k-','linewidth',2);
        h2=plot([0,1],[0,1],'r--','linewidth',2);
        axis([0 1 -inf inf])
        legend([h1 h2],{'monkey performance','ideal'}); 
        set(gca,'fontsize',25,'linewidth',2);
        xlabel('assigned weights','fontsize',35); ylabel('subjective weights','fontsize',35);
        T=sprintf('Individual Days of Subjective Weight Cumulated from %s to %s',startDate,stopDate);
        title(T,'fontsize',45); 
        
   end
    
     if(model_flag)
        figure(); 
        line([0 16.5],[3 3]); hold on;
        h1=bar([BayesFactor' BayesFactorTP']);
        %h2 = bar(BayesFactorTP,'r');
        legend(h1,{'no TP','TP'})
        set(gca,'XTick',1:16,'ylim',[0 10]);
        T=sprintf('Individual Days of Model Comparison Cumulated from %s to %s',startDate,stopDate);
        title(T); 
        set(gca,'fontsize',15,'linewidth',2);
        xlabel('Model number'); ylabel('Bayes factor');
        xlim([0,16.5]);
     end
    
     if(modelCorrFlag)
         
        [rhoOT,~] = getModelCorr2(h.expCond.probDist,CumCueA,CumCueB,CumrespStim);
        if(strcmp(h.monkey,'Lana'))
            rhoOT = rhoOT([4;3;2;1;10;9;7;8;6;5;14;13;12;11;15]);
        end
        
        if(strcmp(h.monkey,'Arya'))
        [rhoTP,~] = getModelCorr2(h.expCond.probDist,CumCueATP,CumCueBTP,CumrespStimTP);
        end
%         if(strcmp(h.monkey,'Lana'))
%             rhoTP = rhoTP([4;3;2;1;10;9;7;8;6;5;14;13;12;11;15]);
%         end 
         
        figure(); 
        line([0 16.5],[3 3],'linewidth',2); hold on;
        h1 = bar([rhoOT' rhoTP']);
        legend(h1,{'no TP','TP'})
        set(gca,'XTick',1:16,'ylim',[0 1]);
        T=sprintf('Individual Days of Model Correlation Cumulated from %s to %s',startDate,stopDate);
        title(T);
        xlabel('Model number'); ylabel('correlation');
        xlim([0,16.5]);
        set(gca,'fontsize',15,'linewidth',2);

    end
end

% if(modelCorrFlag)
%         [rho,p] = getModelCorr2(h.expCond.probDist,cueA,cueB,respStimulus);
%         if(strcmp(h.monkey,'Lana'))
%         % Change model order for Lana. Map Lana model order -> Arya model
%         % order
%             rho = rho([4;3;2;1;10;9;7;8;6;5;14;13;12;11;15]);
%         end
%         rhoAll(n,:) = rho;
% 
% %             [maxRho(n),maxRhoI(n)] = max(rho);
% end


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
            %text(i,1.05*allCorrect(i),sprintf('%i',1000*trialTimeAll(i)));
        end
        
        h1 = plot(1:n,allCorrect,'k.','markersize',10);
        medianCorr = mean(allCorrect);
        h2 = plot(1:0.1:n,medianCorr*ones(length(1:0.1:n)),'r--','linewidth',2);
        set(gca,'xlim',[0,n+1], ...
               'ylim',[0,1],'yTick',[0:0.1:1]);
            %'Xtick',1:length(allDates),...

        set(gca,'fontsize',15,'linewidth',2);
        ylabel('proportion correct','fontsize',25);
        legend(h2,'mean proportion correct');
        xlabel('days ran in order','fontsize',25);
        T=sprintf('Proportion "Correct" from %s to %s',startDate,stopDate);
        title(T,'fontsize',30);         
    end

    if(choice_flag) 
        f2 = figure('Position', [200, 500, 1300, 500]);
        hold on;
        for i = 1:2:length(Matrix)-1
            area([Matrix(i) Matrix(i+1)], [1 1], 'FaceColor', [.9 .9 .9]);           
        end
        for i = 1:n
           % text(i,0.95,sprintf('%i',1000*trialTimeAll(i)));
        end
       
        h1=plot(1:n,allChoiceProb(:,1),'r-o');
        h2=plot(1:n,allChoiceProb(:,2),'b-o');
        h3=plot(1:n,allChoiceProb(:,3),'g-o');
        h4=plot(1:n,allChoiceProb(:,4),'m-o');
        
        plot([0,n],[0.5,0.5],'k');
        set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates), ...
            'ylim',[0.4,1])
        if(strcmp(h.monkey,'Arya'))
        legend([h1 h2 h3 h4],{'Color','Shape','Border','Direction'});
        elseif(strcmp(h.monkey,'Lana'))
        legend([h4 h3 h2 h1],{'Direction','Border','Shape','Color'});
        end
        set(gca,'fontsize',15,'linewidth',2);
        xlabel('days ran in order','fontsize',25);
        ylabel('probability','fontsize',25)
        T=sprintf('Cue Choice Preference from %s to %s',startDate,stopDate);
        title(T,'fontsize',30);  

    end

    if(logReg_flag)
        f3 = figure('Position', [200, 500, 1300, 500]);
        hold on;
        for i = 1:2:length(Matrix)-1
            area([Matrix(i) Matrix(i+1)], [1.1*max(allThresholds) 1.1*max(allThresholds)], 'FaceColor', [.9 .9 .9]);           
        end
        for i = 1:n
            %text(i,1.05*allThresholds(i),sprintf('%i',1000*trialTimeAll(i)));
        end
        Threshx=1:n;
        plot(Threshx,allThresholds,'o');
%         set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),...
%             'ylim',[0,1.1*max(allThresholds)]);

         set(gca,'xlim',[0,n+1],...
             'ylim',[0,1.1*max(allThresholds)]);

        set(gca,'fontsize',25,'linewidth',2);
        xlabel('days ran in order','fontsize',35);
        ylabel('threshold value','fontsize',35);
        T=sprintf('Thresholds from %s to %s',startDate,stopDate);
        title(T,'fontsize',35);        
        hold off
        
        m=1;
        p=1;
        for i=1:n
            if (Mat(i) == 1)
                %text(m,1.05*allThresholds(i),sprintf('%m',1000*trialTimeAll(m)));
                TPthreshold(m) = allThresholds(i);
                TPindValid(m) = allindValid(i);
                m=m+1;
            else
                NPthreshold(p) = allThresholds(i);
                NPindValid(p) = allindValid(i);
                p=p+1;
            end
        end
        TPindValid = TPindValid./max(TPindValid);
        NPindValid = NPindValid./max(NPindValid);
        
        figure();
        hold on
        h1 = plot(1:m-1,TPthreshold,'k--o');
        %h2 = plot(1:m-1,TPindValid,'k-x');
        h3 = plot(1:p-1,NPthreshold,'r-o');
        %h4 = plot(1:p-1,NPindValid,'r-x');
        set(gca,'xlim',[0,m],'Xtick',1:length(m-1),...
            'ylim',[0,1.1*max(TPthreshold)]);
        legend([h1,h3], 'time pressure (TP)','no TP')
        set(gca,'fontsize',15,'linewidth',2);
        xlabel('days ran in order','fontsize',25);
        ylabel('threshold value','fontsize',25)
        T=sprintf('Thresholds from %s to %s',startDate,stopDate);
        title(T,'fontsize',30);        
        hold off
        
        %% testing:
        
        if(test)
        figure();
        TPtrialtime = trialTimeAll(trialTimeAll<2);
        TPthreshold = allThresholds(trialTimeAll<2);
        plot(TPtrialtime,TPthreshold,'k.','markersize',20)
        end
        
        
    end

    if(subWeight_flag)
        f4 = figure('Position', [200, 500, 1300, 500]);
        hold on;
        for i = 1:2:length(Matrix)-1
            area([Matrix(i) Matrix(i+1)], [1.1*max(subweightAll(:)) 1.1*max(subweightAll(:))], 'FaceColor', [.9 .9 .9]);           
        end
     
        h1 = plot(1:n,subweightAll(:,4),'r-o');
        h2=plot(1:n,subweightAll(:,3),'b-o');
        h3=plot(1:n,subweightAll(:,2),'g-o');
        h4=plot(1:n,subweightAll(:,1),'m-o');
        
        %plot([0,n],[0.8,0.8],'b:',[0,n],[0.6,0.6],'r:',[0,n],[0.4,0.4],'g:',...
        %   [0,n],[0.2,0.2],'m:','LineWidth',2);
        set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),...%num2cell(1000*trialTimeAll),...
            'ylim',[-0.3,1.1*max(subweightAll(:))]);
        if(strcmp(h.monkey,'Arya'))
        legend([h1 h2 h3 h4],{'Color','Shape','Border','Direction'});
        elseif(strcmp(h.monkey,'Lana'))
        legend([h4 h3 h2 h1],{'Direction','Border','Shape','Color'});
        end
        set(gca,'fontsize',15,'linewidth',2);
        xlabel('days ran in order','fontsize',25);
        ylabel('subjective weights','fontsize',25);
        T=sprintf('Subjective Cue Weight from %s to %s',startDate,stopDate);
        title(T,'fontsize',30);

        
        m=1;
        p=1;
        for i=1:n
            if (Mat(i) == 1)
                %text(m,1.05*allThresholds(i),sprintf('%m',1000*trialTimeAll(m)));
                TPsubweight(m,1) = subweightAll(i,1);
                TPsubweight(m,2) = subweightAll(i,2);
                TPsubweight(m,3) = subweightAll(i,3);
                TPsubweight(m,4) = subweightAll(i,4);
                m=m+1;
            else
                NPsubweight(p,1) = subweightAll(i,1);
                NPsubweight(p,2) = subweightAll(i,2);
                NPsubweight(p,3) = subweightAll(i,3);
                NPsubweight(p,4) = subweightAll(i,4);
                p=p+1;
            end
        end
        
        figure();
        subplot(2,1,1);
        hold on
        h1 = plot(1:m-1,TPsubweight(:,4),'r-o');
        h2 = plot(1:m-1,TPsubweight(:,3),'b-o');
        h3 = plot(1:m-1,TPsubweight(:,2),'g-o');
        h4 = plot(1:m-1,TPsubweight(:,1),'m-o');
        set(gca,'xlim',[0,m],'Xtick',1:length(m-1),...%num2cell(1000*trialTimeAll),...
            'ylim',[-0.3,1.1*max(subweightAll(:))]);
        if(strcmp(h.monkey,'Arya'))
        legend([h1 h2 h3 h4],{'Color','Shape','Border','Direction'});
       % legend([h5 h6 h7 h8],{'Color','Shape','Border','Direction'});
        elseif(strcmp(h.monkey,'Lana'))
        legend([h4 h3 h2 h1],{'Direction','Border','Shape','Color'});
        %legend([h8 h7 h6 h5],{'Direction','Border','Shape','Color'});
        end
        set(gca,'fontsize',15,'linewidth',2);
        xlabel('days ran in order','fontsize',25);
        ylabel('subjective weights','fontsize',25);
        T=sprintf('Time Pressured Subjective Cue Weight from %s to %s',startDate,stopDate);
        title(T,'fontsize',30);
        hold off
        subplot(2,1,2);
        hold on
        h5 = plot(1:p-1,NPsubweight(:,4),'r--x');
        h6 = plot(1:p-1,NPsubweight(:,3),'b--x');
        h7 = plot(1:p-1,NPsubweight(:,2),'g--x');
        h8 = plot(1:p-1,NPsubweight(:,1),'m--x');
        set(gca,'xlim',[0,p],'Xtick',1:length(p-1),...%num2cell(1000*trialTimeAll),...
            'ylim',[-0.3,1.1*max(subweightAll(:))]);
        if(strcmp(h.monkey,'Arya'))
        %legend([h1 h2 h3 h4],{'Color','Shape','Border','Direction'});
        legend([h5 h6 h7 h8],{'Color','Shape','Border','Direction'});
        elseif(strcmp(h.monkey,'Lana'))
        %legend([h4 h3 h2 h1],{'Direction','Border','Shape','Color'});
        legend([h8 h7 h6 h5],{'Direction','Border','Shape','Color'});
        end
        set(gca,'fontsize',15,'linewidth',2);
        xlabel('days ran in order','fontsize',25);
        ylabel('subjective weights','fontsize',25);
        T=sprintf('No Time Pressured Subjective Cue Weight from %s to %s',startDate,stopDate);
        title(T,'fontsize',30);

        hold off
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
           % text(i,0.95,sprintf('%i',1000*trialTimeAll(i)));
            [~,maxrI(i)] = max(rhoAll(i,:));
        end
        
        h1=plot(1:n,rhoAll(:,5),'k-*',1:n,rhoAll(:,11),'k-o',1:n,rhoAll(:,15),'r-o');
        
         set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates), ...%num2cell(1000*trialTimeAll),...
        'ylim',[0,1]);
        legend(h1,{'Model 5','Model 11','Model 15'});
        set(gca,'fontsize',15,'linewidth',2);
        xlabel('days ran in order','fontsize',25);
        ylabel('Correlation coefficient','fontsize',25);
        T=sprintf('Model Correlation from %s to %s',startDate,stopDate);
        title(T,'fontsize',30);
        ylim([0,1]);
        subplot(2,1,2);
        hist(maxrI,1:15);
        set(gca,'xlim',[0,15],'Xtick',1:15);
        T=sprintf('Max correlation from %s to %s',startDate,stopDate);
        set(gca,'fontsize',15,'linewidth',2);
        xlabel('Model number','fontsize',25)
        ylabel('count','fontsize',25)
        title(T,'fontsize',30);
    end
    
    if(responseTime_flag)
        responseTimeAll(responseTimeAll==0)=NaN;

        f2 = figure('Position', [200, 500, 1800, 500]);
        clf, hold on
        boxplot(responseTimeAll);
        h=findobj(gca,'tag','Outliers');
        delete(h)
        plot(1:n,trialTimeAll,'k--')
        legend('Trial Time');
        ylabel('ms')
        set(gca,'Ytick',0:0.25:1)
        ylim([0.25,2.25])
    end

end



