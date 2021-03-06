function run_analyses(filePath,files,flags)
% Run all analyses
% Designed to run with SatisficingGui, but can be run independently
% filePath - directory path from Analysis folder to desired data folder.
%          e.g. for Arya: '../Data/Arya'
%               for Lana: '../Data/Lana'
% files - cell arya of specific .mat files 
%          e.g. {'Aryasummary_01_Apr_2016fourcueDim_1.mat',...'}
% flags - Specifies which analyses to run and settings. See below

if(nargin == 3)
%% Get flags from GUI
% IMPORTANT: "flags" must match up with appropriate flags in 'satisficing
% GUI.m'
    plotIndiv = flags(1); % plot individual days
    plotAll = flags(2); % plot cumulative data
    choiceFlag = flags(3); % choice probabilities
    logRegFlag = flags(4); % logisitc regression curve
    subWeightFlag = flags(5); % subjective weights
    modelFlag = flags(6); % Bayes model comparison
    responseTimeFlag = flags(7); % response times histogram
    correctFlag = flags(8); % plot proportion correct
    modelCorrFlag = flags(9); % model correlations
    idealFlag = flags(10); % ideal observer
    xdateFlag = flags(11); % label x axis with date vs. number
    dispTTflag = flags(12); % dispplay trial times 
elseif(nargin == 2)
% Else, use defaults
    plotIndiv = 1; % plot individual days
    plotAll = 1; % plot cumulative data
    choiceFlag = 1; % choice probabilities
    logRegFlag = 1; % logisitc regression curve
    subWeightFlag = 1; % subjective weights
    modelFlag = 1; % Bayes model comparison
    responseTimeFlag = 1; % response times histogram
    correctFlag = 1; % plot proportion correct
    modelCorrFlag = 1; % model correlations
    idealFlag = 1; % ideal observer on plots
    xdateFlag = 1; % label x axis with date vs. number
    dispTTflag = 1; % dispplay trial times 
end

%% Initialize main loop
N = length(files);
n=0;
clc;
fprintf('Running analysis on following days:\n');
ignoreList = {};
allRT_G = []; 
allSummaries = [];
for f =1:N
%% Read each file and index according to data
    currentFile = load(fullfile(filePath,files{f}));  
    h = currentFile.header;
    s = currentFile.summary;
    
    %%% -1 => no reponse
    %%% -2 => touched screen, missed stimulus
    %%% -3 => objective weighting of both stimuli were the same
    indResponse = find(s.respLocation ~= -1 &...
                       s.respLocation ~= 0)'; % indeces of all trials where touchscreen was touched
    indMiss = find(s.respLocation == -2)'; % indeces of trials where touchscreen was hit, stimuli were missed
    indValid= find(s.respLocation ~= -2 &...
                   s.respLocation ~= -1 &...
                   s.respLocation ~= 0)'; % indeces of valid trials, where either stimulus was hit
    indCorrect = find(s.correct == 1);
    
    indIgnore = find(s.correctLoc == -3);
    indValid(ismember(indValid,indIgnore)) = []; % exclude trials where both stimuli had equal weighting
    
    indNoGrace = find(s.responseTime<=h.trialTime);
    indWithGrace = find(s.responseTime<=(h.trialTime + h.gracePeriod));
    indValidWithGrace = indValid;
    
    indValid(~ismember(indValid,indNoGrace)) = []; % exclude trials where response came after trial time
    indValidWithGrace(~ismember(indValid,indWithGrace)) = []; % include trials within grace period
    
    nValid = length(indValid);

%% Ignore days with not enough valid responses
if(nValid<30)
    % 30 is a pretty arbitrary amount, however, it appears to be good at
    % preventing warnings from occuring in certain analyses
    ignoreList{end+1} = files{f}; 
else
%% Get date from header and create figure titles
    n=n+1;
    datesOut{n} = h.date;
    % Reformat date for figure titles
    dateI = regexprep(h.date,'_',' ');
    allDates{n} = dateI;
    disp(dateI);
    
    % kesh alert 7/18/2016 - Somehow, experimental code changed from
    % h.subject to h.monkey. Revert back to h.subject to be consistent with
    % the past year of data
    
    if isfield(h,'subject')
        % do nothning
    elseif isfield(h, 'monkey')
        h.subject = h.monkey;
    else
        disp('ruh roh')
    end
    
    titleStr = [h.subject '  ' dateI, sprintf('  (T = %.3f s)',h.trialTime)];
%% Extract data from each specified trial
    % Note, cue probability infor is different for Arya and Lana. Make sure
    % proper file is loaded
    if(strcmp(h.subject,'Arya'))
        cueProbInfo = load('cueProbInfo_Arya.mat');
    elseif(strcmp(h.subject,'Lana'))
        cueProbInfo = load('cueProbInfo_Lana.mat');
    else
        warning('Missing: "cueProbInfo.mat" for specified subject name');
    end
    
    ind = indValid; % Only look at valid trials
    cueA = s.cueA(ind); cueB = s.cueB(ind);
    cueALoc = s.cueALoc(ind); cueBLoc = s.cueBLoc(ind);
    respStimulus = s.respStimulus(ind); 
    correct = s.correct(ind);
    reward = s.reward(ind);
    responseTime = s.responseTime(ind);
    responseTimeWithGrace = s.responseTime(indValidWithGrace);
    trialTime = h.trialTime;
    
    allSummaries(end+1:end+length(ind),:) = s(ind,:);
%% Reponse Time and correctness
    trialTimeAll(n) = trialTime;
    allRT_G(n) = median(responseTimeWithGrace);
    allCorrect(n) = sum(correct)/nValid;

%% Ideal observer responses
respIdeal = zeros(size(respStimulus));
for ri = 1:nValid
    if(reward(ri) == 1)
        respIdeal(ri) = respStimulus(ri);
    else
        temp = [cueA(ri),cueB(ri)];
        respIdeal(ri) = temp(temp~=respStimulus(ri));
    end
end
    allCorrect_i(n) = sum(correct==reward)/nValid;
    
%% Run various analyses
if(choiceFlag)
    [choiceCount,nChoice] = run_cueChoice(cueA,cueB,respStimulus,cueProbInfo);
    choiceProb = choiceCount./nChoice;
    allChoiceProb(n,:) = choiceProb;
    if(strcmp(h.subject,'Lana'))
        allChoiceProb(n,:) = 1-choiceProb;
    end
    
    [choiceCount,nChoice] = run_cueChoice(cueA,cueB,respIdeal,cueProbInfo);
    CP_i = choiceCount./nChoice;
    allChoiceProb_i(n,:) = choiceProb;
    if(strcmp(h.subject,'Lana'))
        allChoiceProb_i(n,:) = 1-CP_i;
    end
    
    assignin('base','allChoiceProb',allChoiceProb);
    assignin('base','allChoiceProb_i',allChoiceProb_i);
end

matTrial = generate_matTrial(cueA,cueB,respStimulus,correct,cueProbInfo);
matTrial_i = generate_matTrial(cueA,cueB,respIdeal,correct,cueProbInfo);

% Psychometric curve
if(logRegFlag)%  || subWeightFlag || modelFlag)
    % Responses
   [sumW,zSumW,threshold] = run_psychCurve(matTrial);
    allThresholds(n) = threshold;
    % Ideal
   [sumW_i,zSumW_i,threshold_i] = run_psychCurve(matTrial);
    allThresholds_i(n) = threshold_i;
    
    assignin('base','allThresholds',allThresholds);
    assignin('base','allThresholds_i',allThresholds_i);
end

% Subjective weights
if(subWeightFlag)
    % Responses
    locMat = [cueALoc,cueBLoc];
    [bNS,bSpatial] = run_subjectiveWeights(matTrial,locMat);
    bNS = bNS([5,4,3,2,1]);
    allSW(n,:) = bNS;
   
    % Ideal
    [bNS_i,bSpatial_i] = run_subjectiveWeights(matTrial_i,locMat);
    bNS_i = bNS_i([5,4,3,2,1]);
    allSW_i(n,:) = bNS_i;
    
    assignin('base','allSW',allSW);
    assignin('base','allSW_i',allSW_i);
end    

% Bayes Model Comparison
if(modelFlag)
    % Note: Order is dependent on how cue weights are set up.
    BF = run_BayesModel(matTrial);
    BF_i = run_BayesModel(matTrial_i);
    if(strcmp(h.subject,'Lana'))
        % Change model order for Lana. Map Lana model order -> Arya model
        % order
        reorderBF = [4;3;2;1;10;9;7;8;6;5;14;13;12;11;15];
        BF = BF(reorderBF);
        BF_i = BF_i(reorderBF);

    end
    allBF(n,:) = BF;
    allBF_i(n,:) = BF_i;
    assignin('base','allBF',allBF);
    assignin('base','allBF_i',allBF_i);

end

% Model Correlations
if(modelCorrFlag)
        [rho,p] = run_modelCorr(h.expCond.probDist,cueA,cueB,respStimulus);
        if(strcmp(h.subject,'Lana'))
        % Change model order for Lana. Map Lana model order -> Arya model
        % order
            rho = rho([4;3;2;1;10;9;7;8;6;5;14;13;12;11;15]);
        end
        allRho(n,:) = rho;   
        assignin('base','allRho',allRho);
end


%% Individual Plots    
if(plotIndiv)
    % Plot cue choice proportions
    if(choiceFlag) 
        figure();
        barLabels = {'Blue','Red';'Circle','Square';'White','Black';'Horizontal','Vertical'};
        barTitles = {'Color','Shape','Border','Orienation'};
        for ci=1:4 % index through cues
            subplot(2,2,ci), hold on
            bar([choiceProb(ci);1-choiceProb(ci)])
            ylim([0 1]);
            set(gca, 'XTick', [1 2], 'XTickLabel', barLabels(ci,:))
            title(barTitles{ci});
            if(ci == 1)
            ylabel('Proportion of Times Picked')
            end
            text(2,0.95,sprintf('n = %d',nChoice(ci)),'fontsize',14);
            if(idealFlag)
              lineH = CP_i(ci);
              lineH = (lineH<0.5).*(1-lineH) + (lineH>=0.5).*lineH;
              plot([0.5 2.5], [lineH lineH], 'k', 'LineWidth', 3)
              if(ci==1), text(1.6,0.8,'-- Ideal Observer'),end;
            else
              lineH = cueProbInfo.expCond.probDist(ci);
              lineH = (lineH<0.5).*(1-lineH) + (lineH>=0.5).*lineH;
              plot([0.5 2.5], [lineH lineH], 'k', 'LineWidth', 3)
              set(gca,'fontsize',13,'linewidth',2);
              if(ci==1), text(1.6,0.8,'-- Assigned Weight'),end;
            end
        end
        st = suptitle(titleStr);
        set(st,'fontsize',15,'fontweight','bold');
    end

   % Log regression 
   if(logRegFlag)
        figure();
        plot(sumW(:,3),sumW(:,1)./sumW(:,2),'b.','MarkerSize',25); hold on;
        plot(sumW(:,3),zSumW,'k','linewidth',2);
        xlabel('evidence for cue A'); ylabel('ratio of cue A choices');
        title(titleStr);
        set(gca,'fontsize',15,'linewidth',2);
        
        if(idealFlag)
            plot(sumW_i(:,3),zSumW_i,'k:','linewidth',2);
        end
   end

    % Subjective Weight
    if(subWeightFlag)
        figure();
        orderedWeights = abs(1-2*sort(h.expCond.probDist));
        % puts weights in appropraite order
        %(e.g. 0.2 -> 0.8 for Arya, 0.8 -> 0.2 for Lana)
        plot(orderedWeights,bLog10([5 4 3 2]),'b.','markersize',40);
        xlim([0 1]); set(gca,'XTick',0.2:0.2:0.8);
        xlabel('assigned weights'); ylabel('subjective weights');
        b = glmfit(0.2:0.2:0.8, bLog10([5 4 3 2]));
        hold on;
        plot([0,1].*(strcmp(h.subject,'Arya'))+[1,0].*(strcmp(h.subject,'Lana')),...
            b(1) + b(2).*[0,1],'k','LineWidth',2)
        plot(0:1,0:1,'k:','LineWidth',2);
        title(titleStr);
        set(gca,'fontsize',15,'linewidth',2);
       
        if(idealFlag)
            plot(orderedWeights,bLog10_i([5 4 3 2]),'r.','markersize',40);
            bi = glmfit(0.2:0.2:0.8, bLog10_i([5 4 3 2]));
            hold on;
            plot([0,1].*(strcmp(h.subject,'Arya'))+[1,0].*(strcmp(h.subject,'Lana')),...
            bi(1) + bi(2).*[0,1],'r','LineWidth',2) 
        end
        
    end

    % Model performance Bayes
    if(modelFlag)
        figure(); 
        line([0 16.5],[3 3],'linewidth',2); hold on;
        bar(allBF(n,:));
        set(gca,'XTick',1:16,'ylim',[0 10]);
        title(titleStr);
        xlabel('Model number'); ylabel('Bayes factor');
        xlim([0,16.5]);
        set(gca,'fontsize',15,'linewidth',2);
     
    end
    
    % Model correlation
    if(modelCorrFlag)
        figure(); 
        if(idealFlag)
            subplot(1,2,1);  
        end
        bar(rho);
        set(gca,'XTick',1:16,'ylim',[0 1]);
        title(titleStr);
        xlabel('Model number'); ylabel('correlation');
        xlim([0,16.5]);
        set(gca,'fontsize',15,'linewidth',2);
        if(idealFlag)
            %KESH ALERT: This part is not working 7/21/2016
            subplot(1,2,2);  
            bar(rho_i);
            set(gca,'XTick',1:16,'ylim',[0 1]);
            xlabel('Model number'); ylabel('correlation');
            xlim([0,16.5]);
            set(gca,'fontsize',15,'linewidth',2);
        end

        
    end
    
    % Response Time
    if(responseTimeFlag)
        figure(); 
        hist(responseTimeWithGrace*1000,100); title(titleStr);
        xlabel('Response time (ms)');
        ymax = max(get(gca,'ylim'));
        l=line(1000*[trialTime,trialTime],[0,ymax]);
        set(l,'linewidth',2,'color','k');
        ylim([0,ymax]);
        xlim(1000*[0,trialTime + h.gracePeriod]);
        set(gca,'fontsize',15,'linewidth',2);
    end
        
    %%%%
    end

end


end
%% Return variables to base workspace
assignin('base','n',n);
assignin('base','allDates',allDates);
allSummaries = mat2dataset(allSummaries,'VarNames',{'trial','cueA','cueB','respStimulus',...
    'respLocation','correct','reward','cueALoc','cueBLoc','correctLoc',...
    'x','y','responseTime','trialTime'});
assignin('base','allSummaries',allSummaries);


%% Cumulative Plots
if(plotAll)
    if(xdateFlag)
        xStr = allDates;
    else
        xStr = num2cell(1:n);
    end
    TPcutoff = 2;
    NPidx = trialTime>=TPcutoff;
    TPidx = trialTime<TPcutoff;
    
    if(correctFlag)
        f1 = figure('Position', [200, 500, 1300, 500]);
        plot(1:n,allCorrect,'k.','LineWidth',2,'MarkerSize',20); hold on;
        if(dispTTflag)
            for ci = 1:n
                text(ci,1.05*allCorrect(ci),sprintf('%i',1000*trialTimeAll(ci)));
            end
        end
        legend('Correct');
        xl=xlabel('Date'); set(xl,'fontsize',25);
        yl=ylabel('Proportion "Correct"'); set(yl,'fontsize',25);

        if(idealFlag);
           % Ideal observer is correct 79% of the time
           plot(1:n,ones(1,n)*0.79,':');
        end
        
        set(gca,'xlim',[0,n+1],'ylim',[0.5,1],'yTick',0.5:0.1:1);
        set(gca,'FontSize',15);
        if(xdateFlag);
            set(gca,'Xtick',1:n,'XtickLabel',allDates);
        end
    end

    if(choiceFlag) 
        f2 = figure('Position', [200, 500, 1300, 500]);
        plot(1:n,allChoiceProb(:,1),'r-o',1:n,allChoiceProb(:,2),'b-o',...
            1:n,allChoiceProb(:,3),'g-o',1:n,allChoiceProb(:,4),'m-o');
        hold on;
        if(dispTTflag)
            for ci = 1:n
                text(ci,0.95,sprintf('%i',1000*trialTimeAll(ci)));
            end
        end
        plot([0,n],[0.5,0.5],'k');
        legend('Color','Shape','Border','Direction');
        xlabel('Trial Time (ms)');
        title('Cue Choice Preference');
        set(gca,'xlim',[0,n+1],'ylim',[0.5,1],'yTick',0.5:0.1:1);
        set(gca,'FontSize',15);
        if(xdateFlag);
            set(gca,'Xtick',1:n,'XtickLabel',allDates);
        end
        
        f2b = figure('Position', [400 200 800 600]);
        barLabels = {'Blue','Red';'Circle','Square';'White','Black';'Horizontal','Vertical'};
        barTitles = {'Color','Shape','Border','Orienation'};
        for ci=1:4 % index through cues
            subplot(2,2,ci), hold on
            bar([mean(allChoiceProb(:,ci));mean(1-allChoiceProb(:,ci))])
            errorbar([mean(allChoiceProb(:,ci));1-mean(allChoiceProb(:,ci))],...
                [std(allChoiceProb(:,ci));std(allChoiceProb(:,ci))],'.');
            ylim([0 1]);
            set(gca, 'XTick', [1 2], 'XTickLabel', barLabels(ci,:))
            title(barTitles{ci});
            if(ci == 1)
            ylabel('Proportion of Times Picked')
            end
            set(gca,'fontsize',13,'linewidth',2);
        end
        title('Averaged cue choice probabilites');
        text(2,0.95,sprintf('n = %d',size(allChoiceProb,1)),'fontsize',14);
    end

    if(logRegFlag)
        f3 = figure('Position', [200, 500, 1300, 500]);
        plot(1:n,allThresholds,'k.','markersize',20); hold on;
        if(dispTTflag)
            for ci = 1:n
                text(ci,1.05*allThresholds(ci),sprintf('%i',1000*trialTimeAll(ci)));
            end
        end
        xl=xlabel('Day'); set(xl,'fontsize',25);
        yl= ylabel('Thresholds');set(yl,'fontsize',25);
        title('Psychometric thresholds');
        if(idealFlag)
           plot(1:n,ones(1,n)*(mean(allThresholds_i)),'k:','linewidth',2);         
        end
        
        set(gca,'xlim',[0,n+1]);
        set(gca,'FontSize',15);
        if(xdateFlag);
            set(gca,'Xtick',1:n,'XtickLabel',allDates);
        end
    end

    if(subWeightFlag)
        f4 = figure('Position', [200, 500, 1300, 500]);
        plot(1:n,allSW(:,4),'r-o',1:n,allSW(:,3),'b-o',...
            1:n,allSW(:,2),'g-o',1:n,allSW(:,1),'m-o');
        if(dispTTflag)
            for ci = 1:n
                text(ci,1.05*max(allSW(:)),sprintf('%i',1000*trialTimeAll(ci)));
            end
        end

        legend('Color', 'Shape', 'Border','Direction');
        xlabel('Trial Time (ms)');
        title('Subjective Cue Weight');
        
        set(gca,'xlim',[0,n+1],'ylim',[-0.3,1.1*max(allSW(:))]);
        set(gca,'FontSize',15);
        if(xdateFlag);
            set(gca,'Xtick',1:n,'XtickLabel',allDates);
        end
    end

    if(responseTimeFlag)
        f7 = figure('Position', [200, 500, 1800, 500]);
        plot(1:n,allRT_G,'k','LineWidth',2); hold on;
        plot(1:n,trialTimeAll,'k:');
        set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),'XtickLabel',xStr);
        legend('Median Response Time','Trial Time');
        xlabel('Date');
        title('Time (s)');
    end
    
    if(modelFlag)
%         f5 = figure('Position', [200, 500, 1800, 800]);
%         subplot(2,1,1);
%         plot(1:n,allBF(:,5),'k-*',1:n,allBF(:,11),'k-o');
%         if(dispTTflag)
%             for i = 1:n
%                 text(i,3.1,sprintf('%i',1000*trialTimeAll(i)));
%                 [~,maxbI(i)] = max(allBF(i,:));
% 
%             end
%         end
%         line([0 n+1],[3 3],'linewidth',2); hold on;
%          set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),'XtickLabel',xStr,...%num2cell(1000*trialTimeAll),...
%         'ylim',[-0.3,1.1*max(allSW(:))]);
%         legend('Model 5','Model 11');
%         title('Bayes Model Comparison');
%         ylabel('Bayes Factor (normalized to 15)');
%         ylim([0,15]);
%         subplot(2,1,2);
%         hist(maxbI,1:15);
%         set(gca,'xlim',[0,15],'Xtick',1:15);
%         title('Max Bayes Model');

        for bi = 1:n
            
            % kesh alert 7/18/2016 - unknown use of 'bf'. Changed to 'allBF'
            bfcount(bi,:) = allBF(bi,:)>3;
            if(sum(bfcount(bi,:))==0)
                bfcount(bi,15) = 1;
            elseif(sum(bfcount(bi,:))>1)
                
               bfcount(bi,:) = allBF(bi,:)== max(allBF(bi,:));
            end
        end

        figure()
        bar(sum(bfcount));
        set(gca,'XTick',1:16);
        title(titleStr);
        xlabel('Strategy model number'); ylabel('Days most significant');
        xlim([0,16.5]);
        set(gca,'fontsize',15,'linewidth',2);
    end 

    if(modelCorrFlag)
        f6 = figure('Position', [200, 500, 1800, 800]);
        subplot(2,1,1);
        plot(1:n,allRho(:,5),'k-*',1:n,allRho(:,11),'k-o',1:n,allRho(:,15),'r-o');
        if(dispTTflag)
            for ci = 1:n
                text(ci,0.95,sprintf('%i',1000*trialTimeAll(ci)));
                [~,maxrI(ci)] = max(allRho(ci,:));
            end
        end
         set(gca,'xlim',[0,n+1],'XTick',1:length(allDates),'XtickLabel',xStr, 'ylim',[0,1]);
    %num2cell(1000*trialTimeAll),...
        legend('Model 5','Model 11','Model 15');
        title('Model Correlation');
        ylabel('Correlation coefficient');
        ylim([0,1]);
        subplot(2,1,2);
        hist(maxrI,1:15);
        set(gca,'xlim',[0,15],'Xtick',1:15);
        title('Max correlation');
    end


end

%% Display list of files that were skipped
if(length(ignoreList) >0)
    fprintf('\nFinished. Ignored following files (nValid<30):\n');
    for ci = 1:length(ignoreList)
        disp(ignoreList{ci});
    end
end
