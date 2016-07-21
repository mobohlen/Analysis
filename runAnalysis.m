function datesOut = runAnalysis(filePath,files,flags)
% Run analysis

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


%%
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

%% Initialize main loop
N = length(files);
n=0;
clc;
fprintf('Running analysis on following days:\n');
ignoreList = {};
allRT_G = []; 
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
    titleStr = [h.monkey '  ' dateI, sprintf('  (T = %.3f s)',h.trialTime)];
%% Extract data from each specified trial
    % Note, cue probability infor is different for Arya and Lana. Make sure
    % proper file is loaded
    if(strcmp(h.monkey,'Arya'))
        cueProbInfo = load('cueProbInfo.mat');
    elseif(strcmp(h.monkey,'Lana'))
        cueProbInfo = load('cueProbInfo_Lana.mat');
    else
        warning('Missing: "cueProbInfo.mat" for specified subject name');
    end
    
    ind = indValid; % Only look at valid trials
    cueA = s.cueA(ind); cueB = s.cueB(ind);
    respStimulus = s.respStimulus(ind); 
    correct = s.correct(ind);
    reward = s.reward(ind);
    responseTime = s.responseTime(ind);
    responseTimeWithGrace = s.responseTime(indValidWithGrace);
    trialTime = h.trialTime;

%% Reponse Time and correctness
    trialTimeAll(n) = trialTime;
    allRT_G(n) = median(responseTimeWithGrace);
    allCorrect(n) = sum(correct)/nValid;

%% Ideal observer
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
%% Calculate cue choice proportions
if(choiceFlag)
    [choiceCount,nChoice] = cueChoice(stimprop,[cueA,cueB],respStimulus);
    choiceProb = choiceCount./nChoice;
    allChoiceProb(n,:) = choiceProb;
    if(strcmp(h.monkey,'Lana'))
        allChoiceProb(n,:) = 1-choiceProb;
    end
    
    if(idealFlag) % Repeat for Ideal observer
        [choiceCount,nChoice] = cueChoice(stimprop,[cueA,cueB],respIdeal);
            CP_i = choiceCount./nChoice;
            allCP_i(n,:) = choiceProb;
        if(strcmp(h.monkey,'Lana'))
            allCP_i(n,:) = 1-CP_i;
        end
    end    
end

%% Logistic regresison, Subjective weight & Model performance
if(logRegFlag || subWeightFlag || modelFlag)
   [sumW,zSumW,threshold,bLog10,BayesFactor,L] = logRegPerformance([cueA,cueB],...
            respStimulus,correct,cueProbInfo);
    allThresholds(n) = threshold;
    subweightAll(n,:) = bLog10([5 4 3 2]);
    if(strcmp(h.monkey,'Lana'))
        % Change model order for Lana. Map Lana model order -> Arya model
        % order
        BayesFactor = BayesFactor([4;3;2;1;10;9;7;8;6;5;14;13;12;11;15]);
    end
    bf(n,:) = BayesFactor;
    Lall(n,:) = L;
    
    if(idealFlag) % Repeat for Ideal observer
           [sumW_i,zSumW_i,threshold_i,bLog10_i,BayesFactor_i] = logRegPerformance([cueA,cueB],...
            respIdeal,correct,cueProbInfo);
        allThresholds_i(n) = threshold_i;
        allSW_i(n,:) = bLog10_i([5 4 3 2]);
        if(strcmp(h.monkey,'Lana'))
            % Change model order for Lana. Map Lana model order -> Arya model
            % order
            BayesFactor_i = BayesFactor_i([4;3;2;1;10;9;7;8;6;5;14;13;12;11;15]);
        end
        bf(n,:) = BayesFactor_i;
    end

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

%% Individual Plots    
if(plotIndiv)
    % Plot cue choice proportions
    if(choiceFlag) 
        figure();
        barLabels = {'Blue','Red';'Circle','Square';'White','Black';'Horizontal','Vertical'};
        barTitles = {'Color','Shape','Border','Orienation'};
        for i=1:4 % index through cues
            subplot(2,2,i), hold on
            bar([choiceProb(i);1-choiceProb(i)])
            ylim([0 1]);
            set(gca, 'XTick', [1 2], 'XTickLabel', barLabels(i,:))
            title(barTitles{i});
            if(i == 1)
            ylabel('Proportion of Times Picked')
            end
            text(2,0.95,sprintf('n = %d',nChoice(i)),'fontsize',14);
%             lineH = h.expCond.probDist(i);
%             % Place line up to heighest weight for each dimension
%             lineH = (lineH<0.5)*(1-lineH) + (lineH>=0.5)*lineH;
            if(idealFlag)
              lineH = CP_i(i);
              lineH = (lineH<0.5)*(1-lineH) + (lineH>=0.5)*lineH;
              plot([0.5 2.5], [lineH lineH], 'k', 'LineWidth', 3)
            end
            set(gca,'fontsize',13,'linewidth',2);

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
        plot([0,1].*(strcmp(h.monkey,'Arya'))+[1,0].*(strcmp(h.monkey,'Lana')),...
            b(1) + b(2).*[0,1],'k','LineWidth',2)
        plot(0:1,0:1,'k:','LineWidth',2);
        %set(gca,'ylim',[0 5],'YTick',0:0.2:1.6);
        title(titleStr);
        set(gca,'fontsize',15,'linewidth',2);
       
        if(idealFlag)
            plot(orderedWeights,bLog10_i([5 4 3 2]),'r.','markersize',40);
            bi = glmfit(0.2:0.2:0.8, bLog10_i([5 4 3 2]));
            hold on;
            plot([0,1].*(strcmp(h.monkey,'Arya'))+[1,0].*(strcmp(h.monkey,'Lana')),...
            bi(1) + bi(2).*[0,1],'r','LineWidth',2) 
        end
        
    end

    % Model performance Bayes
    if(modelFlag)
        figure(); 
        line([0 16.5],[3 3],'linewidth',2); hold on;
        bar(bf(n,:));
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
    %%%
        stims = stimprop(respStimulus,:);
        [r,p] = corr(stims(:,1),stims(:,2));

% KeshAlert - Doesn't look like rNP and rTP are being use.        
%         if(trialTime > 1)
%             rNP(end+1) = r;
%         else
%             rTP(end+1) = r;
%         end
        
        figure(200);
        plot(n,r,'b.','markersize',30); hold on;
        text(n,-0.22,sprintf('%d',1000*trialTime));
        
    %%%%
    end

end


end
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
            for i = 1:n
                text(i,1.05*allCorrect(i),sprintf('%i',1000*trialTimeAll(i)));
            end
        end
        set(gca,'xlim',[0,n+1],'Xtick',0:10:length(allDates),'XtickLabel',0:10:length(allDates),...%xStr,...
            'ylim',[0.5,1],'yTick',[0.5:0.1:1]);
        set(gca,'FontSize',15)
        legend('Correct');
        xl=xlabel('Date'); set(xl,'fontsize',25);
        yl=ylabel('Proportion "Correct"'); set(yl,'fontsize',25);
        b = glmfit(1:n, allCorrect);
        hold on;
        plot(b(1) + b(2)*(1:n),'k','LineWidth',2)
        if(idealFlag);
           %cumCorr = cumsum(allCorrect_i);
           %plot(1:n,ones(1,n)*mean(allCorrect_i),'k:','linewidth',2); 
           plot(1:n,ones(1,n)*0.79);
        end
    end

    if(choiceFlag) 
        f2 = figure('Position', [200, 500, 1300, 500]);
        plot(1:n,allChoiceProb(:,1),'r-o',1:n,allChoiceProb(:,2),'b-o',...
            1:n,allChoiceProb(:,3),'g-o',1:n,allChoiceProb(:,4),'m-o');
        hold on;
        if(dispTTflag)
            for i = 1:n
                text(i,0.95,sprintf('%i',1000*trialTimeAll(i)));
            end
        end
        plot([0,n],[0.5,0.5],'k');
        set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),'XtickLabel',xStr,...
            'ylim',[0.4,1])
        legend('Color','Shape','Border','Direction');
        xlabel('Trial Time (ms)');
        title('Cue Choice Preference');
        
        f2b = figure('Position', [816 179 1016 727]);
        barLabels = {'Blue','Red';'Circle','Square';'White','Black';'Horizontal','Vertical'};
        barTitles = {'Color','Shape','Border','Orienation'};
        for i=1:4 % index through cues
            subplot(2,2,i), hold on
            bar([mean(allChoiceProb(:,i));mean(1-allChoiceProb(:,i))])
            errorbar([mean(allChoiceProb(:,i));1-mean(allChoiceProb(:,i))],...
                [std(allChoiceProb(:,i));std(allChoiceProb(:,i))],'.');
            ylim([0 1]);
            set(gca, 'XTick', [1 2], 'XTickLabel', barLabels(i,:))
            title(barTitles{i});
            if(i == 1)
            ylabel('Proportion of Times Picked')
            end
            set(gca,'fontsize',13,'linewidth',2);
        end
        text(2,0.95,sprintf('n = %d',size(allChoiceProb,1)),'fontsize',14);
        
    end

    if(logRegFlag)
        f3 = figure('Position', [200, 500, 1300, 500]);
        plot(1:n,allThresholds,'k.','markersize',20); hold on;
        if(dispTTflag)
            for i = 1:n
                text(i,1.05*allThresholds(i),sprintf('%i',1000*trialTimeAll(i)));
            end
        end
%         set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),'XtickLabel',xStr,...
%             'ylim',[-0.3,1.1*max(allThresholds)]);
%         xlabel('Day');
%         title('Thresholds');
%         
         set(gca,'xlim',[0,n+1],'Xtick',10:10:n,'XtickLabel',[10:10:n],...
            'ylim',[-0.3,1.1*max(allThresholds)],'fontsize',15);
        xl=xlabel('Day'); set(xl,'fontsize',25);
        yl= ylabel('Thresholds');set(yl,'fontsize',25);
        
        b = glmfit(1:n, allThresholds);
        hold on;
        plot(b(1) + b(2)*(1:n),'k','LineWidth',2)
        Rsq = 1 - sum((allThresholds - (b(1) + b(2)*(1:n))).^2)/sum((allThresholds - mean(allThresholds)).^2);
       % t = text(10,0.4,sprintf('%.2f',Rsq));
        ylim([0,0.6]);
%         f=fit((2:n)',allThresholds(2:end)','power1');
%         hold on; plot(1:n,f.a*(1:n).^(f.b),'k-','linewidth',2);
%         Rsq = 1 - sum((allThresholds(2:end) - f.a*(2:n).^(f.b)).^2)/sum((allThresholds(2:end) - mean(allThresholds(2:end))).^2);
%         t = text(10,1.4,sprintf('%.2f',Rsq));
%         set(t,'fontsize',15);
        if(idealFlag)
           plot(1:n,ones(1,n)*(mean(allThresholds_i)),'k:','linewidth',2);         
        end
    end

    if(subWeightFlag)
        f4 = figure('Position', [200, 500, 1300, 500]);
%         for i = 1:n
%             subweightAll(i,:) = subweightAll(i,:)./(max(subweightAll(i,:)));
%         end
        plot(1:n,subweightAll(:,4),'r-o',1:n,subweightAll(:,3),'b-o',...
            1:n,subweightAll(:,2),'g-o',1:n,subweightAll(:,1),'m-o');
        if(dispTTflag)
            for i = 1:n
                text(i,1.05*max(subweightAll(:)),sprintf('%i',1000*trialTimeAll(i)));
            end
        end

        %plot([0,n],[0.8,0.8],'b:',[0,n],[0.6,0.6],'r:',[0,n],[0.4,0.4],'g:',...
         %   [0,n],[0.2,0.2],'m:','LineWidth',2);
        set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),'XtickLabel',xStr,...
            'ylim',[-0.3,1.1*max(subweightAll(:))]);
        legend('Color', 'Shape', 'Border','Direction');
        xlabel('Trial Time (ms)');
        title('Subjective Cue Weight');
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
%         plot(1:n,bf(:,5),'k-*',1:n,bf(:,11),'k-o');
%         if(dispTTflag)
%             for i = 1:n
%                 text(i,3.1,sprintf('%i',1000*trialTimeAll(i)));
%                 [~,maxbI(i)] = max(bf(i,:));
% 
%             end
%         end
%         line([0 n+1],[3 3],'linewidth',2); hold on;
%          set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),'XtickLabel',xStr,...%num2cell(1000*trialTimeAll),...
%         'ylim',[-0.3,1.1*max(subweightAll(:))]);
%         legend('Model 5','Model 11');
%         title('Bayes Model Comparison');
%         ylabel('Bayes Factor (normalized to 15)');
%         ylim([0,15]);
%         subplot(2,1,2);
%         hist(maxbI,1:15);
%         set(gca,'xlim',[0,15],'Xtick',1:15);
%         title('Max Bayes Model');

        for bi = 1:n
            bfcount(bi,:) = bf(bi,:)>3;
            if(sum(bfcount(bi,:))==0)
                bfcount(bi,15) = 1;
            elseif(sum(bfcount(bi,:))>1)
                
               bfcount(bi,:) = bf(bi,:)== max(bf(bi,:));
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
        plot(1:n,rhoAll(:,5),'k-*',1:n,rhoAll(:,11),'k-o',1:n,rhoAll(:,15),'r-o');
        if(dispTTflag)
            for i = 1:n
                text(i,0.95,sprintf('%i',1000*trialTimeAll(i)));
                [~,maxrI(i)] = max(rhoAll(i,:));
            end
        end
         set(gca,'xlim',[0,n+1],'Xtick',1:length(allDates),'XtickLabel',xStr,...%num2cell(1000*trialTimeAll),...
        'ylim',[0,1]);
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
    for i = 1:length(ignoreList)
        disp(ignoreList{i});
    end
end
