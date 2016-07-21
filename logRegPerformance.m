function [sumW,zSumW,threshold,bLog10,BayesFactor,L] = logRegPerformance(cues,respStimulus,correct,cueProbInfo)
cueA = cues(:,1); cueB = cues(:,2); 
probFeedback = cueProbInfo.probFeedback{cueProbInfo.iC};
stimComb = cueProbInfo.stimComb;
nValid = length(respStimulus);
matTrial = zeros(nValid,8);
for iT = 1 : nValid
    numStimComb = find(ismember(stimComb,[cueA(iT) cueB(iT)],'rows')); % identify the stimulus combination presented
    matTrial(iT,1:6) = probFeedback(numStimComb,:); % associated probability for a given combination
    % determine the cue chosen (A = 1, B = 0)
    if respStimulus(iT) == cueA(iT)
        matTrial(iT,7) = 1;
    end
    matTrial(iT,8) = correct(iT); % flag for response
end
% create regressors for sum of cues (this is for regression curve)
sumWeights = unique(round(matTrial(:,5)*10)); % determine unique sums (integers are easy to handle)
sumW = zeros(length(sumWeights),3);
for iS = 1 : length(sumWeights)
    tempS = find(round(matTrial(:,5)*10) == sumWeights(iS));
    sumW(iS,1) = sum(matTrial(tempS,7) == 1); % number of A choices
    sumW(iS,2) = length(tempS); % total number of observations
    sumW(iS,3) = sumWeights(iS)/10; % sum of weights, y
end

% create regressors for individual cues (this is for subjective cue weights)
regX = zeros(nValid,4);
regX(matTrial(:,1:4) > 0) = 1;
regX(matTrial(:,1:4) < 0) = -1;
regY = matTrial(:,7);

%% run regression: curve based on sum(cueL - cueR)
[b1,dev1,stats1] = glmfit(sumW(:,3),sumW(:,1:2),'binomial','link','logit');
zSumW = 1./(1+exp(-(b1(1)+sumW(:,3).*b1(2))));

xVec = linspace(-2,2, 10000);
zSumWVec = 1./(1+exp(-(b1(1)+xVec.*b1(2))));
x50 = xVec(find(zSumWVec > 00.5, 1,'first'));
x75 = xVec(find(zSumWVec > 0.75, 1,'first'));
if ~isempty(x75 - x50)
    threshold = x75 - x50;
    if threshold < 0.001 
        threshold = nan;
    end
else
    threshold = nan;
end

%% run regression: estimating subjective cue weights
[b2,dev2,stats2] = glmfit(regX,[regY ones(nValid,1)],'binomial','link','logit');
bLog10 = log10(exp(b2)); % transform beta weights to log base 10
% bLog10 = [bias; cue0.8; cue0.6; cue0.4; cue0.2]

%% Strategy model comparison
% model construction
nC = 1; 
for iN = 1 : 4 % all possible combinations
    indCue = nchoosek([1 2 3 4],iN);
    for iC = 1 : size(indCue,1)
        indCueModel{nC} = indCue(iC,:);
        nC = nC + 1;
    end
end

% model comparsion based on variational bayes
% for this analysis, response should be defined as cueA = 1; cueB = -1
yModel = matTrial(:,7);
yModel(matTrial(:,7) == 0) = -1;

w = zeros(size(indCueModel,2)+1,5);
wLog10 = zeros(size(indCueModel,2)+1,5);
diagV = zeros(size(indCueModel,2)+1,5);
for iM = 1 : size(indCueModel,2)+1
    % create input matrix
    if(iM<=15)
    xModel{iM} = regX(:,indCueModel{iM});
    else
%         for i = 1:length(regX)
%             for j = 1:4
%                 if(regX(i,j)~=0)
%                 xModel{16}(i) = regX(i,j);
%                 break;
%                 end
%             end
%         end
%         xModel{16} = xModel{16}';
        xModel{iM} = zeros(size(regX));
        for iT = 1 : nValid
            temp1 = find(regX(iT,:) ~= 0, 1);
            xModel{iM}(iT,temp1) = regX(iT,temp1);
        end
    end
    % run Bayesian logistic regression
    % see Drusgowitsch(2013)
    [wB, V, invV, logdetV, E_a, LB] = bayes_logit_fit([ones(nValid,1) xModel{iM}], yModel);
%     w(iM,[1 indCueModel{iM}+1]) = wB; % cue weights
%     w10(iM,[1 indCueModel{iM}+1]) = log10(exp(wB)); % cue weights with log base 10
%     diagV(iM,[1 indCueModel{iM}+1]) = diag(V); % variance of w
    L(iM) = LB; % log-liklihood of the data given the model (lower bound); larger L = better fit
    clear wB V invV logdetV E_a LB
end
adjL15 = L - L(15); % since L is relative, adjust L based on model 15 (this gives you log Bayes factor)
BayesFactor = exp(adjL15); % a model with Bayes factor > 3 usually means that this model significantly outperforms the optimal model (model 15)
