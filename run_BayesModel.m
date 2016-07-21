function BayesFactor = run_BayesModel(matTrial)
% Calculates subjective weights for a set of trials
%%% Inputs: 
% matTrial - generated from generate_matTrial.m
%%% Outputs:
% BayesFactor - calculated Bayes factor for each model relatvie to model 15

nValid = size(matTrial,1);
% create regressors for individual cues (this is for subjective cue weights)
regX = zeros(nValid,4);
regX(matTrial(:,1:4) > 0) = 1;
regX(matTrial(:,1:4) < 0) = -1;

%% Model construction
nC = 1; 
for iN = 1 : 4 % all possible combinations
    indCue = nchoosek([1 2 3 4],iN);
    for iC = 1 : size(indCue,1)
        indCueModel{nC} = indCue(iC,:);
        nC = nC + 1;
    end
end

%% model comparsion based on variational bayes
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
        xModel{iM} = zeros(size(regX));
        for iT = 1 : nValid
            temp1 = find(regX(iT,:) ~= 0, 1);
            xModel{iM}(iT,temp1) = regX(iT,temp1);
        end
    end
    % run Bayesian logistic regression
    % see Drusgowitsch(2013)
    [wB, V, invV, logdetV, E_a, LB] = bayes_logit_fit([ones(nValid,1) xModel{iM}], yModel);
    L(iM) = LB; % log-liklihood of the data given the model (lower bound); larger L = better fit
    clear wB V invV logdetV E_a LB
end
adjL15 = L - L(15); % since L is relative, adjust L based on model 15 (this gives you log Bayes factor)
BayesFactor = exp(adjL15); % a model with Bayes factor > 3 usually means that this model significantly outperforms the optimal model (model 15)


end

