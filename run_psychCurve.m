function [sumW,zSumW,threshold] = run_psychCurve(matTrial)
% Calculates psychometric curve and 50-75 threshold for a set of trials
%%% Inputs: 
% matTrial - generated from generate_matTrial.m
%%% Outputs:
% sumW - col 1: number of A choices
%        col 2: total number of observations
%        col 3: sum of weights, y
% sumZ - Regression output
%
% threshold - estimated difference in weight needed to increase proportion
%             of cue A choices from 50% to 75%
% bS - subjective weights with spatial component included
%      organized:
%       [cue dim 1 (lowest), cue dim 2, cue dim 3, cue dim 4 (highest),...
%       constant bias, L/R dim, U/D dim]

nValid = size(matTrial,1);
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