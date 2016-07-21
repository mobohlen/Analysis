function [rho,p] = run_modelCorr(stimWeights,cueA,cueB,response)
% clear; clc
% close all
% load('../Data/Arya/Aryasummary_08_Sep_2015fourCueDim_1.mat')
% indValid= find(summary.respLocation ~= -2 &...
%                    summary.respLocation ~= -1 &...
%                    summary.respLocation ~= 0)'; % indeces of valid trials, where either cue was hit
% cueA = summary.cueA(indValid);
% cueB = summary.cueB(indValid);
% response = summary.respStimulus(indValid);
% stimWeights = [0.9,0.8,0.7,0.6];

% cue dimension properties 
% col 1: blue = 1, red = 0
% col 2: circle = 1, square = 0
% col 3: white = 1, black = 0
% col 4: horiz = 1, verti = 0
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

weightedStims(:,1) = stimWeights(1).*stimprop(:,1) + (1-stimWeights(1)).*(1-stimprop(:,1));
weightedStims(:,2) = stimWeights(2).*stimprop(:,2) + (1-stimWeights(2)).*(1-stimprop(:,2));
weightedStims(:,3) = stimWeights(3).*stimprop(:,3) + (1-stimWeights(3)).*(1-stimprop(:,3));
weightedStims(:,4) = stimWeights(4).*stimprop(:,4) + (1-stimWeights(4)).*(1-stimprop(:,4));

Model(1,:) = [1 0 0 0]; % 1
Model(2,:) = [0 1 0 0]; % 2
Model(3,:) = [0 0 1 0]; % 3
Model(4,:) = [0 0 0 1]; % 4
Model(5,:) = [1 1 0 0]; % 5
Model(6,:) = [1 0 1 0]; % 6
Model(7,:) = [1 0 0 1]; % 7
Model(8,:) = [0 1 1 0]; % 8
Model(9,:) = [0 1 0 1]; % 9
Model(10,:) = [0 0 1 1]; % 10
Model(11,:) = [1 1 1 0]; % 11
Model(12,:) = [1 1 0 1]; % 12
Model(13,:) = [1 0 1 1]; % 13
Model(14,:) = [0 1 1 1]; % 14
Model(15,:) = [1 1 1 1]; % 15


for m = 1:length(Model)
    subWeightA = sum(repmat(Model(m,:),length(cueA),1) .* weightedStims(cueA,:),2);
    subWeightB = sum(repmat(Model(m,:),length(cueA),1) .* weightedStims(cueB,:),2);
    
    idxDiffW = subWeightA ~= subWeightB; % is it a problem that weights are much 
    idxSameW = subWeightA == subWeightB;
    % more likely to be the same for model 1 vs model 4????????
    dec = zeros(length(cueA),1);
    dec(idxDiffW) = subWeightA(idxDiffW) > subWeightB(idxDiffW);
    dec(idxSameW) = 0.5;
    monkeyDec = (response == cueA);
    
    [rho(m),p(m)] = corr(dec,monkeyDec);
end
