function [bNS,bS] = run_subjectiveWeights(matTrial,locMat)
% Calculates subjective weights for a set of trials
%%% Inputs: 
% matTrial - generated from generate_matTrial.m
% locMat - col 1: cueALoc
%          col 2: cueBLoc
%%% Outputs:
% bNS - subjective weights with no spatial component
% bS - subjective weights with spatial component included
%      organized:
%       [cue dim 1 (lowest), cue dim 2, cue dim 3, cue dim 4 (highest),...
%       constant bias, L/R dim, U/D dim]

nValid = size(matTrial,1);
cueALoc = locMat(:,1); cueBLoc = locMat(:,2);
%% create regressors for individual cues 
regX = zeros(nValid,4);
regX(matTrial(:,1:4) > 0) = 1;
regX(matTrial(:,1:4) < 0) = -1;
regY = matTrial(:,7);

for iT = 1:size(matTrial,1);
    % stimulus positions 1-4; 1 = leftUp, 2 = leftDown, 3 = rightUp, 4 = rightDown
    % Col 1 = L/R     1 = right
    locX(iT,1) = ( 0*((cueALoc(iT) == 1) || (cueALoc(iT) == 2)) +...
                     1*((cueALoc(iT) == 3) || (cueALoc(iT) == 4)) ) -...  
                  ( 0*((cueBLoc(iT) == 1) || (cueBLoc(iT) == 2)) +...
                     1*((cueBLoc(iT) == 3) || (cueBLoc(iT) == 4)) );
    % Col 2 = U/D     1 = down         
    locX(iT,2) = ( 1*((cueALoc(iT) == 2) || (cueALoc(iT) == 4)) +...
                     0*((cueALoc(iT) == 1) || (cueALoc(iT) == 3)) ) -...  
                  ( 1*((cueBLoc(iT) == 2) || (cueBLoc(iT) == 4)) +...
                     0*((cueBLoc(iT) == 1) || (cueBLoc(iT) == 3)) );
end

%% No spatial component
[b1,dev,stats] = glmfit(regX,[regY ones(nValid,1)],'binomial','link','logit');
bNS = log10(exp(b1)); % transform beta weights to log base 10
    
%% With spatial component
[b2,dev2,stats2] = glmfit([regX,locX],[regY ones(nValid,1)],'binomial','link','logit');
bS = log10(exp(b2)); % transform beta weights to log base 10

end

