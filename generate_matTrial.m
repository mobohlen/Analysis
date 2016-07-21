function matTrial = generate_matTrial(cueA,cueB,respStimulus,correct,cueProbInfo)
% Creates matrix used in calculating several perfomance metrics including
% subjective weights, psychometric curves, and Bayes model analysis
%%% Inputs
% cueA, cueB - list of stimulus numbers from 1-16, with each entry
%              corresponding to the two presented stimuli on a given trial
% respStimulus - the stumlus number 1-16 that was selected for each trial
% correct - Whether or not the given stimulus was the higher weighted
%           stimulus fore each trial
% cueProbInfo - contains information about the associated cue weights and
%               probabilities for the given set of trials (cueProbInfo.mat)
%%% Outpus
% matTrial - rows: trial
%            cols: 1-4 -> difference in weights for each possible
%                         combination of cue dimensions
%                  5 -> sum(cols 1-4)
%                  6 -> Conditional probability of A being correct given
%                       both stimuli present, where
%                       Pr(A | A&B) = 10^[col5]/(1+10^[col5])

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
end

