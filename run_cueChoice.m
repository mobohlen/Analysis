function [choiceCount,nChoice] = run_cueChoice(cueA,cueB,respStimulus,cueProbInfo)
% Gets cue choice counts for each dimension for each trial
%%% Inputs
% cueA, cueB - list of stimulus numbers from 1-16, with each entry
%              corresponding to the two presented stimuli on a given trial
% respStimulus - the stumlus number 1-16 that was selected for each trial
% cueProbInfo - contains information about the associated cue weights and
%               probabilities for the given set of trials (cueProbInfo.mat)
%%% Outputs
% choiceCount: Choices of cue option 1, when that dimension differs.
%       choiceCount(1) = Most informative dim (Arya - red over blue, Lana -
%                        vertical over horizontal)
%       choiceCount(2) 
%       choiceCount(3) 
%       choiceCount(4) = Least informative dim (Arya - horizontal over
%                       vertical, Lana - Blue over red)
%        
% nChoice : totals for each cue choice

stimprop = 2 - cueProbInfo.stimMat; % (convert 2 -> 0, 1 -> 1)
nValid = length(respStimulus);
nChoice = zeros(1,4);
choiceCount = zeros(1,4); 
for i = 1:nValid
    if(stimprop(cueA(i),1)~=stimprop(cueB(i),1))
        nChoice(1) = nChoice(1) +1;
        choiceCount(1) = choiceCount(1) + stimprop(respStimulus(i),1);
    end
    if(stimprop(cueA(i),2)~=stimprop(cueB(i),2))
        nChoice(2) = nChoice(2) +1;
        choiceCount(2) = choiceCount(2) + stimprop(respStimulus(i),2);
    end
    if(stimprop(cueA(i),3)~=stimprop(cueB(i),3)) 
        nChoice(3) = nChoice(3) +1;
        choiceCount(3) = choiceCount(3) + stimprop(respStimulus(i),3);
    end
    if(stimprop(cueA(i),4)~=stimprop(cueB(i),4)) 
        nChoice(4) = nChoice(4) +1;
        choiceCount(4) = choiceCount(4) + stimprop(respStimulus(i),4);
    end
end

end

