function [DetectorSen,TargetPopResp,ComboPopResp,ModelDecision,histdest]=KidsCrowdNoiseModelFun(expt,ModelSD,EarlyNoise,CrowdNoise)
%[DetectorSen,TargetPopResp,CrowdPopResp,ModelDecision,histdest]=KidsCrowdNoiseModelFun(ModelSD,EarlyNoise,CrowdNoise)
%v2.1 - some tweaks by JG to rename params and remove unused lines
%new in v2.1 is variable bin widths through the 'expt' structure
%
%code by Alexandra Kalpadakis Smith
%modified J Greenwood September 2021

%orientation / experiment parameters
% OrientAxis    = -180:1:180;
% TargetOrient  = 0;
% FlankerOrient = [NaN 30 90]; %uncrowded / crowded-30 / crowded-90
% binmidpoint   = -180:10:180;
NumTrials     = [1000 1000 1000];

%model parameters: fixed
ModelPeakOrientDiff = 1;
OrientAxis          = -180-(round(expt.BinWidth/2)-1):ModelPeakOrientDiff:180+(round(expt.BinWidth/2)-1); %setup additional detectors outside range to avoid wrapping issues
ModelPeakVals       = OrientAxis;
NumDetectors        = numel(ModelPeakVals);
PeakVal             = 1;
BaseVal             = 0;

%% generate the population

for dd=1:NumDetectors
    DetectorSen(dd,:) = DrawGaussian(OrientAxis,ModelPeakVals(dd),ModelSD,PeakVal-BaseVal,BaseVal);
end

%% generate a response to stimuli and make the crowding happen (if it does)

TargetInd = find(OrientAxis==expt.TargetOrient);

for cc=1:expt.NumCrowdConds
    TargetPopResp{cc} = repmat((DetectorSen(:,TargetInd)'),[NumTrials(cc) 1])+(randn(NumTrials(cc),NumDetectors)*EarlyNoise);
    if cc>1
        ComboPopResp{cc} = TargetPopResp{cc}+(randn(NumTrials(cc),NumDetectors)*CrowdNoise); %add crowding noise if it's a crowded trial
    else %uncrowded
        ComboPopResp{cc} = TargetPopResp{cc};
    end
end

%% pull out the orientation decisions and plot the distribution

for cc=1:expt.NumCrowdConds
    [~,MaxInd] = max(ComboPopResp{cc},[],2);
    ModelDecision{cc}   = OrientAxis(MaxInd);
    
    [histdest(cc,:)] = hist(ModelDecision{cc},expt.binmidpoint);
    histdest(cc,:) = histdest(cc,:)/sum(histdest(cc,:));
end

end
