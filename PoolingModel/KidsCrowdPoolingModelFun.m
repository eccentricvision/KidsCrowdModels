function [DetectorSen,TargetPopResp,FlankerPopResp,ComboPopResp,ModelDecision,histdest]=KidsCrowdPoolingModelFun(expt,ModelSD,EarlyNoise,CrowdNoise,FlankerWeight_30,FlankerWeight_90)

%[DetectorSen,TargetPopResp,FlankerPopResp,ComboPopResp,ModelDecision,histdest]=KidsCrowdPoolingModelFun(expt,ModelSD,NoiseMax,CrowdNoise, FlankWeight_30, FlankWeight_90,FitParams)
%v2.1 - crowdnoise now added to flanker population only, flipped to model flanker weights instead of target weight
%new in v2.1 is variable bin widths through the 'expt' structure
%
%code by Alexandra Kalpadakis Smith
%modified J Greenwood September 2021

%orientation / experiment parameters - mostly given by the expt structure now
% OrientAxis    = -180:1:180;
% TargetOrient  = 0;
% FlankerOrient = [NaN 30 90]; %uncrowded / crowded-30 / crowded-90
% binmidpoint   = -180:10:180;
% NumCrowdConds = 3; %uncrowded / crowded-30 / crowded-90
NumTrials     = [1000 1000 1000];%[120 240 240]*20;

%model parameters: fixed
ModelPeakOrientDiff = 1;
OrientAxis          = -180-(round(expt.BinWidth/2)-1):ModelPeakOrientDiff:180+(round(expt.BinWidth/2)-1); %setup additional detectors outside range to avoid wrapping issues
ModelPeakVals       = OrientAxis;
NumDetectors        = numel(ModelPeakVals);
PeakVal             = 1;
BaseVal             = 0;

FlankerWeightCond   = [0 FlankerWeight_30 FlankerWeight_90];    %uncrowded / crowded-30 / crowded-90
TargetWeightCond    = 1-FlankerWeightCond;%[1 TargetWeight_30 TargetWeight_90]; %uncrowded / crowded-30 / crowded-90

%% generate the population

for dd=1:NumDetectors
    DetectorSen(dd,:) = DrawGaussian(OrientAxis,ModelPeakVals(dd),ModelSD,PeakVal-BaseVal,BaseVal);
end

%% generate a response to stimuli

TargetInd = find(OrientAxis==expt.TargetOrient);

for cc=1:expt.NumCrowdConds
    TargetPopResp{cc} = repmat((DetectorSen(:,TargetInd)'),[NumTrials(cc) 1])+(randn(NumTrials(cc),NumDetectors)*EarlyNoise);
    if cc>1
        FlankerInd = find(OrientAxis==expt.FlankerOrient(cc));
        FlankerPopResp{cc} = repmat((DetectorSen(:,FlankerInd)'),[NumTrials(cc) 1])+(randn(NumTrials(cc),NumDetectors)*CrowdNoise);
    else %uncrowded
        FlankerPopResp{cc} = TargetPopResp{cc};
    end
end

%% make the crowding happen

for cc=1:expt.NumCrowdConds
    TrialTargetWeights{cc} = (TargetWeightCond(cc)*ones(1,NumTrials(cc)));
    TrialFlankerWeights{cc} = (FlankerWeightCond(cc)*ones(1,NumTrials(cc)));
    
    if cc>1
        for tt=1:NumTrials(cc)
            ComboPopResp{cc}(tt,:) = ((TargetPopResp{cc}(tt,:).*TrialTargetWeights{cc}(tt))+(FlankerPopResp{cc}(tt,:).*TrialFlankerWeights{cc}(tt))); 
        end
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
