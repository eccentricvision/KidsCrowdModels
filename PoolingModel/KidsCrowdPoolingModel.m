%KidsCrowdPoolingModel
%
%generates the output of a population pooling model of crowding
%target and flanker responses are generated via a population of
%orientation-selective detectors and combined/pooled with weights to
%determine the relative magnitude of each population response
%formerly OrientCrowdPoolNoiseModelPlotGrp
%v2.1 - Uses v2.1 of OrientCrowdPopModelFun 
%
%Alexandra Kalpadakis-Smith & John Greenwood
%April 2017
%modified J Greenwood September 2021

clear all;
modelver = 2.1;

%For saving params later
thisFile = 'KidsCrowdPoolingModel.m';
ThisDirectory=which(thisFile); ThisDirectory=ThisDirectory(1:end-length(thisFile));

RootFolder    = ThisDirectory;%strcat(MatlabRoot,'Models/PacManKidsPopulationModel/PoolModel/'); %where the code and parameters can be found (in subdirs)
DataFolder    = strcat(ThisDirectory(1:end-13),'GroupData/'); %where the data files are (separate from Pool & Noise Model folders)

%% select subject and group

%which group
for WhichGroup=1:6 %WhichGroup = 1;%DefInput('Which group do you want to model? 1=amblyope kids, 2=control kids, 3-6=adult periphery', 1); %0 = controls, 1=amblyopes
    
    WhichFit  = 3; %use bootstrapped values (3) if possible, or fine fit (2) or coarse (1)
    
    %Select and load data
    if WhichGroup == 1
        load(strcat(DataFolder,'AmbRawGrpData.mat')); %load('/Users/alex/Documents/MATLAB/PacMan2/PacAlex/PoolModel/Data/AmbRawIndData.mat');
        %CompData = histdest;%squeeze(AmbRawIndData);
        grplab={'amb'};
    elseif WhichGroup ==2
        load(strcat(DataFolder,'ConRawGrpData.mat'));
        %CompData = squeeze(ConRawIndData);
        grplab={'con'};
    elseif WhichGroup == 3
        load(strcat(DataFolder,'AdultRawGrpData_2.mat'));
        %CompData = histdest;
        grplab={'adult2.5'};
    elseif WhichGroup ==4
        load(strcat(DataFolder,'AdultRawGrpData_5.mat'));
        %CompData = histdest;
        grplab={'adult5'};
    elseif WhichGroup == 5
        load(strcat(DataFolder,'AdultRawGrpData_10.mat'));
        %CompData = histdest;
        grplab={'adult10'};
    elseif WhichGroup == 6
        load(strcat(DataFolder,'AdultRawGrpData_15.mat'));
        %CompData = histdest;
        grplab={'adult15'};
    end
    CompData = histdest; %same for all groups
    
    CondLabels = {'Unflanked','30deg flankers','90deg flankers'};
    
    %% get parameters and plot
    
    %orientation / experiment parameters
    expt.OrientAxis    = -180:1:180;
    expt.TargetOrient  = 0;
    expt.FlankerOrient = [NaN 30 90]; %uncrowded / crowded-30 / crowded-90
    expt.BinWidth      = 10;
    expt.binmidpoint   = -180:expt.BinWidth:180;
    expt.NumCrowdConds = 3; %uncrowded / crowded-30 / crowded-90
    
    NumDataPoints = numel(expt.binmidpoint);
    SmoothData    = 0;
    NumTrials     = [1000 1000 1000];
    
    %model parameters: free
    [param,boot,error,FitType] = KidsCrowdPoolingModel_LookupParams(WhichGroup); %load current best-fitting parameters
    disp(grplab); disp(param.EarlyNoise); disp(param.CrowdNoise); disp(param.FlankerWeight_30); disp(param.FlankerWeight_90);
    
    param.NumFreeP = 4; %if ModelSD is fixed then PoolNoise model has 4 free parameters
    
    %% Convert to frequency of response data
    
    for po =1:3 %T-F orient diff
        CompData(po,:) = CompData(po,:)./sum(CompData(po,:));
    end
    
    %% run the model function
    
    if ~isfield(boot,'AllResps') %if the bootstrapping has not been run - just plot a single instantiation of the model
        clear error; %need to delete the NaN entry loaded above if not actually run
        [DetectorSen,TargetPopResp,FlankerPopResp,ComboPopResp,ModelDecision,histdest]=KidsCrowdPoolingModelFun(expt,param.ModelSD,param.EarlyNoise, param.CrowdNoise, param.FlankerWeight_30,param.FlankerWeight_90); %,LateNoise(lns)
        
        error.LSE    = sum(All(CompData-histdest).^2);
        error.LSEsum = sum(error.LSE);
        
        [error.AIC,error.AICc] = ComputeAIC(error.LSEsum,NumDataPoints*expt.NumCrowdConds,param.NumFreeP);
        disp('Single model run');
    else %model, bootstrapping and AIC have already been run, but get a model output to plot example distributions (ignoring histdest)
        [DetectorSen,TargetPopResp,FlankerPopResp,ComboPopResp,ModelDecision,~]=OrientCrowdPoolNoiseModelFun(expt,param.ModelSD,param.EarlyNoise, param.CrowdNoise, param.FlankerWeight_30,param.FlankerWeight_90); %,LateNoise(lns)
        disp('Bootstrapped runs - mean used for AIC')
    end
    
    disp(' ');disp(strcat('Least Squared Error: ',num2str(error.LSEsum)));
    disp(strcat('AIC: ',num2str(error.AIC)));disp(' ');
    
    %get mean distributions
    for cc=1:expt.NumCrowdConds
        MeanTargetPopResp(cc,:)  = mean(TargetPopResp{cc},1);
        MeanFlankerPopResp(cc,:) = mean(FlankerPopResp{cc},1);
        MeanComboPopResp(cc,:)   = mean(ComboPopResp{cc},1);
    end
    ModelPeakOrientDiff = 1;
    ModelPeakVals       = -180-(round(expt.BinWidth/2)-1):ModelPeakOrientDiff:180+(round(expt.BinWidth/2)-1);
    
    %% plot
    if WhichGroup<3 %kids
        ymaxval = 0.25;
    else
        ymaxval = 0.4;
    end
    
    %plot responses
    if ~isfield(boot,'AllResps') %if the bootstrapping has not been run - just plot a single instantiation of the model
        %Histogram of Combined Response
        figure
        clear h;
        for cc=1:expt.NumCrowdConds
            subplot(1,expt.NumCrowdConds,cc)
            h(1) = plot(expt.binmidpoint,histdest(cc,:),'b-');
            hold on;
            h(2) = plot(expt.binmidpoint, CompData(cc,:), 'ro');
            title(strcat(grplab,' group, ',CondLabels{cc}));
            if cc==3
                legend(h, {'PoolingModel', 'Data'});
            end
            axis square;
            xlim([-180 180]);
            xtick(-180:45:180);
            ylim([0 ymaxval]);
            ytick(0:0.05:ymaxval);
            ax = gca; % current axes
            ax.TickDir = 'out';
        end
    else %bootstrapping run - plot full distributions
        
        figure
        clear h;
        for cc=1:expt.NumCrowdConds
            subplot(1,3,cc)
            h(1)=plot(expt.binmidpoint, CompData(cc,:),'o');
            hold on;
            h(2)=plot(expt.binmidpoint,boot.HistMean(cc,:),'k-');
            plot(expt.binmidpoint, boot.UpperRange(cc,:),'k--');
            plot(expt.binmidpoint, boot.LowerRange(cc,:),'k--');
            xlim([-180 180]);
            xtick(-180:45:180);
            ylim([0 ymaxval]);
            ytick(0:0.05:ymaxval);
            plot([0 0],[0 ymaxval],'r--'); %target location
            if cc>1
                plot([expt.FlankerOrient(cc) expt.FlankerOrient(cc)],[0 ymaxval],'r--'); %flanker location
            end
            axis square;
            box off;
            legend(h,{'Data', 'PoolingModel'});
            title(strcat(grplab,' group',CondLabels{cc}));
            ax = gca; % current axes
            ax.TickDir = 'out';
        end
    end
end
