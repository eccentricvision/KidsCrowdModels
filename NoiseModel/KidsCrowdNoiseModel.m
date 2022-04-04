%KidsCrowdNoiseModel
%
%generates the output of a population noise model of crowding
%target responses are generated via a population of orientation-selective detectors
%and corrupted by both early noise and late noise when flankers are
%present (but no population responses to the flankers are generated)
%formerly OrientCrowdPoolNoiseModelPlotGrp
%Plots one group against model output
%v2.1 - Uses v2.1 of OrientCrowdNoiseModelFun
%
%Alexandra Kalpadakis-Smith & John Greenwood
%April 2017
%modified J Greenwood September 2021

clear all;
modelver = 2.1;

%For saving params later
thisFile= 'KidsCrowdNoiseModel.m';
ThisDirectory=which(thisFile); ThisDirectory=ThisDirectory(1:end-length(thisFile));

RootFolder    = ThisDirectory;%strcat(MatlabRoot,'Models/PacManKidsPopulationModel/PoolModel/'); %where the code and parameters can be found (in subdirs)
DataFolder    = strcat(ThisDirectory(1:end-11),'GroupData/'); %where the data files are (separate from Pool & Noise Model folders)

%% select subject and group

%which group
for WhichGroup = 1:6%WhichGroup=1;%DefInput('Which group do you want to model? 1=amblyope kids, 2=control kids, 3-6=adult periphery', 1); %0 = controls, 1=amblyopes

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
        grplab={'adu2.5'};
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
    disp(grplab);

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
    SmoothData    = 1;
    NumTrials     = [1000 1000 1000];

    %model parameters: free
    [param,boot,error,FitType] = KidsCrowdNoiseModel_LookupParams(WhichGroup); %load current best-fitting parameters
    disp(param.EarlyNoise); disp(param.CrowdNoise);
    param.NumFreeP = 2; %if ModelSD is fixed then Noise model has 2 free parameters

    %% Convert to frequency of response data

    for po =1:3 %T-F orient diff
        CompData(po,:) = CompData(po,:)./sum(CompData(po,:));
    end

    %% run the model function

    if ~isfield(boot,'AllResps') %if the bootstrapping has not been run - just plot a single instantiation of the model
        clear error; %need to delete the NaN entry loaded above if not actually run
        [DetectorSen,TargetPopResp,ComboPopResp,ModelDecision,histdest]=KidsCrowdNoiseModelFun(expt,param.ModelSD,param.EarlyNoise,param.CrowdNoise); %,LateNoise(lns)
        error.LSE    = sum(All(CompData-histdest).^2);
        error.LSEsum = sum(error.LSE);

        [error.AIC,error.AICc] = ComputeAIC(error.LSEsum,NumDataPoints*expt.NumCrowdConds,param.NumFreeP);
        disp('Single model run');
    else %do nothing - the model, bootstrapping and AIC have already been run
        disp('Bootstrapped runs - mean used for AIC')
    end

    disp(' ');disp(strcat('Least Squared Error: ',num2str(error.LSEsum)));
    disp(strcat('AIC: ',num2str(error.AIC)));disp(' ');

    %% plot

    if WhichGroup<3 %kids
        ymaxval = 0.25;
    else
        ymaxval = 0.4;
    end

    if ~isfield(boot,'AllResps') %if the bootstrapping has not been run - just plot a single instantiation of the model

        %Histogram of Combined Response
        figure
        for cc=1:expt.NumCrowdConds
            subplot(1,expt.NumCrowdConds,cc)
            h(1) = plot(expt.binmidpoint,histdest(cc,:),'b-');
            hold on;
            h(2) = plot(expt.binmidpoint, CompData(cc,:), 'ro');
            title(strcat(grplab,' group, ',CondLabels{cc}));
            if cc==3
                legend(h, {'Noise Model', 'Data'});
            end
            axis square;
            xlim([-180 180]);
            xtick(-180:45:180);
            ylim([0 ymaxval]);%max(All([histdest CompData(:,:)]))+0.05]);
        end
    else %bootstrapping run - plot full distributions

        figure
        for cc=1:expt.NumCrowdConds
            subplot(1,3,cc)
            h(1)=plot(expt.binmidpoint, CompData(cc,:),'o');
            hold on;
            h(2)=plot(expt.binmidpoint,boot.HistMean(cc,:),'k-');
            plot(expt.binmidpoint, boot.UpperRange(cc,:),'k--');
            plot(expt.binmidpoint, boot.LowerRange(cc,:),'k--');
            xlim([-180 180]);
            xtick(-180:45:180);
            ylim([0 ymaxval]);%max(All([histdest CompData(:,:)]))+0.05]);
            plot([0 0],[0 max(All([boot.HistMean(:,:) CompData]))+0.05],'r--'); %target location
            if cc>1
                plot([expt.FlankerOrient(cc) expt.FlankerOrient(cc)],[0 max(All([boot.HistMean(:,:) CompData]))+0.05],'r--'); %flanker location
            end
            axis square;
            box off;
            legend(h,{'Data', 'Noise Model'});
            title(strcat(grplab,' group, ',CondLabels{cc}));

        end
    end
end
