function [paramOut,boot,error,FitType] = KidsCrowdNoiseModel_LookupParams(grp)
%[paramOut,boot,error,FitType] = KidsCrowdNoiseModel_LookupParams(grp
%function to load current best-fitting parameters for either individual or group
%input grp=1-6 where 1=amblyope kids, 2=control kids, 3-6=adult periphery,
%returns params as a structure, including param.ModelSD, param.EarlyNoise, param.CrowdNoise
%
%v2.1 J Greenwood September 2021

thisFile= 'KidsCrowdNoiseModel_LookupParams.m';
ThisDirectory=which(thisFile); ThisDirectory=ThisDirectory(1:end-length(thisFile));

switch grp
    case 1 %amblyopia kids
        grplab = 'amb';
    case 2 %controls
        grplab = 'con';
    case 3 %adult 2.5
        grplab = 'adult2.5';
    case 4 %adult 5
        grplab = 'adult5';
    case 5 %adult 10
        grplab = 'adult10';
    case 6 %adult 15
        grplab = 'adult15';
end

%set up empty variables to return if not yet done
paramOut = NaN;
boot     = NaN;
error    = NaN;
FitType  = NaN;
%bootstrapped fits
%load bootstrap values and AIC
if exist(strcat(ThisDirectory,'NoiseModelBootstrapValues/',grplab,'-group-NoiseFine-BootAIC_v2.1.mat'),'file')
    load(strcat(ThisDirectory,'NoiseModelBootstrapValues/',grplab,'-group-NoiseFine-BootAIC_v2.1.mat')); %loads structures boot and error
    FitType    = 2;%1=coarse fit, 2=fine fit
    paramOut   = param;
else
    paramOut   = NaN;
    boot    = NaN;
    error   = NaN;
    FitType = NaN;
end
