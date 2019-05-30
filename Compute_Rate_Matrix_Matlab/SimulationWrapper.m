function [] = SimulationWrapper(trialFolder, modelFile)
% This wrapper should call all necessary scripts to calculate the RateMatrix for a given model
%function [] = SimulationWrapper()
%trialFolder = "../outputs_tmp/Trial_0001";
%modelFile = "Compute_RateMatrix_MISAEx";
% Saving paths to Map object corresponding to bash variable names, see directoryPaths.keys and .values. Adding model folder path.
addpath('models/', '-end');

% Created dynamic function handle
rateMatrixCalc = str2func(modelFile);

% Loading parameters file
parametersFile = strcat(trialFolder,'/','paramValues.csv');
parameters = readtable(parametersFile);
paramNames = parameters.Properties.VariableNames;
tic
for row = 1:size(parameters,1)
    parameter_assignments = containers.Map(paramNames, parameters(row,:).Variables);
	paramSetNum = int16(row); 

    %Value check in prompt
	%disp(['Values: ', num2str(paramSetNum), ', ' ,num2str(parameter_assignments('ha')), ',', num2str(parameter_assignments('hr')), ',', num2str(parameter_assignments('fa')),',', num2str(parameter_assignments('fr'))]) 
    %disp('Values hardcoded')
	
	% Step 1: Specify the model and compute the Rate Matrix:
    % Everything about the model is specified in this script, outputs RateMatrix.mat and StatesList.mat
    [ModelName, RateMatrix, Dimensions] = rateMatrixCalc(paramSetNum, parameter_assignments, trialFolder);

	% Step 2: Calculate ProbVec and Prob2D
	Calc_Prob2D (RateMatrix, Dimensions, paramSetNum, trialFolder);

end
toc
end
