function Calc_Prob2D( RateMatrix, Dimensions, paramSetNum, trialFolder )
% This function takes a RateMatrix as input, computes the eigenvalues and associated timescales.

%   paramSetNum: counter to keep track of the simulation parameters
%   trialFolder: folder path to save results from simulations

% Formatting parameter set number
paramSetNumFormatted = sprintf('set_%05d',paramSetNum);

% Generating paths to save files
ProbVecFile=strcat(trialFolder, '/ProbVec/', paramSetNumFormatted, '.mat');
Prob2DFile=strcat(trialFolder, '/Prob2D/', paramSetNumFormatted, '.mat');
EigenValuesFile=strcat(trialFolder, '/EigenValues/', paramSetNumFormatted, '.mat');
TimeScalesFile=strcat(trialFolder, '/TimeScales/', paramSetNumFormatted, '.mat');

%Calculating NSig eigenvalues and saving eigenvalues and timescales
NSig=15;
sigma=1E-12;
dim=size(RateMatrix);
[RightEigenVectors, EigenValues] = eigs(RateMatrix,NSig,sigma,'StartVector',ones(dim(1),1));

% Calculating and saving probvec
ProbVec=RightEigenVectors(:,1)/sum(RightEigenVectors(:,1));

% Saving entropy to text file for checking against python version
Svalue = real(-sum(ProbVec.*log(ProbVec)));
EntropyFile = strcat(trialFolder,'/','sValues.txt');
dlmwrite(EntropyFile,Svalue,'-append');
% #####

% Calculating TimeScales
EigenVals_r=real(diag(EigenValues));
TimeScales=-1./EigenVals_r(2:end);

NDims=numel(Dimensions);
NStates=prod(Dimensions); %total number of states
Prob2D=zeros(Dimensions);
Prob2D(1:NStates)=ProbVec;

ListDims=1:NDims;
ProjectDims=[1,2];
OtherDims=setdiff(ListDims,ProjectDims);
RemainDims=ListDims;
for ii=1:numel(OtherDims);
    SumInd=find(ListDims==OtherDims(ii));
    Prob2D=sum(Prob2D,SumInd);
    RemainDims=setdiff(RemainDims,OtherDims(ii));
end;

save(TimeScalesFile,'TimeScales');
save(ProbVecFile,'ProbVec');
save(Prob2DFile, 'Prob2D');
save(EigenValuesFile,'EigenValues');
end

