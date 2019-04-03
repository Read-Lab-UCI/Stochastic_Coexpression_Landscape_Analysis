function [ ModelName, RateMatrix, Dimensions ] = Compute_RateMatrix_MISAInc( paramSetNum, parameters, trialFolder)
% We are interested in the stationary solution to the Master Equation, in
% vector-matrix form, dP/dt=AP, where P(t) is a time-dependent vector
% containing the probability for the system to be in each state. We can
% solve the stationary "probability landscape" by (1) enumerating the
% reaction rate matrix A and computing the stochastic reaction propensities
% and (2) solving the eigenproblem, A*phi=lambda*phi. P(t->inf.)=0 is the
% eigenvector phi corresponding to lambda=0

% ha is the binding rate of activators
% hr is the binding rate of repressors
% fa is the off-rate of activators
% fr is the off-rate of repressors
% N is the max copy-number of transcription factors (TFs). It should be set to something that safely exceeds g1/k
% g0 is the basal rate of production of TFs
% g1 is the acvtivated rate of production of TFs
% kd is the degradation/dilution rate of TFs (keep this at 1, because the other parameters can be scaled accordingly

% paramSetNum and trialFolder are for keeping track during multiple simulations.

Model='MISAInc'; % Used for filenames 
N=parameters('N');
g0=parameters('g0');
g1=parameters('g1');
kd=parameters('kd');
ha=parameters('ha');
hr=parameters('hr');
fa=parameters('fa');
fr=parameters('fr');

ModelName=strcat(Model,'_N',num2str(N));

% Below: Enumerating the state-space for the model. 
% There are two genes; A and B. 
% The state of the system is specified by:
% [n_A,n_B,GeneState_A,GeneStateB] 
% n_A and n_B are copy numbers of A and B, currently both equal to N
% GeneState_A and GeneStateB are the total number of possible gene states 

low = 0;
A=low:N;
B=low:N;

GeneA_00=[0,1]; GeneA_01=[0,1]; GeneA_10=[0,1]; GeneA_11=[0,1];
GeneB_00=[0,1]; GeneB_01=[0,1]; GeneB_10=[0,1]; GeneB_11=[0,1];

GeneA = 0:1;
GeneB = 0:1;
BoundGeneA = [0,1];
BoundGeneB = [0,1];

% Total number of states:
NS=numel(A)*numel(B)*4*4;

% Range of accessible states:
Smalls=[A(1),B(1),GeneA_00(1),GeneA_01(1),GeneA_10(1),GeneA_11(1),GeneB_00(1),GeneB_01(1),GeneB_10(1),GeneB_11(1)];
Bigs=[A(end),B(end),GeneA_00(end),GeneA_01(end),GeneA_10(end),GeneA_11(end),GeneB_00(end),GeneB_01(end),GeneB_10(end),GeneB_11(end)];

% Specifying the Reaction Network--16 reactions
NumRxn=26;  % Number of reactions
NumSpec=10;  % Number of species: A,B,GeneA_00,GeneA_01,GeneA_10,GeneB_00,GeneB_01,GeneB_10

% Reactions and Stoich represented in Rxn.Law and Rxn.Stoich are listed in MISA_Reactions_NonCompetitive.txt
Parameters=[g0;g1;g0;g0;ha;hr;hr;ha;fa;fr;fa;fr;kd;g0;g1;g0;g0;ha;hr;hr;ha;fa;fr;fa;fr;kd];
Rxn.Par=Parameters; % Rxn is a struct datatype that holds the model information

% Reaction rate laws, number of each species involved in the reaction
Rxn.Law = zeros(NumRxn,NumSpec);
Rxn.Law(1,3)=1;
Rxn.Law(2,4)=1;
Rxn.Law(3,5)=1;
Rxn.Law(4,6)=1;
Rxn.Law(5,1)=2; Rxn.Law(5,3)=1;
Rxn.Law(6,2)=2; Rxn.Law(6,3)=1;
Rxn.Law(7,2)=2; Rxn.Law(7,4)=1;
Rxn.Law(8,1)=2; Rxn.Law(8,5)=1;
Rxn.Law(9,4)=1;
Rxn.Law(10,5)=1;
Rxn.Law(11,6)=1;
Rxn.Law(12,6)=1;
Rxn.Law(13,1)=1;

Rxn.Law(14,7)=1;
Rxn.Law(15,8)=1;
Rxn.Law(16,9)=1;
Rxn.Law(17,10)=1;
Rxn.Law(18,2)=2; Rxn.Law(18,7)=1;
Rxn.Law(19,1)=2; Rxn.Law(19,7)=1;
Rxn.Law(20,1)=2; Rxn.Law(20,8)=1;
Rxn.Law(21,2)=2; Rxn.Law(21,9)=1;
Rxn.Law(22,8)=1;
Rxn.Law(23,9)=1;
Rxn.Law(24,10)=1;
Rxn.Law(25,10)=1;
Rxn.Law(26,2)=1;


% Reaction stoichiometry, change in species resulting from reaction
Rxn.Stoich = zeros(NumRxn,NumSpec);
% Gene A Production (1-4), Binding (5-8), Unbinding (9-12), and Degradation (13)
Rxn.Stoich(1,1) = 1; 
Rxn.Stoich(2,1) = 1; 
Rxn.Stoich(3,1) = 1;
Rxn.Stoich(4,1) = 1;
Rxn.Stoich(5,3) = -1; Rxn.Stoich(5,1) = -2; Rxn.Stoich(5,4) = 1;
Rxn.Stoich(6,3) = -1; Rxn.Stoich(6,2) = -2; Rxn.Stoich(6,5) = 1;
Rxn.Stoich(7,4) = -1; Rxn.Stoich(7,2) = -2; Rxn.Stoich(7,6) = 1;
Rxn.Stoich(8,5) = -1; Rxn.Stoich(8,1) = -2; Rxn.Stoich(8,6) = 1;
Rxn.Stoich(9,4) = -1; Rxn.Stoich(9,1) = 2; Rxn.Stoich(9,3) = 1;
Rxn.Stoich(10,5) = -1; Rxn.Stoich(10,2) = 2; Rxn.Stoich(10,3) = 1;
Rxn.Stoich(11,6) = -1; Rxn.Stoich(11,1) = 2; Rxn.Stoich(11,5) = 1;
Rxn.Stoich(12,6) = -1; Rxn.Stoich(12,2) = 2; Rxn.Stoich(12,4) = 1;
Rxn.Stoich(13,1) = -1;

% Gene B Production (14-17), Binding (18-21), Unbinding (22-25), and Degradation (26)
Rxn.Stoich(14,2) = 1;
Rxn.Stoich(15,2) = 1;
Rxn.Stoich(16,2) = 1;
Rxn.Stoich(17,2) = 1;
Rxn.Stoich(18,7) = -1; Rxn.Stoich(18,2) = -2; Rxn.Stoich(18,8) = 1;
Rxn.Stoich(19,7) = -1; Rxn.Stoich(19,1) = -2; Rxn.Stoich(19,9) = 1;
Rxn.Stoich(20,8) = -1; Rxn.Stoich(20,1) = -2; Rxn.Stoich(20,10) = 1;
Rxn.Stoich(21,9) = -1; Rxn.Stoich(21,2) = -2; Rxn.Stoich(21,10) = 1;
Rxn.Stoich(22,8) = -1; Rxn.Stoich(22,2) = 2; Rxn.Stoich(22,7) = 1;
Rxn.Stoich(23,9) = -1; Rxn.Stoich(23,1) = 2; Rxn.Stoich(23,7) = 1;
Rxn.Stoich(24,10) = -1; Rxn.Stoich(24,2) = 2; Rxn.Stoich(24,9) = 1;
Rxn.Stoich(25,10) = -1; Rxn.Stoich(25,1) = 2; Rxn.Stoich(25,8) = 1;
Rxn.Stoich(26,2) = -1;

%% Defining the reaction rate matrix, A. The elements will be computed by
%% looping over all states in state-space and checking allowed reactions
%% that can lead to that state from all other states.

% Initialize reaction rate matrix: 
RateMatrix=sparse(NS,NS);

GeneA_States=[1,0,0,0;
              0,1,0,0;
              0,0,1,0;
              0,0,0,1];
GeneB_States=[1,0,0,0;
              0,1,0,0;
              0,0,1,0;
              0,0,0,1];

StatesList=zeros(NS,NumSpec);
Dimensions=[numel(A),numel(B),length(GeneA_States),length(GeneB_States)];
n=0;
for ii=1:numel(A);
    for jj=1:numel(B);
        for kk=1:length(GeneA_States);
            for ll=1:length(GeneB_States);
                n=n+1;
                % Find the current state of the system:
                Cur=[A(ii),B(jj),GeneA_States(kk,:),GeneB_States(ll,:)];
                % Find the current state index:
                CurInd=sub2ind([length(A),length(B),length(GeneA_States),length(GeneB_States)],ii,jj,kk,ll);
                % List all states:
                StatesList(CurInd,:)=Cur;
                % Loop over all possible reactions:
                for mm=1:NumRxn;
                    % Check if the destination state is in allowed range,
                    % using reaction stoichiometries
                    TestDest=Cur+Rxn.Stoich(mm,:);
                    TestinRange=[TestDest>=Smalls,TestDest<=Bigs];
                    if prod(+TestinRange)>0             
                        Dest=TestDest;
                        GeneA_Dest=TestDest(3:6);
                        GeneB_Dest=TestDest(7:10);
                        [c,GAind]=find(GeneA_Dest);
                        [c,GBind]=find(GeneB_Dest);
                        Aind=Dest(1)+1;
                        Bind=Dest(2)+1;
                        inds=[Aind,Bind,GAind,GBind];
                        DestInd=sub2ind([length(A),length(B),length(GeneA_States),length(GeneB_States)],Aind,Bind,GAind,GBind);
                        par=Rxn.Par(mm);
                        law=Rxn.Law(mm,:);
                        % Compute the reaction propensity using parameters and rate laws
                        rate0=ones(NumSpec,1);
                        rate0=(Cur.^law-max(law-1,0).*Cur.^(max(law-1,0)))./max(1,law);
                        Rate=par*prod(rate0);   
                        % Place the computed propensity in the Rate Matrix
                        RateMatrix(DestInd,CurInd)=RateMatrix(DestInd,CurInd)+Rate;
                    end;
                end;
            end;
        end;
    end;
end;
% Making columns in rate matrix sum to 0
RateMatrix=RateMatrix-diag(sum(RateMatrix));

paramSetNumFormatted = sprintf('set_%05d',paramSetNum);
paramSetNumFormatted = sprintf('set_%05d',paramSetNum);
MatrixFile=strcat(trialFolder, '/RateMatrix/', paramSetNumFormatted, '.mat');
%save('RateMatrix_ml.mat', 'RateMatrix') 
save(MatrixFile,'RateMatrix','Dimensions','parameters')

end
