clear all
% clear variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% This is a code demo for exactly simulating qudits measurement procedures
%
% For more information, see README.md
%
% % Q Huang, J Sun, and xxx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Global parameters and Paths included
global Nq dim psi0
global observ num_observ ObservMat
global pauli_vec % alias for "observ", due to past code convention
global GGM Veigens Deigens
global prob_expm_accum

% for LBCS
global beta_dis

addpath(genpath('./utils'));
addpath(genpath('./measPrescribe'));
addpath(genpath('./LBCS'));
addpath(genpath('./Unbiased'));
addpath(genpath('./DerandPy'));

%% Parameter settings and Initialisation

% ------ To Set for Every Run ------
Nq = 4; % number of qudits
dim = 2; % qudit dimension. 2 for qubit, 3 for qutrit, ...

% Initialise observables
% fid = 'paulivec_Nq4_dim2.txt';
% fid = ['paulivec_Nq',num2str(Nq),'_dim', num2str(dim),'.txt'];
fid = 'simpleHamil.txt'; % txt file for Hamiltonian, or other set of observables

% Data generation mode
mode_data = 1;
% 0: Experiment mode (TODO), read & load quantum state info from experimental data
% 1: Theoretical mode, derive from ideal quantum pure state (in our case, GHZ)

% -----------------------------------

[GGM, Veigens, Deigens] = get_GellMann_Matrix(dim); % Initialise Generalised Gell-Mann Matrices

psi0 = get_GHZ(); % Initialise psi0 as a pure GHZ state, i.e. a column vector

observ = load(fid); % observables, following {{observ_matlab}} format, similar to
% e.g. {{observ.txt}}
% 
% 0  3  7  0  0.2452102298125
% 0  0  4  8  0.1412847218385
% ...
%
num_observ = size(observ,1);
pauli_vec = observ; % alias

[observ_py, weight_py] = hamil_manip(observ); % python related manipulation

if mode_data > 0.5
    % quantum state considered is an ideally pure GHZ state
    ObservMat = zeros(dim^Nq, dim^Nq, num_observ);
    for i = 1:num_observ
        tmp = 1;
        for k = 1:Nq
            tmp = sparse(kron(tmp, sparse(GGM(:,:,observ(i,k)+1))));
        end
        ObservMat(:,:,i) = tmp; 
    end
    
    % Generate exact measurement outcomes
    Ideal = zeros(num_observ,1); % exact outcomes for every observable
    En_ideal = 0; % exact overall outcome
    for j = 1:num_observ
        Ideal(j) = psi0' * sparse(ObservMat(:,:,j)) * psi0;
        En_ideal = En_ideal + observ(j, end) * Ideal(j);
    end

    % Probabilistic measurement prescription
    prob_expm_accum = zeros((dim^2-1)^Nq, dim^Nq);
    fid_measPrescribe = ['./measPrescribe/measPrescribeGHZ_Nq', num2str(Nq), '_dim', num2str(dim), '.mat'];
    load(fid_measPrescribe); % stored in prob_expm_accum, automatically normalised

    % generated by:
%     for i = 1:(dim^2-1)^Nq
%         for j = 1:dim^Nq % i -> Basis, j -> Outcome
%             prob_expm_accum(i,j) = abs(psi0' * get_product_state(i-1,j-1))^2;
%         end
%     end

elseif mode_data < 0.5

    % To write later. - Q Huang

end

%% Error vs number of Samples

% ------ To Set for Every Run ------
Sample_list = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000];
% Sample_list = [200, 400, 600, 800, 1000];
% Sample_list = [2000,];
% floor((2:2:20) * 10^2); % List of samples

L = 20; % independent repetitions

% LBCS optimisation
beta_dis = get_beta(1000, 0.01, dim); % Minimises Variance
% outputs a (Nq * (dim^2-1)) matrix of beta's

% -----------------------------------

% TODO: add variance analyses to shadow methods, while not necessary.


% Currently available:
% Unbiased Shadow, Locally-biased Shadow, DerandQudit, Locally-biased DerandQudit
% strictly in this order

MaxErrStat = zeros(L, length(Sample_list),4);
AllErrStat = zeros(L, length(Sample_list),4);

MaxErrMean = zeros(length(Sample_list),4);
AllErrMean = zeros(length(Sample_list),4);

MaxErrVar = zeros(length(Sample_list),4);
AllErrVar = zeros(length(Sample_list),4);


for rep = 1:L
    fprintf('Now in the %d-th round, out of a total of %d.\n', rep, L);

    Error_Unbiased = [];
    Max_Error_Unbiased = [];
    
    Error_LBCS = [];
    Max_Error_LBCS = [];
    
    Derand_Sample_num = [];
    Error_Derand = [];
    Max_Error_Derand = [];
    
    DerandLB_Sample_num = [];
    Error_DerandLB = [];
    Max_Error_DerandLB = [];
    
    fprintf('\n\n Sample \t Unbiased \t LBCS \t Derand \t LB Derand\n');
    fprintf('----------------------------------------------------------\n\n');

    % Exact Simulation of Measurement Procedures
    for Sample = Sample_list
    
        % Exact Unbiased Shadow procedure
        basis_list = floor((dim^2-1)*rand(Sample,Nq)+1); % Exclude Identity
        fweight = get_fweight(basis_list,pauli_vec, dim);
    
        mu = Eig_Measure(basis_list, dim);

        [En,Ol] = postShadow_Hamil(mu, Sample, 1:num_observ, fweight);
    
        Max_Error_Unbiased = [Max_Error_Unbiased, max(abs(real(Ol') - real(Ideal)))];
        Error_Unbiased = [Error_Unbiased, En - En_ideal];
    
    
        % Exact LBCS procedure
        basis_list = get_basis_list(Sample, dim); % Get measurement basis
        fweight = get_fweight_random(basis_list, pauli_vec); % Compute f-weights for LBCShadow
    
        mu = Eig_Measure(basis_list, dim); % get measurement outcomes

        [En,Ol] = postShadow_Hamil(mu, Sample, 1:num_observ, fweight);


        Max_Error_LBCS = [Max_Error_LBCS, max(abs(real(Ol') - real(Ideal)))]; % Max Error within Observables
        Error_LBCS = [Error_LBCS, En - En_ideal]; % LBCS Error for Energy
    
    
        % Qudit Derandomisation, with a definite measurement budget
        measurement_budget_py = py.int(Sample);
        % num_of_measurements_per_observable_py = py.int(floor(Sample / num_observ));
        system_size_py = py.int(Nq);
        dim_py = py.int(dim);
    
        [all_observables_py, weight_py] = hamil_manip(observ);
    
        basis_list_py = py.derandTools.derandomized_classical_shadow_qudit(all_observables_py, ...
            measurement_budget_py, system_size_py, pyargs('dim', dim_py));
        
        basis_list = basis_manip(basis_list_py);

        mu = Eig_Measure(basis_list, dim);
    
        [En, Ol] = postDerand_Hamil(mu, size(basis_list,1), 1:num_observ, basis_list);
    
            % Derand_Sample_num = [Derand_Sample_num, size(basis_list,1)];
    
        Max_Error_Derand = [Max_Error_Derand, max(abs(real(Ol') - real(Ideal)))]; % Max Error within Obersvables
        Error_Derand = [Error_Derand, En - En_ideal]; % Energy
    
        
        % Qudit Derandomisation with Local Bias, with definite measurement budget
        local_bias_py = bias_manip(beta_dis);
    
        basis_list_py = py.derandTools.derandomized_classical_shadow_qudit(all_observables_py, ...
            measurement_budget_py, system_size_py, pyargs('dim', dim_py, 'local_bias', local_bias_py));
    
        basis_list = basis_manip(basis_list_py);

        mu = Eig_Measure(basis_list, dim);

        [En, Ol] = postDerand_Hamil(mu, size(basis_list,1), 1:num_observ, basis_list);

            % DerandLB_Sample_num = [DerandLB_Sample_num, size(basis_list,1)];
    
        Max_Error_DerandLB = [Max_Error_DerandLB, max(abs(real(Ol') - real(Ideal)))]; % Max Error within Obersvables
        Error_DerandLB = [Error_DerandLB, En - En_ideal]; % Energy

        fprintf('%d \t %f \t %f \t %f \t %f \n', Sample, Error_Unbiased(end), Error_LBCS(end), Error_Derand(end), Error_DerandLB(end))
    end
    MaxErrStat(rep,:,1) = Max_Error_Unbiased';
    MaxErrStat(rep,:,2) = Max_Error_LBCS';
    MaxErrStat(rep,:,3) = Max_Error_Derand';
    MaxErrStat(rep,:,4) = Max_Error_DerandLB';

    AllErrStat(rep,:,1) = abs(Error_Unbiased');
    AllErrStat(rep,:,2) = abs(Error_LBCS');
    AllErrStat(rep,:,3) = abs(Error_Derand');
    AllErrStat(rep,:,4) = abs(Error_DerandLB');
    
end

MaxErrMean = reshape(mean(MaxErrStat,1), [length(Sample_list), 4]);
AllErrMean = reshape(mean(AllErrStat,1), [length(Sample_list), 4]);

for method = 1:size(MaxErrStat, 3)
    for sample = 1:size(MaxErrStat, 2)
        MaxErrVar(sample,method) = get_stat_variance(MaxErrStat(:,sample,method));
        AllErrVar(sample,method) = get_stat_variance(AllErrStat(:,sample,method));
    end
end

close all
