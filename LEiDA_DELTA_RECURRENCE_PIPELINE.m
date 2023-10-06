%% LEiDA_DELTA_RECURRENCE_PIPELINE
%
% 
% 1 - Read the BOLD data from the folders and computes the BOLD phases
%   - Calculate the instantaneous BOLD synchronization matrix
%   - Compute the Leading Eigenvector at each frame from all fMRI scans
% 2 - Cluster the Leading Eigenvectors into recurrent Functional Networks
% 3 - Analyse the Clustering results
%     For every fMRI scan calculate Fractional occupancy (P) and lifetimes
%     (LT) of each state c.
% 4 - Test statistical significance of within-subject changes in P and LT
% 5 - Visualize the BOLD activity in the 2 decoupled communities of a given state.
%
% Code from Joana Cabral March 2018
% joana.cabral@psych.ox.ac.uk
% adapted by sonsoles.alonso@cfin.au.dk June 2023

%% Set parameters and load data
% we start by characterizing the dynamics of the recurrent group
groupID = 5; % 5: relapse; 4: no-relapse

% Load time series and subject information
load ts_data.mat group group_labels tcT0 tcT8

Index_Patients = find(group==groupID);
n_Patients = length(Index_Patients);

tc=[tcT0(Index_Patients),tcT8(Index_Patients)]'; % concatenate all timeseries
Index_t0 = 1:n_Patients;
Index_t8 = n_Patients+1:n_Patients*2;
n_Scans = length(tc);% n_Scans = n_Patients * 2 scans
for i=1:n_Scans
    tc_80{i,1} = tc{i,1}.tc; % stores ts as matrix (regions*volumes) in a cell (n_Scans*1)
end
[N_areas,Tmax] = size(tc_80{1,1});

clear tcT0 tcT8 i tc

%% 1 - Compute the Leading Eigenvectors from the BOLD datasets
disp('Processing the eigenvectors from BOLD data') 

% Preallocate variables to save FC patterns and associated information
V1_all   = zeros((Tmax-2)*n_Scans,N_areas); % All leading eigenvectors
t_all=0; % Index of time (starts at 0 and will be updated until n_Scans*(Tmax-2))
Time_all= zeros((Tmax-2)*n_Scans,1); % Vector that links each frame to a subject

TR = 2;
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.01;                   % lowpass frequency of filter (Hz)
fhi = 0.1;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
clear fnq Wn k

for s=1:n_Scans

    % USER: Adapt to load here the BOLD matrix (NxT) from each scan
    BOLD = tc_80{s,1};

    % Get the BOLD phase using the Hilbert transform
    Phase_BOLD=zeros(N_areas,Tmax);
    for seed=1:N_areas
        ts = demean(detrend(BOLD(seed,:)));
        signal_filt =filtfilt(bfilt,afilt,ts);
        Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
    end


    % Slide over time discarding the first and last epochs
    for t=2:Tmax-1

        % Calculate the Instantaneous FC (BOLD Phase Synchrony)
        iFC=zeros(N_areas);
        for n=1:N_areas
            for p=1:N_areas
                iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
            end
        end

        % Get the leading eigenvector
        [V1,~]=eigs(iFC,1);
        if sum(V1)>0
            V1=-V1;
        end

        % Save V1 from all frames in all fMRI sessions
        t_all=t_all+1; % Update time
        V1_all(t_all,:)=V1;
        Time_all(t_all)=s;
    end
end


%% 2 - Cluster the Leading Eigenvectors

disp('Clustering the eigenvectors into')
% V1_all is a matrix containing all the eigenvectors:
% Collumns: N_areas are brain areas (variables)
% Rows: (Tmax-2)*n_Scans are all time points (independent observations)

% USER: Set maximum/minimum number of clusters
% There is no fixed number of FC states that the brain can display
% Keep the range small for the first trials
% Extend the range depending on the hypothesis of each work

% Set parameters
maxk=20;
mink=2;
rangeK=mink:maxk;
dist = 'cosine';

if groupID == 5
    % Main analysis on the patients with a relapse at 2nd scan
    Kmeans_results=cell(size(rangeK));
    for k=1:length(rangeK)
        disp(['- ' num2str(rangeK(k)) ' PL states'])
        [IDX, C, SUMD, D]=kmeans(V1_all,rangeK(k),'Distance',dist,'Replicates',100,'MaxIter',1000,'Display','off');   
        [~, ind_sort]=sort(hist(IDX,1:rangeK(k)),'descend');
        [~,idx_sort]=sort(ind_sort,'ascend');
        Kmeans_results{k}.IDX=idx_sort(IDX);   % Cluster time course - numeric collumn vectors
        Kmeans_results{k}.C=C(ind_sort,:);     % Cluster centroids (FC patterns) V
        Kmeans_results{k}.SUMD=SUMD(ind_sort); % Within-cluster sums of point-to-centroid distances
        Kmeans_results{k}.D=D(:,ind_sort);     % Distance from each point to every centroid 
    end

elseif groupID == 4
    % Validation step using data from the patients whithout a relapse at 2nd scan
    load LEiDA_results.mat Kmeans_results
    centroids = Kmeans_results;
    Kmeans_results=cell(size(rangeK));
    for k=1:length(rangeK)
        disp(['- ' num2str(rangeK(k)) ' PL states'])
        [IDX, C, SUMD, D]=kmeans(V1_all,rangeK(k),'Distance',dist,'MaxIter',1,'Display','off','Start',centroids{k}.C);
        Kmeans_results{k}.IDX=IDX;   % Cluster time course - numeric collumn vectors
        Kmeans_results{k}.C=centroids{k}.C;     % Cluster centroids (FC patterns) V
        Kmeans_results{k}.SUMD=SUMD; % Within-cluster sums of point-to-centroid distances
        Kmeans_results{k}.D=D;     % Distance from each point to every centroid 
    end
end

%% 3. - Analyse the Clustering results
% For every fMRI scan calculate Fractional occupancy (P) and lifetimes of each state c.
P=zeros(n_Scans,maxk-mink+1,maxk);
LT=zeros(n_Scans,maxk-mink+1,maxk);

for k=1:length(rangeK)
    for s=1:n_Scans

        % Select the time points representing this subject and task
        T=(Time_all==s);
        Ctime=Kmeans_results{k}.IDX(T);

        for c=1:rangeK(k)
            % Fractional Occupancy (P)
            P(s,k,c)=mean(Ctime==c);

            % Mean Lifetime
            Ctime_bin=Ctime==c;

            % Detect switches in and out of this state
            a=find(diff(Ctime_bin)==1);
            b=find(diff(Ctime_bin)==-1);

            % We discard the cases where state sarts or ends ON
            if length(b)>length(a)
                b(1)=[];
            elseif length(a)>length(b)
                a(end)=[];
            elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                b(1)=[];
                a(end)=[];
            end
            if ~isempty(a) && ~isempty(b)
                C_Durations=b-a;
            else
                C_Durations=0;
            end
            LT(s,k,c)=mean(C_Durations)*TR; %time in seconds
        end
    end
end

%% 4. - Tests statistical significance

% Compare Fractional Occupancy (P)
P_pval=zeros(length(rangeK),max(rangeK));
P_pval_bh=zeros(length(rangeK),max(rangeK));
P_cohenD=zeros(length(rangeK),max(rangeK));
disp('Testing between-scan differences in State Fractional Occupancy')
for k=1:length(rangeK)
    disp(['Now running statistics for ' num2str(rangeK(k)) ' PL states'])
    for c=1:rangeK(k)
        a=squeeze(P(Index_t0,k,c))';  % Vector containing Prob of c at baseline
        b=squeeze(P(Index_t8,k,c))';  % Vector containing Prob of c st follow-up
        stats=permutation_htest_np_paired([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],5000,0.05,'ttest');
        P_pval(k,c)=min(stats.pvals);
        P_cohenD(k,c)= computeCohen_d(a, b, 'paired'); 
    end
    [~,~,~, P_pval_bh(k,1:rangeK(k))]=fdr_bh(P_pval(k,1:rangeK(k)),0.05,'pdep','yes');   
end

% Compare Lifetime or Duration
LT_pval=zeros(length(rangeK),max(rangeK));
LT_pval_bh=zeros(length(rangeK),max(rangeK));
LT_cohenD=zeros(length(rangeK),max(rangeK));
disp('Testing between-scan differences in State Lifetime')
for k=1:length(rangeK)
    disp(['Now running statistics for ' num2str(rangeK(k)) ' PL states'])
    for c=1:rangeK(k)
        a=squeeze(LT(Index_t0,k,c))';  % Vector containing Prob of c at baseline
        b=squeeze(LT(Index_t8,k,c))';  % Vector containing Prob of c at follow-up
        stats=permutation_htest_np_paired([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],5000,0.05,'ttest');
        LT_pval(k,c)=min(stats.pvals);
        LT_cohenD(k,c)= computeCohen_d(a, b, 'paired'); 
    end
    [~,~,~, LT_pval_bh(k,1:rangeK(k))]=fdr_bh(LT_pval(k,1:rangeK(k)),0.05,'pdep','yes');
end

% save results
if groupID == 5
    save LEiDA_results.mat N_areas TR Tmax n_Scans  n_Patients Index_t0 Index_t8 subjectsID ...
     tc_80 dist rangeK mink maxk V1_all Time_all Kmeans_results ...
     P   P_pval  P_pval_bh P_cohenD LT LT_pval LT_pval_bh LT_cohenD
elseif groupID == 4
    save LEiDA_results_validation.mat N_areas TR Tmax n_Scans  n_Patients Index_t0 Index_t8 subjectsID ...
     tc_80 dist rangeK mink maxk V1_all Time_all Kmeans_results ...
     P   P_pval  P_pval_bh P_cohenD LT LT_pval LT_pval_bh LT_cohenD
end

%% 5. - Visualize BOLD activity of the two decoupled communities of a given state
% see get_phaseShift.m

