%%%%%%%%%%%%%%
%
% Code to visualize the BOLD activity in the 2 decoupled communities 
% detected by the LEiDA patterns.
% The BOLD time series of all scans are plotted in 3D (using complex signal representation)
% Signals from brain areas in the positive community are plotted in red 
% and the others in blue.
% Green patches are added to the plot to highlight the moments when that
% LEiDA state is dominant. Lighter patches are added when the BOLD phase
% pattern is close to the centroid, even if not dominant.
%
% Statistics are performed between Baseline and Follow-up
% And Relapsed vs Remitted
%
%
% Code by Joana Cabral
% October 2020
% joanacabral@med.uminho.pt
%
%%%%%%%%%%%%%%%%%%%%%

% Set here the network of interest
k=find(rangeK==18); 
c=10;


load LEiDA_results.mat Kmeans_results Time_all
IC_time_course=Kmeans_results{k}.IDX;
Centroid_dist=Kmeans_results{k}.D;
VCentroids=Kmeans_results{k}.C;

clear Kmeans_results


VCentroid=round(VCentroids(c,:),3);

Ang_diff = [];
for s=1:n_Scans  
        
        BOLD = tc_80{s,1};
        Complex_BOLD=zeros(size(BOLD));
        
        for seed=1:N_areas
            ts = detrend(BOLD(seed,:));
            BOLD(seed,:) =filtfilt(bfilt,afilt,ts);
            Complex_BOLD(seed,:) = hilbert(BOLD(seed,:));
        end
        
        Phase_BOLD=angle(Complex_BOLD(:,2:end-1));
        Amplitude_BOLD=abs(Complex_BOLD(:,2:end-1));
        
        Ang_diff(s)=mean(abs(angdiff(mean(Phase_BOLD(VCentroid>=0,:)),mean(Phase_BOLD(VCentroid<0,:)))))*360/(2*pi);

        figure('Name',['Scan ' num2str(s)])
        plot3(0:TR:(length(BOLD)-1-2)*TR,Amplitude_BOLD(VCentroid<0,:).*sin(Phase_BOLD(VCentroid<0,:)),Amplitude_BOLD(VCentroid<0,:).*cos(Phase_BOLD(VCentroid<0,:)),'b','LineWidth',1);
        hold on
        plot3(0:TR:(length(BOLD)-1-2)*TR,Amplitude_BOLD(VCentroid>0,:).*sin(Phase_BOLD(VCentroid>0,:)),Amplitude_BOLD(VCentroid>0,:).*cos(Phase_BOLD(VCentroid>0,:)),'r','LineWidth',1);
        xlabel('Time (seconds)')
        ylabel('Imag')
        zlabel('Complex BOLD')
        xlim([0 (length(BOLD)-1-2)*TR])
        grid on
        view(0,0)
    
        % PATCH WHEN DISTANCE TO CENTROID IS SMALL
        state_ON=find(Centroid_dist(Time_all==s,c)<0.24);
        lim=max(abs(BOLD(:)));
        for t_low=1:length(state_ON)
            t_ON=state_ON(t_low);
            x1=(t_ON-1)*TR;
            x2=(t_ON)*TR;
            x = [x1 x2 x2 x1];
            y = [0 0 0 0];
            z = [-lim -lim lim lim];
            p=patch(x,y,z,'g');
            set(p,'LineStyle','none','FaceColor',[0 1 0],'FaceAlpha',0.1);
            p=patch(x,z,y,'g');
            set(p,'LineStyle','none','FaceColor',[0 1 0],'FaceAlpha',0.1);
        end
    
        % DARKER PATCH WHEN STATE IS DOMINANT
        state_ON=find(IC_time_course(Time_all==s)==c);
        lim=max(abs(BOLD(:)));
        for t_ON=state_ON
            x1=(t_ON-1)*TR;
            x2=(t_ON)*TR;
            x = [x1 x2 x2 x1];
            y = [0 0 0 0];
            z = [-lim -lim lim lim];
            p=patch(x,y,z,'g');
            set(p,'LineStyle','none','FaceColor',[0 1 0],'FaceAlpha',0.5);
            p=patch(x,z,y,'g');
            set(p,'LineStyle','none','FaceColor',[0 1 0],'FaceAlpha',0.5);
        end
    
        zlim([-lim lim])
        ylim([-lim lim])
            
end

% Compute perm-based paired-ttest statistics
a=Ang_diff(Index_t0);
b=Ang_diff(Index_t8);  
stats=permutation_htest_np_paired([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],5000,0.05,'ttest');
p_angmean=min(stats.pvals);

% Plot results
figure
for sub=1:length(a)
    plot([1 2],[a(sub) b(sub)],'-','Color',[.8 .8 .8],'MarkerSize',6);hold on
end  
box off
errorbar([1 2], [ mean(a) mean(b)],[ std(a)/sqrt(numel(a)) std(b)/sqrt(numel(b))],'LineStyle','-','Color','r','LineWidth',1)
set(gca,'XTick',[1,2],'XTickLabel',{'Baseline','Followup'},'XLim',[0.9 2.1] );
ylabel('BOLD phase shift (in degrees)')
title('Difference between Baseline and Follow-up')

