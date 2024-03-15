%% POTT Pref - Selective satiety session analyses
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2022.12
%- Related to: Stoll & Rudebeck, Neuron, 2024

clear

mk = 'Mimic';

%- locate the files
path2go = ('/Users/Fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data-final/')
%path2go = ('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\')

%- list all the files
if strcmp(mk(1:2),'Mi')
    list = dir([path2go 'X*e_behav.mat']);
end

showfig = false;

%- reorganize list by date
days = [];
for ss = 1 : length(list)
    days = [days;datenum(list(ss).name(2:7),'mmddyy') , ss];
end
date_order = sortrows(days);
list = list(date_order(:,2));

for sess = 1 : length(list) % for every sessions
    filename = [path2go list(sess).name];

    [ALL_temp,param] = behav_model_bins_deval(filename,false,showfig); %- no permutations
    ALL(sess,:) = ALL_temp';

    dumm = load(filename);
    tim(sess,1).t_FPfix = dumm.t_FPfix;    
end
save([path2go mk '_behav_norm_DEVAL_prevJuice.mat'],'ALL','param','tim','-v7.3')


%% POST-PROCESSING OF THE DEVAL SESSIONS IN MONKEY X

clear
path2go = ('/Users/Fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data-final/')

cd(path2go)
load('Mimic_behav_norm_DEVAL_prevJuice.mat')

%- which juice was given (from notebook)
bolus_type = [0 0 1 1 0 0]; %- 0 = juice2 / 1 = juice1

%- extract data
for s = 1 : length(ALL)
    pre{s} = ALL(s,1).pref_bins; %- remove 'bins' to stop before bolus and avoid effect of smoothing
    post{s} = ALL(s,2).pref_bins;
    pre_IT{s} = ALL(s,1).ft_bins; %- remove 'bins' to stop before bolus and avoid effect of smoothing
    post_IT{s} = ALL(s,2).ft_bins;
    pre_RT{s} = ALL(s,1).rt_bins; %- remove 'bins' to stop before bolus and avoid effect of smoothing
    post_RT{s} = ALL(s,2).rt_bins;
end

%- inverse so that it is express in juice given preference
pre_all =NaN(6,500);
post_all=NaN(6,500);
for s = 1 : length(ALL)
    if bolus_type(s)==1
        pre_all(s,end-length(pre{s})+1:end) = 1-pre{s};
        post_all(s,1:length(post{s})) = 1-post{s};
    else
        pre_all(s,end-length(pre{s})+1:end) = pre{s};
        post_all(s,1:length(post{s})) = post{s};
    end
end

%- compute the trend statistic before and afer
for i = 1 : length(pre)
    if bolus_type(i)==1 %- reverse to express in juice given
        [Z_trend(i),p_trend(i)]=Mann_Kendall(1-pre{i},0.01);
        [Z_trend_post(i),p_trend_post(i)]=Mann_Kendall(1-post{i},0.01);
    else
        [Z_trend(i),p_trend(i)]=Mann_Kendall(pre{i},0.01);
        [Z_trend_post(i),p_trend_post(i)]=Mann_Kendall(post{i},0.01);
    end
end
% invert the Z values when negative slope during pre
% Forced to have a positive trend during PRE. 
invert = sign(Z_trend)==-1;
Z_trend(invert)=-Z_trend(invert);
Z_trend_post(invert)=-Z_trend_post(invert);

[p,h,st] = signrank(Z_trend,Z_trend_post,'method','approximate')
figure;subplot(1,3,1);plot([Z_trend;Z_trend_post],'.-','Color','k','MarkerSize',20)
xlim([.5 2.5])
set(gca,'XTick',[1 2],'XTickLabel',{'PRE' 'POST'},'FontSize',16)
ylabel('Z-value (trend in Juice Preference)')
text(1.5,double(round(max(Z_trend))+1),['p=' num2str(p)],'HorizontalAlignment','center','FontSize',16)
%- this means that there is an overall significant decrease in the trend after
%- being sated compared to before. Note however than in 2 cases, he was
%developping a negative trend before bolus (like the other more), which
%became less negative after. 

%- compute the trend statistic before and afer for IT
clear Z_trend Z_trend_post
for i = 1 : length(pre)
    [Z_trend(i),p_trend(i)]=Mann_Kendall(pre_IT{i},0.01);
    [Z_trend_post(i),p_trend_post(i)]=Mann_Kendall(post_IT{i},0.01);
end

[p,h,st] = signrank(Z_trend,Z_trend_post,'method','approximate')
subplot(1,3,2);plot([Z_trend;Z_trend_post],'.-','Color','k','MarkerSize',20)
xlim([.5 2.5])
set(gca,'XTick',[1 2],'XTickLabel',{'PRE' 'POST'},'FontSize',16)
ylabel('Z-value (trend in IT)')
text(1.5,double(round(max(Z_trend))+1),['p=' num2str(p)],'HorizontalAlignment','center','FontSize',16)

%- compute the trend statistic before and afer for RT
clear Z_trend Z_trend_post
for i = 1 : length(pre)
    [Z_trend(i),p_trend(i)]=Mann_Kendall(pre_RT{i},0.01);
    [Z_trend_post(i),p_trend_post(i)]=Mann_Kendall(post_RT{i},0.01);
end

[p,h,st] = signrank(Z_trend,Z_trend_post,'method','approximate')
subplot(1,3,3);plot([Z_trend;Z_trend_post],'.-','Color','k','MarkerSize',20)
xlim([.5 2.5])
set(gca,'XTick',[1 2],'XTickLabel',{'PRE' 'POST'},'FontSize',16)
ylabel('Z-value (trend in RT)')
text(1.5,double(round(max(Z_trend))+1),['p=' num2str(p)],'HorizontalAlignment','center','FontSize',16)

%- plot the pref changes for all sessions
figure;
for s = 1:length(ALL)
    subplot(1,6,s)
    if bolus_type(s)
        plot(1:length(pre{s}),1-pre{s},'k');hold on;
        plot(length(pre{s})+1:(length(pre{s})*2),1-post{s},'r');
    else
        plot(1:length(pre{s}),pre{s},'k');hold on;
        plot(length(pre{s})+1:(length(pre{s})*2),post{s},'r');
    end
    ylim([0 1]);
    line([1 101],[.5 .5],'Color','k')
    ylabel('Preference (Given Juice when > 0.5)');
    xlabel('Trials');
end

%- plot the pref changes for all sessions
figure;x = 0;
for s = 1:length(ALL)
    x = x + 1;
    subplot(3,6,x)
    if bolus_type(s)
        plot(1:length(pre{s}),1-pre{s},'k');hold on;
        plot(length(pre{s})+1:(length(pre{s})*2),1-post{s},'r');
    else
        plot(1:length(pre{s}),pre{s},'k');hold on;
        plot(length(pre{s})+1:(length(pre{s})*2),post{s},'r');
    end
    ylim([0 1]);
    line([1 101],[.5 .5],'Color','k')
    ylabel('Preference (Given Juice when > 0.5)');
    xlabel('Trials');

    subplot(3,6,x+6)
    plot(1:length(pre{s}),pre_IT{s},'k');hold on;
    plot(length(pre{s})+1:(length(pre{s})*2),post_IT{s},'r');
    ylabel('IT');
    xlabel('Trials');

    subplot(3,6,x+12)
    plot(1:length(pre{s}),pre_RT{s},'k');hold on;
    plot(length(pre{s})+1:(length(pre{s})*2),post_RT{s},'r');
    ylabel('RT');
    xlabel('Trials');
end

% plot choice vs proba difference for all sessions
prob = [10:10:90];
figure;
for s = 1 : length(ALL)
    %- retrieve some variables
    %choice = abs(double(ALL(s).T.choice)-1);
    if bolus_type(s)==1 %- inverse to express in given vs non given juice
        probaG = [ALL(s,1).T.probJ2 ; ALL(s,2).T.probJ2] ;
        probaU = [ALL(s,1).T.probJ1 ; ALL(s,2).T.probJ1];
        choice = [abs(double(ALL(s,1).T.choice)-2) ; abs(double(ALL(s,2).T.choice)-2) ];
    else
        probaG = [ALL(s,1).T.probJ1 ; ALL(s,2).T.probJ1] ;
        probaU = [ALL(s,1).T.probJ2 ; ALL(s,2).T.probJ2];
        choice = [abs(double(ALL(s,1).T.choice)-1) ; abs(double(ALL(s,2).T.choice)-1) ];
    end
    %- unnormalize the proba
    prob_norm = unique(probaG);
    for ii = 1 : length(prob)
        probaG(probaG==prob_norm(ii))=prob(ii);
        probaU(probaU==prob_norm(ii))=prob(ii);
    end
    
    %- Pre/post and quartile groups
    prepost = [zeros(1,size(ALL(s,1).T,1)) ones(1,size(ALL(s,2).T,1)) ];
    quarti = [0:1/(size(ALL(s,1).T,1)-1):1 , 0:1/(size(ALL(s,2).T,1)-1):1 ] ;

    bins_pb = -80 : 10 : 80;
    clear probaCh_bins tot_bins nb_bins
    for i = 1 : length(bins_pb)
        samePb =   ((probaG-probaU)==bins_pb(i))';
        %- only consider the last half of the PRE bolus period
        chSamePb = [sum(choice(samePb & prepost==0 & quarti>=0.5)) length(choice(samePb & prepost==0 & quarti>=0.5)) ;
                    sum(choice(samePb & prepost==1 )) length(choice(samePb & prepost==1 )) ];
                
        probaCh_bins(:,i) = chSamePb(:,1)./chSamePb(:,2);
        tot_bins(:,i) = chSamePb(:,2);
        nb_bins(:,i) = chSamePb(:,1);
    end
    pre_pref(s) = nanmean(pre{s}(round(length(pre{s})/2):end)); %- only take the 2nd half of the pre

    weights = ones(numel(bins_pb),1);
    coeffs_1 = glmfit(bins_pb, [nb_bins(1,:)', tot_bins(1,:)'], 'binomial','link','logit');
    coeffs_2 = glmfit(bins_pb, [nb_bins(2,:)', tot_bins(2,:)'], 'binomial','link','logit');
    
    % Create a new xAxis with higher resolution
    fineX = min(bins_pb):1:max(bins_pb);
    
    % Generate curve from fit
    curve_1 = glmval(coeffs_1, fineX, 'logit');
    curve_1 = [fineX', curve_1];
    curve_2 = glmval(coeffs_2, fineX, 'logit');
    curve_2 = [fineX', curve_2];
    
    subplot(2,6,s)
    plot(bins_pb,probaCh_bins(1,:),'.','Color',[.4 .4 .4],'MarkerSize',30);hold on
    plot(bins_pb,probaCh_bins(2,:),'.','Color',[255 100 100]/255,'MarkerSize',30)
    plot(curve_1(:,1),curve_1(:,2),'k','LineWidth',2)
    plot(curve_2(:,1),curve_2(:,2),'r','LineWidth',2);
    xlabel('Proba difference (Giv-Ungiv)');
    ylabel('GivenJuice Choice Probability')
    xlim([-80 80])
    subplot(2,6,s+6)
    plot(curve_1(:,1),curve_2(:,2)-curve_1(:,2),'k','LineWidth',2)
    diff_curve(s,:) = (curve_2(:,2)-curve_1(:,2))';
    diff_proba(s) = nansum(probaCh_bins(1,:)-probaCh_bins(2,:));
end

%- Plot change in pref and pre preference
ch = sum(diff_curve')./length(diff_curve');
avg_ch = median(ch);
sem_ch = std(ch)/sqrt(length(diff_curve(:,1)));
avg_pr = median(pre_pref);
sem_pr = std(pre_pref)/sqrt(length(diff_curve(:,1)));

figure;plot(pre_pref,ch,'.k','MarkerSize',20)
hold on
line([0 1],[0 0],'Color','k')
plot(avg_pr,avg_ch,'.r','MarkerSize',30)
line([avg_pr-sem_pr avg_pr+sem_pr ],[avg_ch avg_ch],'Color','r')
line([avg_pr avg_pr],[avg_ch-sem_ch avg_ch+sem_ch],'Color','r')
xlim([0 1])
ylim([-0.05 0.1])
set(gca,'FontSize',16)
xlabel('Preference before Bolus')
ylabel('Difference in choices (POST-PRE)')

[r,p] = corrcoef(pre_pref,ch);
text(.05,.08,['r=' num2str(r(1,2))],'FontSize',16)
text(.05,.07,['p=' num2str(p(1,2))],'FontSize',16)

%- MAIN FIGURE!!!!!
expl = 2;

%- plot the pref changes 
figure;
for s = expl
    subplot(2,4,[1 2 5 6])
    if bolus_type(s)
        plot(1:length(pre{s}),1-pre{s},'k');hold on;
        plot(length(pre{s})+1:(length(pre{s})*2),1-post{s},'r');
    else
        plot(1:length(pre{s}),pre{s},'k');hold on;
        plot(length(pre{s})+1:(length(pre{s})*2),post{s},'r');
     end
    ylim([0 1]);xlim([0 101])
    line([1 100],[.5 .5],'Color','k')
    ylabel('Preference (Given Juice when > 0.5)');
    xlabel('Trials');
end

% plot choice vs proba difference for all sessions
prob = [10:10:90];
for s = expl
    %- retrieve some variables
    %choice = abs(double(ALL(s).T.choice)-1);
    if bolus_type(s)==1 %- inverse to express in given vs non given juice
        probaG = [ALL(s,1).T.probJ2 ; ALL(s,2).T.probJ2] ;
        probaU = [ALL(s,1).T.probJ1 ; ALL(s,2).T.probJ1];
        choice = [abs(double(ALL(s,1).T.choice)-2) ; abs(double(ALL(s,2).T.choice)-2) ];
    else
        probaG = [ALL(s,1).T.probJ1 ; ALL(s,2).T.probJ1] ;
        probaU = [ALL(s,1).T.probJ2 ; ALL(s,2).T.probJ2];
        choice = [abs(double(ALL(s,1).T.choice)-1) ; abs(double(ALL(s,2).T.choice)-1) ];
    end
    %- unnormalize the proba
    prob_norm = unique(probaG);
    for ii = 1 : length(prob)
        probaG(probaG==prob_norm(ii))=prob(ii);
        probaU(probaU==prob_norm(ii))=prob(ii);
    end
    
    %- Pre/post and quartile groups
    prepost = [zeros(1,size(ALL(s,1).T,1)) ones(1,size(ALL(s,2).T,1)) ];
    quarti = [0:1/(size(ALL(s,1).T,1)-1):1 , 0:1/(size(ALL(s,2).T,1)-1):1 ] ;

    
    bins_pb = -80 : 10 : 80;
    clear probaCh_bins tot_bins nb_bins
    for i = 1 : length(bins_pb)
        samePb =   ((probaG-probaU)==bins_pb(i))';
        %- only consider the last half of the PRE bolus period
        chSamePb = [sum(choice(samePb & prepost==0 & quarti>=0.5)) length(choice(samePb & prepost==0 & quarti>=0.5)) ;
                    sum(choice(samePb & prepost==1 )) length(choice(samePb & prepost==1 )) ];
        
        % chSamePb = [sum(choice(samePb & prepost==0 & quarti<=0.25)) length(choice(samePb & prepost==0 & quarti<=0.25)) ;
        %             sum(choice(samePb & prepost==0 & quarti>=0.75)) length(choice(samePb & prepost==0 & quarti>=0.75)) ;
        %             sum(choice(samePb & prepost==1 & quarti<=0.25)) length(choice(samePb & prepost==1 & quarti<=0.25)) ;
        %             sum(choice(samePb & prepost==1 & quarti>=0.75)) length(choice(samePb & prepost==1 & quarti>=0.75))];
        
        probaCh_bins(:,i) = chSamePb(:,1)./chSamePb(:,2);
        tot_bins(:,i) = chSamePb(:,2);
        nb_bins(:,i) = chSamePb(:,1);
    end
    pre_pref(s) = nanmean(pre{s}(round(length(pre{s})/2):end)); %- only take the 2nd half of the pre

    weights = ones(numel(bins_pb),1);
    coeffs_1 = glmfit(bins_pb, [nb_bins(1,:)', tot_bins(1,:)'], 'binomial','link','logit');
    coeffs_2 = glmfit(bins_pb, [nb_bins(2,:)', tot_bins(2,:)'], 'binomial','link','logit');
    
    % Create a new xAxis with higher resolution
    fineX = min(bins_pb):1:max(bins_pb);
    
    % Generate curve from fit
    curve_1 = glmval(coeffs_1, fineX, 'logit');
    curve_1 = [fineX', curve_1];
    curve_2 = glmval(coeffs_2, fineX, 'logit');
    curve_2 = [fineX', curve_2];
    
    subplot(2,4,[3 4])
    plot(bins_pb,probaCh_bins(1,:),'.','Color',[.4 .4 .4],'MarkerSize',30);hold on
    plot(bins_pb,probaCh_bins(2,:),'.','Color',[255 100 100]/255,'MarkerSize',30)
    plot(curve_1(:,1),curve_1(:,2),'k','LineWidth',2)
    plot(curve_2(:,1),curve_2(:,2),'r','LineWidth',2);
    xlabel('Proba difference (Giv-Ungiv)');
    ylabel('GivenJuice Choice Probability')
    xlim([-80 80])
    subplot(2,4,[7 8])
    plot(curve_1(:,1),curve_2(:,2)-curve_1(:,2),'k','LineWidth',2)
    diff_curve(s,:) = (curve_2(:,2)-curve_1(:,2))';
    diff_proba(s) = nansum(probaCh_bins(1,:)-probaCh_bins(2,:));
    xlim([-80 80])
end
