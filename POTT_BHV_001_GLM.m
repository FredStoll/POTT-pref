%% Main script for BEHAV analysis of the POTT dataset
%- 
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2022.12

clear

%- locate the files
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/';
% path2go = 'C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\';

skip = true; %- re-run the models or just run on saved models! 

if ~skip
    mk = 'Mimic'; % 'Morbier'

    %- list all the files
    if strcmp(mk(1:2),'Mi')
        list = dir([path2go 'X*a_behav.mat']);
    else
        list = dir([path2go 'M*a_behav.mat']);
    end

    perm = false;
    nPerm = 50;
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

        [ALL(sess,1),param] = behav_model(filename,false,showfig); %- no permutations

        if perm
            for p = 1 : nPerm
                [ALLperm(sess,p),~] = behav_model(filename,perm,showfig);
            end
        end

    end
    if perm
        save([path2go mk '_behav_norm_ALL_prevJuice_ft.mat'],'ALL','ALLperm','param','-v7.3')
    else
        save([path2go mk '_behav_norm_ALL_prevJuice_ft.mat'],'ALL','param','-v7.3')
    end
end

%% POTT - load models for behavioral analyses and Figures (1B/4A-C/5)

M1 = load([path2go 'Morbier_behav_bins.mat'],'ALL','param')
M2 = load([path2go 'Mimic_behav_bins.mat'],'ALL','param')

%- remove sessions before threshold for Morbier, if needed
subset = {'052318' '010122'} ; %- for anything after last change
if ~isempty(subset)
    lim = datenum(subset,'mmddyy');
    for ss = 1 : length(M1.ALL)
        days(ss) = datenum(M1.ALL(ss).name(end-16:end-11),'mmddyy');
    end
    takeme = (days>=lim(1) & days<lim(2)) ;
else
    takeme = true(1,length(M1.ALL));
end
M1.ALL = M1.ALL(takeme);
M1.param = M1.param(takeme);

%- number trials per sessions (all set, not only converging ones)
x=0;
nTr_tot=[];
for m = 1 : 2
    eval(['ALL = M' num2str(m) '.ALL;']);
    for s = 1 : length(ALL)
        x = x + 1;
        nTr_tot(x,:) = [size(ALL(s).T,1) m];
    end
end
[min(nTr_tot(nTr_tot(:,2)==1,1)) median(nTr_tot(nTr_tot(:,2)==1,1)) max(nTr_tot(nTr_tot(:,2)==1,1)) ]
[min(nTr_tot(nTr_tot(:,2)==2,1)) median(nTr_tot(nTr_tot(:,2)==2,1)) max(nTr_tot(nTr_tot(:,2)==2,1)) ]

all_p_trend = [M1.ALL(:).p_trend M2.ALL(:).p_trend ];
all_converge = [M1.ALL(:).converge M2.ALL(:).converge ];
all_p_trend(all_converge==0)=[];
[h_sig, crit_p, adj_p]=fdr_bh(all_p_trend,0.05,'pdep','yes');
thr_trend = crit_p;

%- save a matrix with the names of the considered sessions + converge model
clear behav_sess
x=0;
for m = 1 : 2
    clear conv
    eval(['ALL = M' num2str(m) '.ALL;']);
   
    takeme = [ALL(:).converge]==1;
    sessions = cat(1,ALL(takeme).name);
    for s = 1 : size(sessions,1)
        x = x + 1;
        behav_sess{x,1} = sessions(s,end-17:end-10);
    end
end

% save('Sessions.mat','behav_sess')


%% Fig 1B - Behav model R squared + Estimates
figure;
colors = cbrewer('qual', 'Paired', 10);

subplot(1,6,1)
%- plot adjusted R2
for m = 1 : 2
    clear converge r_adj  keep
    eval(['ALL = M' num2str(m) '.ALL;']);
    for i = 1 : length(ALL)
        r_adj(:,i) = ALL(i).mdl.Rsquared.Adjusted;
        converge(:,i) = ALL(i).converge;
    end
    keep = converge==1 ;
    if m == 1
        col = colors(3:4,:);
        mm = -0.2;
    else
        col = colors(9:10,:);
        mm = .2;
    end
    X= r_adj(1,keep);
    yl=1+mm;
    wdth = .35;
    boxplot_ind(X,yl,wdth,col)
end

set(gca,'view',[90 -90],'color','none','FontSize',16);
set(gca,'YTick',[])
xlabel('R-adj')

xlim([0 1])

subplot(1,6,3:6)
%- plot the estimates
nbSigSess = [];
for m = 1 : 2
    clear converge tStat pVal Estimate r_adj LogLik nTr keep 
    eval(['ALL = M' num2str(m) '.ALL;']);
    for i = 1 : length(ALL)
        tStat(:,i) = ALL(i).mdl.Coefficients.tStat;
        pVal(:,i) = ALL(i).mdl.Coefficients.pValue;
        Estimate(:,i) = ALL(i).mdl.Coefficients.Estimate;
        r_adj(:,i) = ALL(i).mdl.Rsquared.Adjusted;
        LogLik(:,i) = ALL(i).mdl.LogLikelihood;
        converge(:,i) = ALL(i).converge;
        nTr(:,i) = height(ALL(i).T);
    end
    
    keep = converge==1  ;

    pVal = pVal(:,keep);
    for pp = 1 : length(pVal(:,1))
        [h_sig, crit_p, adj_p]=fdr_bh(pVal(pp,:),0.05,'pdep','yes');
        nbSigSess(m,pp) = sum(h_sig) ;
        nbTotSess(m,pp) = length(h_sig) ;
    end

    if m == 1
        col = colors(3:4,:);
        mm = -0.2;
        line([0 0],[0 6],'Color','k');hold on
    else
        col = colors(9:10,:);
        mm = .2;
    end
    for i  = 1 : length( ALL(1).mdl.CoefficientNames)-1
        X= Estimate(i+1,keep);
        yl=i+mm;
        wdth = .35;
        boxplot_ind(X,yl,wdth,col)
    end
    set(gca,'view',[90 -90],'color','none','FontSize',16);
    set(gca,'YTick',1:length(ALL(1).mdl.CoefficientNames)-1,'YTickLabel',ALL(1).mdl.CoefficientNames(2:end),'YTickLabelRotation',25)
    disp([sum(keep) length(keep)])
    if m == 1
        text(19.5,4,['mk M'],'FontSize',16,'Color',col(2,:))
    else
        text(18,4,['mk X'],'FontSize',16,'Color',col(2,:))
    end
end
xlabel('Estimates')
xlim([-20 20])
ylim([0 length( ALL(1).mdl.CoefficientNames)])

nbSigSess(:,2:end)./nbTotSess(:,2:end)
nbSigSess(:,2:end)

%% look at choice consistency in same juice trials for the different preferences
both_pref = [];
figure;
for m = 1 : 2
    subplot(1,2,m)
    clear converge tStat pVal Estimate r_adj LogLik nTr keep
    eval(['ALL = M' num2str(m) '.ALL;']);
    pref = [ALL(:).pref];
    p_trend = [ALL(:).p_trend];
    p_trend_bins = [ALL(:).p_trend_bins];
    choice_consistency_J1 = [ALL(:).choice_consistency_J1];
    choice_consistency_J2 = [ALL(:).choice_consistency_J2];
    converge = [ALL(:).converge];

    plot(pref,choice_consistency_J1,'o'); hold on
    plot(pref,choice_consistency_J2,'or');
    histogram(pref,[0:0.15:1],'Normalization','Probability')
    xlim([0 1])

    %- proportion of sessions with trend in residuals
    disp(['Monkey ' num2str(m)])
    disp('     nSig      nTot     %sig')
    disp([sum(p_trend(converge==1)<=thr_trend) length(p_trend(converge==1)) 100*(sum(p_trend(converge==1)<=thr_trend)/length(p_trend(converge==1)))])

    both_pref = [both_pref , pref];
end

%% plot full correlation matrix (Fig 5C)
clear corrmat_both
for m = 1 : 2
    clear corrmat_all
    eval(['ALL = M' num2str(m) '.ALL;']);

    converge = [ALL(:).converge]==1;

    for d = 1 : size(ALL,1)
        corrmat_all(d,:,:)= ALL(d).corrmat;
    end
    dumm = squeeze(nanmean(abs(corrmat_all(converge,:,:))));

    if m == 1 ; corrmat_both(m,:,:) = dumm;
    elseif m ==2 ; corrmat_both(m,:,:) = dumm';
    end
end
corrmat_both = squeeze(nanmean(corrmat_both));

color_corr = cbrewer('seq', 'Greys', 64);
figure;
h = imagesc(corrmat_both,[0 0.6]);colormap(color_corr);colorbar;
for p1 = 1 : length(corrmat_both)
    for p2 = 1 : length(corrmat_both)
        if p1>p2
            text(p2,p1,['\it' num2str(round(100*corrmat_both(p1,p2))/100)],'HorizontalAlignment','center','Color',colors(4,:),'FontWeight','bold','FontSize',12)
        elseif p1<p2
            text(p2,p1,['\it' num2str(round(100*corrmat_both(p1,p2))/100)],'HorizontalAlignment','center','Color',colors(10,:),'FontWeight','bold','FontSize',12)
        end
    end
end
set(gca,'XTick',1:length(corrmat_both),'XTickLabel',M1.param(1).mat_name,'XTickLabelRotation',25,...
        'YTick',1:length(corrmat_both),'YTickLabel',M1.param(1).mat_name,'YTickLabelRotation',25,...
        'FontSize',14);
hold on;
line([length(corrmat_both)+1 0],[length(corrmat_both)+1 0],'Color',[0.6 0.6 0.6],'LineWidth',2)
h = pixelgrid(h);
h.Children(1).Color = [0.6 0.6 0.6];
%h.Children(2).Color = [31 120 180]/255;


%% Correlation btw Zval and preference - Fig 5B

%- https://en.wikipedia.org/wiki/Riemann_sum
figure;
for m = 1:2
    subplot(2,1,m);
    line([0 0],[-15 15],'Color','k');
    hold on
    line([-0.5 0.5],[0 0],'Color','k'); box on;
    % ref_length = 50; %- use 50 when bin size 30
    ref_length = 10;
    colors = cbrewer('qual', 'Paired', 10);

    clear devia_ref Z_trend p_trend converge
    eval(['ALL = M' num2str(m) '.ALL;']);

    converge = [ALL(:).converge]==1;
    p_trend = [ALL(:).p_trend];

    for d = 1 : size(ALL,1)
        Z_trend(d) = ALL(d).Z_trend;
        devia_ref(d) = mean(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)));
    end

    %- remove session that did not converge
    devia_ref=devia_ref(converge);
    p_trend=p_trend(converge);
    Z_trend=Z_trend(converge);

    if m == 1
        col= colors(3:4,:);
    else
        col= colors(9:10,:);
    end

    plot(devia_ref(p_trend>thr_trend),Z_trend(p_trend>thr_trend),'.','Color',col(1,:),'MarkerSize',15);hold on
    plot(devia_ref(p_trend<=thr_trend),Z_trend(p_trend<=thr_trend),'.','Color',col(2,:),'MarkerSize',15);hold on
    sum(p_trend<=thr_trend)
    xx = devia_ref;
    yy = Z_trend;
    ww = ones(size(Z_trend)); %- not weighted!
    [R,p] = corrcoef(xx,yy);
    slope = R(2,1) * ( std(repelem(yy,ww)) / std(repelem(xx,ww)) );
    intercept = mean(repelem(yy,ww)) - ( slope * mean(repelem(xx,ww))); %- or that, the same thing : sum(W.*yy') - ( slope *  sum(W.*xx'));
    corrline = (xx * slope) + intercept;
    plot(xx,corrline,'-','LineWidth',2,'Color',col(2,:))

    ylim([-15 15]);
    %xlim([-.35 .35]);
    ylabel('Z-value (Trend in residuals)')
    xlabel('Integral of Juice preference ')
    set(gca,'FontSize',12)
    text(-.325,-8.5,['R=' num2str(round(100*R(2,1))/100)],'FontSize',12,'Color',col(2,:))
    text(-.325,-10.5,['p=' num2str(p(2,1))],'FontSize',12,'Color',col(2,:))
end

%% Other correlation with RT and IT (Supp fig)

figure;
for m = 1:2
    subplot(2,1,m);
    line([0 0],[-15 15],'Color','k');
    hold on
    line([-0.5 0.5],[0 0],'Color','k'); box on;
    colors = cbrewer('qual', 'Paired', 10);

    clear devia_ref Z_trend p_trend converge
    eval(['ALL = M' num2str(m) '.ALL;']);

    converge = [ALL(:).converge]==1;
    p_trend = [ALL(:).p_trend];

    for d = 1 : size(ALL,1)
        Z_trend(d) = ALL(d).Z_trend;
        devia_ref(d) = mean(ALL(d).ft_bins-mean(ALL(d).ft_bins(1:ref_length)));
    end

    %- remove session that did not converge
    devia_ref=devia_ref(converge);
    p_trend=p_trend(converge);
    Z_trend=Z_trend(converge);

    if m == 1
        col= colors(3:4,:);
    else
        col= colors(9:10,:);
    end

    plot(devia_ref(p_trend>thr_trend),Z_trend(p_trend>thr_trend),'.','Color',col(1,:),'MarkerSize',15);hold on
    plot(devia_ref(p_trend<=thr_trend),Z_trend(p_trend<=thr_trend),'.','Color',col(2,:),'MarkerSize',15);hold on

    xx = devia_ref;
    yy = Z_trend;
    ww = ones(size(Z_trend)); %- not weighted!
    [R,p] = corrcoef(xx,yy);
    slope = R(2,1) * ( std(repelem(yy,ww)) / std(repelem(xx,ww)) );
    intercept = mean(repelem(yy,ww)) - ( slope * mean(repelem(xx,ww))); %- or that, the same thing : sum(W.*yy') - ( slope *  sum(W.*xx'));
    corrline = (xx * slope) + intercept;
    plot(xx,corrline,'-','LineWidth',2,'Color',col(2,:))

    ylim([-15 15]);
    xlim([-.2 .2]);
    ylabel('Z-value (Trend in residuals)')
    xlabel('Integral of IT ')
    set(gca,'FontSize',12)
    text(-.125,-8.5,['R=' num2str(round(100*R(2,1))/100)],'FontSize',12,'Color',col(2,:))
    text(-.125,-10.5,['p=' num2str(p(2,1))],'FontSize',12,'Color',col(2,:))


end

figure;
for m = 1:2
    subplot(2,1,m);
    line([0 0],[-15 15],'Color','k');
    hold on
    line([-0.5 0.5],[0 0],'Color','k'); box on;
    colors = cbrewer('qual', 'Paired', 10);

    clear devia_ref Z_trend p_trend converge
    eval(['ALL = M' num2str(m) '.ALL;']);

    converge = [ALL(:).converge]==1;
    p_trend = [ALL(:).p_trend];

    for d = 1 : size(ALL,1)
        Z_trend(d) = ALL(d).Z_trend;
        devia_ref(d) = mean(ALL(d).rt_bins-mean(ALL(d).rt_bins(1:ref_length)));
    end

    %- remove session that did not converge
    devia_ref=devia_ref(converge);
    p_trend=p_trend(converge);
    Z_trend=Z_trend(converge);

    if m == 1
        col= colors(3:4,:);
    else
        col= colors(9:10,:);
    end

    plot(devia_ref(p_trend>thr_trend),Z_trend(p_trend>thr_trend),'.','Color',col(1,:),'MarkerSize',15);hold on
    plot(devia_ref(p_trend<=thr_trend),Z_trend(p_trend<=thr_trend),'.','Color',col(2,:),'MarkerSize',15);hold on

    xx = devia_ref;
    yy = Z_trend;
    ww = ones(size(Z_trend)); %- not weighted!
    [R,p] = corrcoef(xx,yy);
    slope = R(2,1) * ( std(repelem(yy,ww)) / std(repelem(xx,ww)) );
    intercept = mean(repelem(yy,ww)) - ( slope * mean(repelem(xx,ww))); %- or that, the same thing : sum(W.*yy') - ( slope *  sum(W.*xx'));
    corrline = (xx * slope) + intercept;
    plot(xx,corrline,'-','LineWidth',2,'Color',col(2,:))

    ylim([-15 15]);
    xlim([-.2 .2]);
    ylabel('Z-value (Trend in residuals)')
    xlabel('Integral of RT ')
    set(gca,'FontSize',12)
    text(-.125,-8.5,['R=' num2str(round(100*R(2,1))/100)],'FontSize',12,'Color',col(2,:))
    text(-.125,-10.5,['p=' num2str(p(2,1))],'FontSize',12,'Color',col(2,:))
end

figure;
for m = 1:2
    subplot(2,1,m);
    line([0 0],[-15 15],'Color','k');
    hold on
    line([-0.5 0.5],[0 0],'Color','k'); box on;
    colors = cbrewer('qual', 'Paired', 10);

    clear devia_ref Z_trend p_trend converge
    eval(['ALL = M' num2str(m) '.ALL;']);

    converge = [ALL(:).converge]==1;
    p_trend = [ALL(:).p_trend];

    for d = 1 : size(ALL,1)
        Z_trend(d) = mean(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)));
        devia_ref(d) = mean(ALL(d).ft_bins-mean(ALL(d).ft_bins(1:ref_length)));
    end

    %- remove session that did not converge
    devia_ref=devia_ref(converge);
    p_trend=p_trend(converge);
    Z_trend=Z_trend(converge);

    if m == 1
        col= colors(3:4,:);
    else
        col= colors(9:10,:);
    end

    plot(devia_ref(p_trend>thr_trend),Z_trend(p_trend>thr_trend),'.','Color',col(1,:),'MarkerSize',15);hold on
    plot(devia_ref(p_trend<=thr_trend),Z_trend(p_trend<=thr_trend),'.','Color',col(2,:),'MarkerSize',15);hold on

    xx = devia_ref;
    yy = Z_trend;
    ww = ones(size(Z_trend)); %- not weighted!
    [R,p] = corrcoef(xx,yy);
    slope = R(2,1) * ( std(repelem(yy,ww)) / std(repelem(xx,ww)) );
    intercept = mean(repelem(yy,ww)) - ( slope * mean(repelem(xx,ww))); %- or that, the same thing : sum(W.*yy') - ( slope *  sum(W.*xx'));
    corrline = (xx * slope) + intercept;
    plot(xx,corrline,'-','LineWidth',2,'Color',col(2,:))

    ylim([-.35 .35]);
    xlim([-.2 .2]);
    ylabel('Integral of Juice Pref')
    xlabel('Integral of IT ')
    set(gca,'FontSize',12)
    text(-.125,-0.2,['R=' num2str(round(100*R(2,1))/100)],'FontSize',12,'Color',col(2,:))
    text(-.125,-0.3,['p=' num2str(p(2,1))],'FontSize',12,'Color',col(2,:))
end

figure;
for m = 1:2
    subplot(2,1,m);
    line([0 0],[-15 15],'Color','k');
    hold on
    line([-0.5 0.5],[0 0],'Color','k'); box on;
    colors = cbrewer('qual', 'Paired', 10);

    clear devia_ref Z_trend p_trend converge
    eval(['ALL = M' num2str(m) '.ALL;']);

    converge = [ALL(:).converge]==1;
    p_trend = [ALL(:).p_trend];

    for d = 1 : size(ALL,1)
        Z_trend(d) = mean(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)));
        devia_ref(d) = mean(ALL(d).rt_bins-mean(ALL(d).rt_bins(1:ref_length)));
    end

    %- remove session that did not converge
    devia_ref=devia_ref(converge);
    p_trend=p_trend(converge);
    Z_trend=Z_trend(converge);

    if m == 1
        col= colors(3:4,:);
    else
        col= colors(9:10,:);
    end

    plot(devia_ref(p_trend>thr_trend),Z_trend(p_trend>thr_trend),'.','Color',col(1,:),'MarkerSize',15);hold on
    plot(devia_ref(p_trend<=thr_trend),Z_trend(p_trend<=thr_trend),'.','Color',col(2,:),'MarkerSize',15);hold on

    xx = devia_ref;
    yy = Z_trend;
    ww = ones(size(Z_trend)); %- not weighted!
    [R,p] = corrcoef(xx,yy);
    slope = R(2,1) * ( std(repelem(yy,ww)) / std(repelem(xx,ww)) );
    intercept = mean(repelem(yy,ww)) - ( slope * mean(repelem(xx,ww))); %- or that, the same thing : sum(W.*yy') - ( slope *  sum(W.*xx'));
    corrline = (xx * slope) + intercept;
    plot(xx,corrline,'-','LineWidth',2,'Color',col(2,:))

    ylim([-.35 .35]);
    xlim([-.2 .2]);
    ylabel('Integral of Juice Pref')
    xlabel('Integral of RT ')
    set(gca,'FontSize',12)
    text(-.125,-0.2,['R=' num2str(round(100*R(2,1))/100)],'FontSize',12,'Color',col(2,:))
    text(-.125,-0.3,['p=' num2str(p(2,1))],'FontSize',12,'Color',col(2,:))
end

%% Example session with change in preference (Mimic) - Fig. 5A

ALL = M2.ALL;

figure

d = 12; %145 15
colors = cbrewer('div', 'RdBu', 128);
norm = @(data) -1+((data-min(data))*2)/(max(data)-min(data)) ;

probas = norm([10:1:90]);

clear sensitivity_bins pref_bins radj_bins ft_bins rt_bins mdl_bins resid
x = 0;
nb_bins = 3;
dumm = round(height(ALL(d).T)/nb_bins);
bins = [(1 : dumm : height(ALL(d).T)) height(ALL(d).T)]
for b = 1 : nb_bins
    disp(b)
    warning off
    mdl_bins = fitglm(ALL(d).T(bins(b):bins(b+1),:),M1.param(1).modelspec,'Distribution','binomial','Link','logit');
    ft = nanmedian(table2array(ALL(d).T(:,'ft')));
    rt = nanmedian(table2array(ALL(d).T(:,'rt')));

    for i = 1 : length(probas)
        tab = table(probas(i)*ones(length(probas),1),probas',nanmedian(ft)*ones(length(probas),1),nanmedian(rt)*ones(length(probas),1),-ones(length(probas),1),-ones(length(probas),1),'VariableNames',{'probJ1','probJ2','ft','rt','prevJuice','prevJuiceReceived'});
        tab.prevJuice = categorical(tab.prevJuice);
        tab.prevJuiceReceived = categorical(tab.prevJuiceReceived);
        if i == 1
            tab_all = tab;
        else
            tab_all = vertcat(tab_all,tab);
        end
    end
    newf_bins = [];
    try
        [newf_bins , ~] = predict(mdl_bins,tab_all);
        newf_bins = reshape(newf_bins,length(probas),length(probas));
    end
    x = x + 1;
    subplot(4,nb_bins,x)
    imagesc([10:1:90],[10:1:90],newf_bins);axis xy;colormap(colors)
end
subplot(4,nb_bins,[x+1 : 4*nb_bins])

bar(ALL(d).pref_bins,'FaceColor',[.8 .8 .8],'BaseValue',mean(ALL(d).pref_bins(1:ref_length)))
hold on
plot(ALL(d).pref_bins,'Color','k','LineWidth',1)

ylim([0 1]);xlim([0 length(ALL(d).pref_bins)])
title(num2str(d))
line([0 length(ALL(d).pref_bins)],[.5 .5],'Color','k')
line([0 length(ALL(d).pref_bins)],[mean(ALL(d).pref_bins(1:ref_length)) mean(ALL(d).pref_bins(1:ref_length))],'Color','r')
for i = 1 : length(bins)-2
    line([bins(i+1) bins(i+1)],[0 1],'Color',[.7 .7 .7])
end
xlabel('Trials')
ylabel('Juice Preference')
set(gca,'FontSize',14)

Z_trend(d)
devia_ref(d)

%% Example pref for Fig 4A

ALL = M1.ALL;

figure
dd = randperm(length(ALL),9); %145 15
x = 0;
for d = 1 : length(dd)
    if ALL(dd(d)).converge
        colors = cbrewer('div', 'RdBu', 128);
        norm = @(data) -1+((data-min(data))*2)/(max(data)-min(data)) ;

        probas = norm([10:1:90]);

        clear sensitivity_bins pref_bins radj_bins ft_bins rt_bins mdl_bins resid

        warning off
        mdl_bins = fitglm(ALL(dd(d)).T,M1.param(1).modelspec,'Distribution','binomial','Link','logit');
        ft = nanmedian(table2array(ALL(dd(d)).T(:,'ft')));
        rt = nanmedian(table2array(ALL(dd(d)).T(:,'rt')));

        for i = 1 : length(probas)
            tab = table(probas(i)*ones(length(probas),1),probas',nanmedian(ft)*ones(length(probas),1),nanmedian(rt)*ones(length(probas),1),-ones(length(probas),1),-ones(length(probas),1),'VariableNames',{'probJ1','probJ2','ft','rt','prevJuice','prevJuiceReceived'});
            tab.prevJuice = categorical(tab.prevJuice);
            tab.prevJuiceReceived = categorical(tab.prevJuiceReceived);
            if i == 1
                tab_all = tab;
            else
                tab_all = vertcat(tab_all,tab);
            end
        end
        newf_bins = [];
        try
            [newf_bins , ~] = predict(mdl_bins,tab_all);
            newf_bins = reshape(newf_bins,length(probas),length(probas));
        end
        x = x + 1;
        subplot(3,3,x)
        imagesc([10:1:90],[10:1:90],newf_bins);axis xy;colormap(colors)
        title([num2str(dd(d)) ' ' num2str(ALL(dd(d)).pref)])
    end
end

%% TRANSITIVITY (FIG 4B)

all_pref = [M1.ALL(:).pref M2.ALL(:).pref ]';
all_converge = [M1.ALL(:).converge M2.ALL(:).converge ]';
all_mk = [repmat({'M'},length(M1.ALL),1) ; repmat({'X'},length(M2.ALL),1)];
all_days =cell(size(all_pref));
x=0;
for m = 1 : 2
    eval(['ALL = M' num2str(m) '.ALL;'])
    for i = 1 : length(ALL)
        x = x + 1;
        all_days{x,:} = [ALL(i).name(end-17:end-11) 'a'];
    end
end

modeldata = table(all_pref,all_converge,all_mk,a ll_days,'VariableNames',{'pref' 'converge' 'mk' 'session'})
modeldata(modeldata.converge==0,:)=[];


%- load juice names
load('/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/Juice_pairs.mat')
rmv = ismember(juice_pairs(:,2),'st');
juice_pairs(rmv,:)=[];

juices = {};
for i = 1 : length(modeldata.session)

    t = find( ismember(juice_pairs(:,1),modeldata.session(i)) );
    if length(t)> 1
        disp(['check ' juice_pairs(t,2)'])
        juices(i,1) = juice_pairs(t(1),2);
    elseif length(t)==0
        disp(['missing session' modeldata.session(i)])
        juices(i,1) = {'XX'};
    else
        juices(i,1) = juice_pairs(t,2);
    end
end

pairs_all = {'AC' 'AG' 'AO' 'AP' 'CG' 'CO' 'CP' 'GO' 'GP' 'OP'}
pref_pairs = [];
for i = 1 : length(pairs_all)
    take = ismember(juices,pairs_all(i));
    take_inv = ismember(juices,pairs_all{i}([2 1]));
    %- the 'unique' is not great here....
    pref_pairs{i,1} = [unique(modeldata.pref(take & ismember(modeldata.mk,'M'))) ; unique(1-modeldata.pref(take_inv  & ismember(modeldata.mk,'M')))];
    pref_pairs{i,2} = [unique(modeldata.pref(take & ismember(modeldata.mk,'X'))) ; unique(1-modeldata.pref(take_inv  & ismember(modeldata.mk,'X')))];
end

figure;
colo = cbrewer('qual', 'Set2', 15);
for m = 1 : 2
    mm = 0;
    subplot(1,2,m)
    for i = 1 : length(pairs_all)
        yl=i+mm;
        if ~isempty(pref_pairs{i,m})

            X = pref_pairs{i,m};
            wdth = .5;
            boxplot_ind(X,yl,wdth,[colo(i,:) ; colo(i,:)])
        end
    end
    set(gca,'view',[90 -90],'color','none','FontSize',16);
    set(gca,'YTick',1:length(pairs_all),'YTickLabel',pairs_all,'YTickLabelRotation',30)
    ylim([0 length(pairs_all)+1]);
    xlim([0 1])
    line([.5 .5],[0 length(pairs_all)+1],'Color','k')
end

colo = cbrewer('qual', 'Paired', 10);
colo1 = colo([3 9],:);
colo2 = colo([4 10],:);
figure;
P_all =[];
datamat=[]
nkeep = 40;
xval_all =[];
yval_all =[];
for m = 1 : 2

    ind_all = 'ACGOP'
    if m==1
        all_pps = {'ACO' 'ACP' }
    else
        all_pps = {'ACO' 'ACP' 'ACG' 'AGO' 'AGP' 'AOP' 'CGO' 'CGP' 'COP' 'GOP'}
    end

    for pa = 1 : length(all_pps)
        pp = all_pps{pa};
        dum = [{pp(1:2)} {pp(2:3)} {pp([1 3])} ];
        if ~isempty(ismember(pairs_all,dum(1)))
            p1 = pref_pairs{ismember(pairs_all,dum(1)),m};
        else
            p1 = 1-pref_pairs{ismember(pairs_all,dum{1}([2 1])),m};
        end
        if ~isempty(ismember(pairs_all,dum(2)))
            p2 = pref_pairs{ismember(pairs_all,dum(2)),m};
        else
            p2 = 1-pref_pairs{ismember(pairs_all,dum{2}([2 1])),m};
        end
        if ~isempty(ismember(pairs_all,dum(3)))
            p3 = pref_pairs{ismember(pairs_all,dum(3)),m};
        else
            p3 = 1-pref_pairs{ismember(pairs_all,dum{3}([2 1])),m};
        end
        nbPermut = max([length(p1) length(p2) length(p3) ])

        xval = [];
        yval = [];
        xx = 0;
        for i = 1 : length(p1)
            for i2 = 1 : length(p2)
                for i3 = 1 : length(p3)
                    p1sel = p1(i);
                    p2sel = p2(i2);
                    p3sel = p3(i3);
                    r1 = p1sel/(1-p1sel);
                    r2 = p2sel/(1-p2sel);
                    r3 = p3sel/(1-p3sel);
                    xx = xx + 1;
                    xval(xx) = r1*r2;
                    yval(xx) = r3;
                end
            end
        end
        xval_all = [xval_all , mean(xval)];
        yval_all = [yval_all , mean(yval)];

        %  datamat = [datamat ; xval' , yval' , repmat(m,size(xval))' , repmat(pa,size(xval))'];

        %  if length(xval)>nkeep
        %      take = randperm(length(xval),nkeep);
        %  else
        %      take = 1 : length(xval);
        %  end
        %  [P,H,STATS] = signtest(xval(take),yval(take),'Method','approximate');
        [P,H,STATS] = signtest(xval,yval,'Method','approximate');
        P_all = [P_all; P STATS.zval m pa];

        %plot(xval,yval,'.','Color',colo(pa,:));hold on
        plot(mean(xval),mean(yval),'.','Color',colo1(m,:),'MarkerSize',20);hold on
        line([mean(xval)-std(xval) mean(xval)+std(xval)],[mean(yval) mean(yval)],'Color',colo1(m,:),'LineWidth',1);hold on
        line([mean(xval) mean(xval)],[mean(yval)-std(yval) mean(yval)+std(yval)],'Color',colo1(m,:),'LineWidth',1);hold on
    end
    plot(mean(xval_all),mean(yval_all),'.','Color',colo2(m,:),'MarkerSize',30);hold on
    line([mean(xval_all)-std(xval_all) mean(xval_all)+std(xval_all)],[mean(yval_all) mean(yval_all)],'Color',colo2(m,:),'LineWidth',1);hold on
    line([mean(xval_all) mean(xval_all)],[mean(yval_all)-std(yval_all) mean(yval_all)+std(yval_all)],'Color',colo2(m,:),'LineWidth',1);hold on

    %- stat test (TO DO!!!!!!)
    % [P,H,STATS] = signtest(xval_all-yval_all);
    % P_all = [P_all; P m];

end
line([-1 4],[-1 4],'Color','k')
xlim([-.5 4]);ylim([-.5 4])
set(gca,'color','none','FontSize',16);
P_all

[P,H,STATS] = signtest(xval_all,yval_all,'Method','approximate')

%- show juice preference for both monkey (Fig 4A)
colors = cbrewer('qual', 'Paired', 10);
figure;
for m = 1 : 2
    clear converge tStat pVal Estimate r_adj LogLik nTr keep
    eval(['ALL = M' num2str(m) '.ALL;']);
    pref = [ALL(:).pref];
    converge = [ALL(:).converge];
    pref(converge==0)=[];
    if m == 1
        col = colors(3:4,:);
        mm = -0.2;
    else
        col = colors(9:10,:);
        mm = .2;
    end
    yl=1+mm;
    wdth = .35;
    boxplot_ind(pref,yl,wdth,col);hold on
    set(gca,'view',[90 -90],'color','none','FontSize',16);

end
xlim([0 1])
set(gca,'YTick',[])
xlabel('Juice Preference')


%% some extra stuff - not needed 

figure;
for m = 1 : 2
    clear corrmat_all p_trend
    eval(['ALL = M' num2str(m) '.ALL;']);
    for d = 1 : size(ALL,1)
        corrmat_all(d,:,:)= ALL(d).corrmat;
        p_trend(d) = ALL(d).p_trend;
    end
    res_corr = squeeze(corrmat_all(:,2,1));
    
    if m == 1
        col = colors(3:4,:);
        mm = -0.2;
        line([0 0],[0 6],'Color','k');hold on
    else
        col = colors(9:10,:);
        mm = .2;
    end
    
        ns = abs(res_corr(p_trend>thr_trend,1));
        sig = abs(res_corr(p_trend<=thr_trend,1));
        [H, pValue(m), KSstat(m)] = kstest2(sig,ns);

        
        wdth = .35;
        boxplot_ind(ns(~isnan(ns)),1+mm,wdth,col)
        boxplot_ind(sig(~isnan(sig)),1+1+mm,wdth,col)
    set(gca,'view',[90 -90],'color','none','FontSize',16);
    set(gca,'YTick',1:2,'YTickLabel',{'ns' 'Trend<0.05'},'YTickLabelRotation',25)
    xlabel('Correlation Juice Pref and smoothed Residuals')
    if m == 1
        text(.1,1.25,['mk M - p(KS)=' num2str(round(100000*pValue(m))/100000)],'FontSize',16,'Color',col(2,:))
    else
        text(.05,1.25,['mk X - p(KS)=' num2str(round(100000*pValue(m))/100000)],'FontSize',16,'Color',col(2,:))
    end
end
ylim([0 3])

figure;
dd = [4 10];
for m = 1 : 2
    clear corrmat_all p_trend Z_trend
    eval(['ALL = M' num2str(m) '.ALL;']);
    for d = 1 : size(ALL,1)
        corrmat_all(d,:,:)= ALL(d).corrmat;
        p_trend(d) = ALL(d).p_trend;
        Z_trend(d) = ALL(d).Z_trend;
    end
    res_corr = squeeze(corrmat_all(:,2,1));
    [r,p]=corrcoef(abs(Z_trend),abs(res_corr))
    plot(abs(Z_trend),abs(res_corr),'.','Color',colors(dd(m),:)); hold on
end
xlabel('Z trend');ylabel('Corr Juice Pref vs Residuals')


%% extract period of preference
changed_all = [];
changed_all_venn = [];

figure;
thr_sd_all = [1 2 2.5 3 3.5 4 4.5 5];
binsz_fact_all = [3 4 5 6 7 8 9 10];
for bb = 3%1 : length(binsz_fact_all)
    for th = 3%1 : length(thr_sd_all)
        all_idxs = [];
        for m = 2
            clear corrmat_all p_trend Z_trend
            eval(['ALL = M' num2str(m) '.ALL;']);
            changed=[];

            for d = 1 : size(ALL,1)
                if ALL(d).converge

                    binsz = floor(length(ALL(d).pref_bins)/binsz_fact_all(bb));

                    ref = ALL(d).pref_bins(1:binsz);

                    pref_zsc = (ALL(d).pref_bins - mean(ref))/std(ref);
                    pref_zsc_smo = smooth(pref_zsc,10);
                    pref_zsc_smo = pref_zsc';
                    [idx,idxs_p] = findenough(pref_zsc_smo',thr_sd_all(th),binsz,'>=');
                    [idx,idxs_n] = findenough(pref_zsc_smo',-thr_sd_all(th),binsz,'<=');
                    if ~isempty(idxs_n) | ~isempty(idxs_p)
                        changed(d,:)=[true ALL(d).p_trend<0.05];
                        all_idxs = [all_idxs , idxs_p idxs_n];
                    else
                        changed(d,:)=[false ALL(d).p_trend<0.05];
                    end
                    subplot(10,20,d);plot(ALL(d).pref_bins);hold on ;
                    plot(idxs_p,ALL(d).pref_bins(idxs_p),'.')
                    plot(idxs_n,ALL(d).pref_bins(idxs_n),'.')
                    ylim([0 1])
                else
                    changed(d,:)=[NaN NaN];
                end
            end

        end
        changed_all = [changed_all ; nanmean(changed)];
        changed_all_venn = [changed_all_venn ; sum(sum(changed,2)==2)/sum(changed(:,2)==true) sum(changed(:,1)==true)/sum(changed(:,2)==false)]
    end
end
plot(changed_all)
figure;plot(changed_all_venn)
figure;histogram(all_idxs) % important plot : show that bins with < or >Xsd are mostly at the end (although the 10 requirement makes it drop by defaults on the last bins)


%% extract period of preference
changed_all = [];
changed_all_venn = [];

thr_sd = 2.5 ;

for m = 1 : 2
    ttt=[];vvv=[];uuu=[];
    all_idxs = [];
    x = 0;pref_ch={};
    clear corrmat_all p_trend Z_trend
    eval(['ALL = M' num2str(m) '.ALL;']);
    eval(['param = M' num2str(m) '.param;']);
    converg = [ALL(:).converge];
    ALL(~converg)=[];
    param(~converg)=[];
    changed=[];
figure;
    for d = 1 : size(ALL,1)
        if ALL(d).converge

            binsz = floor(length(ALL(d).pref_bins)/5);
            ref = ALL(d).pref_bins(1:binsz);

            pref_zsc = (ALL(d).pref_bins - mean(ref))/std(ref);
            [idx,idxs_p] = findenough(pref_zsc,thr_sd,binsz,'>=');
            [idx,idxs_n] = findenough(pref_zsc,-thr_sd,binsz,'<=');
            if ~isempty(idxs_n) | ~isempty(idxs_p)
                changed(d,:)=[true ALL(d).p_trend<0.05];

                %- for the very few sessions where pref goes both direction
                %(~2 across both monkeys), keep the strongest avg change
                if ~isempty(idxs_n) & ~isempty(idxs_p)
                    if mean(pref_zsc(idxs_p)) > abs(mean(pref_zsc(idxs_n)))
                        idxs_n=[];
                    else
                        idxs_p=[];
                    end
                end
                all_idxs = [all_idxs , idxs_p idxs_n];

                %- reconstruct the trials consider for each bins (start/end)
                bins_tr = round(1:(height(ALL(d).T)-param(d).binSize)/50:height(ALL(d).T)-param(d).binSize);
                bins_tr(2,:) = bins_tr(1,:) + param(d).binSize;

                %- find trials for reference period
                temp = bins_tr(:,1:binsz);
                tr_temp = [];
                for i = 1 : length(temp)
                    tr_temp =  [tr_temp temp(1,i):temp(2,i)];
                end
                trref_temp = unique(tr_temp);

                if ~isempty(idxs_n)
                    pref_dir = [mean(ALL(d).pref_bins(1:binsz)) mean(ALL(d).pref_bins(idxs_n))];
                    rt_dir = [mean(ALL(d).rt_bins(1:binsz)) mean(ALL(d).rt_bins(idxs_n))];
                    ft_dir = [mean(ALL(d).rt_bins(1:binsz)) mean(ALL(d).ft_bins(idxs_n))];
                    temp = bins_tr(:,idxs_n);
                else
                    pref_dir = [mean(ALL(d).pref_bins(1:binsz)) mean(ALL(d).pref_bins(idxs_p))];
                    rt_dir = [mean(ALL(d).rt_bins(1:binsz)) mean(ALL(d).rt_bins(idxs_p))];
                    ft_dir = [mean(ALL(d).rt_bins(1:binsz)) mean(ALL(d).ft_bins(idxs_p))];
                    temp = bins_tr(:,idxs_p);
                end
                tr_temp = [];
                for i = 1 : length(temp)
                    tr_temp =  [tr_temp temp(1,i):temp(2,i)];
                end
                trpost_temp = unique(tr_temp);
                
                %- check if overlap in trial and remove then from both bins otherwise
                trpost_temp_old = trpost_temp;
                trpost_temp(ismember(trpost_temp,trref_temp))=[];
                trref_temp(ismember(trref_temp,trpost_temp_old))=[];

                %- put all that in matrix!
                x = x + 1;
                pref_ch{x,1} = ALL(d).name(end-17:end-10) ; %- session
                pref_ch{x,2} = trref_temp; %- ref trials
                pref_ch{x,3} = trpost_temp; %- post trials
                pref_ch{x,4} = pref_dir; %- pref in ref bins and post bins
                pref_ch{x,5} = rt_dir; %- pref in ref bins and post bins
                pref_ch{x,6} = ft_dir; %- pref in ref bins and post bins
                ttt(x) = pref_dir(1)-pref_dir(2);
                uuu(x) = rt_dir(1)-rt_dir(2);
                vvv(x) = ft_dir(1)-ft_dir(2);
            else
                changed(d,:)=[false ALL(d).p_trend<0.05];
            end
            %            subplot(5,5,d);plot(pref_zsc);hold on ; plot(pref_zsc_smo)
            %             plot(idxs_p,pref_zsc_smo(idxs_p),'o')
            %             plot(idxs_n,pref_zsc_smo(idxs_n),'o')
            subplot(10,20,d);plot(ALL(d).pref_bins);hold on ;
            plot(idxs_p,ALL(d).pref_bins(idxs_p),'.')
            plot(idxs_n,ALL(d).pref_bins(idxs_n),'.')
         %   ylim([0 1])
        else
            changed(d,:)=[NaN NaN];
        end
    end

all_idxs_mk{m} = all_idxs;
pref_ch_mk{m} = pref_ch;
nb_sess_mk(m) = sum([ALL(:).converge]);
end

colors = cbrewer('qual', 'Paired', 10);
colors = colors([4 10],:)
figure;
for m = 1: 2
    [a,b] = hist(all_idxs_mk{m},[1:1:50]) % important plot : show that bins with < or >Xsd are mostly at the end (although the 10 requirement makes it drop by defaults on the last bins)
    plot(b/50,100*(a/nb_sess_mk(m)),'Color',colors(m,:),'LineWidth',2);hold on
    ylabel('Percent of sessions')
    xlabel('Normalized time in session')
end
set(gca,'FontSize',16)
%figure;histogram(ttt,[-.5:.05:.5])

save('Pref_bins.mat','pref_ch_mk')

for m = 1 : 2
    for s = 1 : length(pref_ch_mk{m}(:,1))
        dumm = sign(pref_ch_mk{m}{s,4}-.5);
        if dumm(1)==dumm(2)
            pref_dirch{m}(s,1) = false;
        else
            pref_dirch{m}(s,1) = true;
        end
    end
end


%% extract period of ft change
changed_all = [];
changed_all_venn = [];

thr_sd = 2.5 ;

for m = 1 : 2
    ttt=[];vvv=[];uuu=[];
    all_idxs = [];
    x = 0;pref_ch={};
    clear corrmat_all p_trend Z_trend
    eval(['ALL = M' num2str(m) '.ALL;']);
    eval(['param = M' num2str(m) '.param;']);
    converg = [ALL(:).converge];
    ALL(~converg)=[];
    param(~converg)=[];
    changed=[];
figure;
    for d = 1 : size(ALL,1)
        if ALL(d).converge

            binsz = floor(length(ALL(d).ft_bins)/5);
            ref = ALL(d).ft_bins(1:binsz);

            pref_zsc = (ALL(d).ft_bins - mean(ref))/std(ref);
            [idx,idxs_p] = findenough(pref_zsc,thr_sd,binsz,'>=');
            [idx,idxs_n] = findenough(pref_zsc,-thr_sd,binsz,'<=');
            if ~isempty(idxs_n) | ~isempty(idxs_p)
                changed(d,:)=[true ALL(d).p_trend<0.05];

                %- for the very few sessions where pref goes both direction
                %(~2 across both monkeys), keep the strongest avg change
                if ~isempty(idxs_n) & ~isempty(idxs_p)
                    if mean(pref_zsc(idxs_p)) > abs(mean(pref_zsc(idxs_n)))
                        idxs_n=[];
                    else
                        idxs_p=[];
                    end
                end
                all_idxs = [all_idxs , idxs_p idxs_n];

                %- reconstruct the trials consider for each bins (start/end)
                bins_tr = round(1:(height(ALL(d).T)-param(d).binSize)/50:height(ALL(d).T)-param(d).binSize);
                bins_tr(2,:) = bins_tr(1,:) + param(d).binSize;

                %- find trials for reference period
                temp = bins_tr(:,1:binsz);
                tr_temp = [];
                for i = 1 : length(temp)
                    tr_temp =  [tr_temp temp(1,i):temp(2,i)];
                end
                trref_temp = unique(tr_temp);

                if ~isempty(idxs_n)
                    pref_dir = [mean(ALL(d).pref_bins(1:binsz)) mean(ALL(d).pref_bins(idxs_n))];
                    rt_dir = [mean(ALL(d).rt_bins(1:binsz)) mean(ALL(d).rt_bins(idxs_n))];
                    ft_dir = [mean(ALL(d).rt_bins(1:binsz)) mean(ALL(d).ft_bins(idxs_n))];
                    temp = bins_tr(:,idxs_n);
                else
                    pref_dir = [mean(ALL(d).pref_bins(1:binsz)) mean(ALL(d).pref_bins(idxs_p))];
                    rt_dir = [mean(ALL(d).rt_bins(1:binsz)) mean(ALL(d).rt_bins(idxs_p))];
                    ft_dir = [mean(ALL(d).rt_bins(1:binsz)) mean(ALL(d).ft_bins(idxs_p))];
                    temp = bins_tr(:,idxs_p);
                end
                tr_temp = [];
                for i = 1 : length(temp)
                    tr_temp =  [tr_temp temp(1,i):temp(2,i)];
                end
                trpost_temp = unique(tr_temp);
                
                %- check if overlap in trial and remove then from both bins otherwise
                trpost_temp_old = trpost_temp;
                trpost_temp(ismember(trpost_temp,trref_temp))=[];
                trref_temp(ismember(trref_temp,trpost_temp_old))=[];

                %- put all that in matrix!
                x = x + 1;
                pref_ch{x,1} = ALL(d).name(end-17:end-10) ; %- session
                pref_ch{x,2} = trref_temp; %- ref trials
                pref_ch{x,3} = trpost_temp; %- post trials
                pref_ch{x,4} = pref_dir; %- pref in ref bins and post bins
                pref_ch{x,5} = rt_dir; %- pref in ref bins and post bins
                pref_ch{x,6} = ft_dir; %- pref in ref bins and post bins
                ttt(x) = pref_dir(1)-pref_dir(2);
                uuu(x) = rt_dir(1)-rt_dir(2);
                vvv(x) = ft_dir(1)-ft_dir(2);
            else
                changed(d,:)=[false ALL(d).p_trend<0.05];
            end
            %            subplot(5,5,d);plot(pref_zsc);hold on ; plot(pref_zsc_smo)
            %             plot(idxs_p,pref_zsc_smo(idxs_p),'o')
            %             plot(idxs_n,pref_zsc_smo(idxs_n),'o')
            subplot(10,20,d);plot(ALL(d).ft_bins);hold on ;
            plot(idxs_p,ALL(d).ft_bins(idxs_p),'.')
            plot(idxs_n,ALL(d).ft_bins(idxs_n),'.')
            %ylim([0 1])
        else
            changed(d,:)=[NaN NaN];
        end
    end

all_idxs_mk{m} = all_idxs;
pref_ch_mk{m} = pref_ch;
nb_sess_mk(m) = sum([ALL(:).converge]);
end

colors = cbrewer('qual', 'Paired', 10);
colors = colors([4 10],:)
figure;
for m = 1: 2
    [a,b] = hist(all_idxs_mk{m},[1:1:50]) % important plot : show that bins with < or >Xsd are mostly at the end (although the 10 requirement makes it drop by defaults on the last bins)
    plot(b/50,100*(a/nb_sess_mk(m)),'Color',colors(m,:),'LineWidth',2);hold on
    ylabel('Percent of sessions')
    xlabel('Normalized time in session')
end
set(gca,'FontSize',16)
%figure;histogram(ttt,[-.5:.05:.5])

save('FT_bins.mat','pref_ch_mk')

