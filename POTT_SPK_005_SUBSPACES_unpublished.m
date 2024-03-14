%% SUBSPACE ANALYSIS - FIGURE 3

%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.03

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/'; %- can take 'MORBIER' or 'MIMIC'
list = dir([path2go '*a_SPKpool.mat']);

disp(['Computing LDA on ' num2str(length(list)) ' sessions'])
areas = utils_POTT_areas;

area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'AMG' };

minNeurons = 5;
overwrite = 'n';

thresh = 'n';
thr = 1; % 1Hz threshold for neurons
thr_len = 500; % for 200 ms
rmv = 100; % remove xxx ms on each side on every events (avoid overlaps between bins and smoothing problems)
% times_evts = {'FixFP_onset' 'Stim_onset' 'Resp_onset' 'FixResp_onset' 'FB_onset' 'Rew_onset' 'FB_offset'};
bins4decoding=[2]; %- perform decoding on subset on bins (stim and go period here)
window = [200 800]; %- empty bracket if every time points, otherwise will do LDA on average FR during [t_start t_end] in ms

minnorm  = @(data) (data-min(min(data)))/(max(max(data))-min(min(data)));
x=0;
modeldata = table();
for s = 1 : length(list)
    name = list(s).name;

    load([path2go name]); disp(s)

    bins2remove = (rmv/subsp)-1;

    n_evt = length(times_evts);
    lin=length(neurons_rank)*n_evt;
    nTimes = [0 ; (sum(abs(times(:,:)),2)/subsp)];
    col = sum(nTimes);
    bins = NaN(1,col);
    time = NaN(1,col);

    for b = 1 : length(times_evts)
        time(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) =  (times(b,1):subsp:times(b,2)-(subsp/2));
        bins(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) = b;
        bins([sum(nTimes(1:b))+1:sum(nTimes(1:b))+1+bins2remove   sum(nTimes(1:b+1))-bins2remove:sum(nTimes(1:b+1))    ]) = NaN; %- remove 100 ms each side (avoid overlaps between bins and smoothing problems)
    end

    if ~isempty(SPK_INS)

        %- cut out the overlapping bins
        SPK_INS_cut = SPK_INS;
        SPK_INS_cut(:,~ismember(bins,bins4decoding))=[];

        time(~ismember(bins,bins4decoding)) = [];

        if strcmp(thresh,'y')%- reject neurons with FR too low
            reject = [];
            for n = 1 : length(neurons_rank)
                if mean(mean(SPK_INS_cut(Tr_Clust_INS(:,2)==n,:)))<thr  %- put NaNs when FR is below thr
                    reject = [reject , n];
                end
            end

            %- remove them from main matrices
            SPK_INS_cut(   ismember(Tr_Clust_INS(:,2),reject)   ,:) = [];
            Tr_Clust_INS(   ismember(Tr_Clust_INS(:,2),reject)   ,:) = [];
        end

        %- normalize FR
        SPK_INSnorm = SPK_INS;
        units = unique(Tr_Clust_INS(:,2));
        for u = 1 : length(units)
            temp = SPK_INS(Tr_Clust_INS(:,2)==units(u),~isnan(bins)); %- could normalize only on the bin2use, but more likely to get some 0 FR?? so NaNs!
            SPK_INSnorm(Tr_Clust_INS(:,2)==units(u),~isnan(bins)) = reshape(minnorm(temp(:)),size(temp,1),size(temp,2));
        end


        SPK_INS_cut = SPK_INSnorm;
        SPK_INS_cut(:,~ismember(bins,bins4decoding))=[];

        %- create data and factor matrix for the decoding
        data = SPK_INS_cut;
        unit_nb = Tr_Clust_INS(:,2);

        %- trials to considerfor instrum
        when = @(arg,cd) TrialType_INS(strcmp(TrialType_header,arg),:)==cd ;
        diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
        diff_juice_all = repmat(diff_juice',length(unique(unit_nb)),1);

        dumm = TrialType_INS(ismember(TrialType_header,{'nTr' 'I_chosenproba' 'I_chosenside' }),:)';
        % dumm = [dumm(:,1) dumm(:,2)+(100*dumm(:,3))]
        dumm = [dumm(:,1) dumm(:,3)+(100*dumm(:,2))];

        predic_all = repmat(dumm(:,2),length(unique(unit_nb)),1);
        factor = [predic_all , unit_nb];

        %- only keep a subset of proba

        keep = [130 150 170 190 230 250 270 290];
        remove = ~ismember(factor(:,1),keep) | ~diff_juice_all ;
        data(remove,:) = [];
        factor(remove,:)=[];
        unit_nb(remove)=[];

        tr2cons = factor(factor(:,2)==factor(1,2),1);
        for k = 1 : length(keep)
            nTr(k)=sum(tr2cons==keep(k));
        end
        tr2take = min(nTr);

        %  factor(factor(:,1)==10  ,1)=1;
        %  factor(factor(:,1)==90  ,1)=2;
        % data = data(factor(:,1)==1 | factor(:,1)==2 ,:);
        % factor = factor(factor(:,1)==1 | factor(:,1)==2 ,:);

        %   ref_juice = [ true(5,5) false(5,5) ; false(5,5) true(5,5)];
        %   ref_proba = [ diag(true(5,1)) diag(true(5,1)) ; diag(true(5,1)) diag(true(5,1))];
        %   ref_both = diag(true(10,1));

        %- for every area with enough neurons... min is 5
        for ar = 1 : length(area2test)
            eval(['takeit = find(ismember(neurons_area,areas.' area2test{ar} '));'])
            if length(takeit)>= minNeurons

                %% Initialize

                data_sub = nanmean(data(ismember(unit_nb,takeit),time>=window(1) & time<window(2)),2);
                factor_sub = factor(ismember(unit_nb,takeit),:);

                %- remove unwanted probas (10)
                all_pbs = [130 150 170 190 230 250 270 290];
                probas = unique(mod(all_pbs,100));
                remove = ~ismember(factor_sub(:,1),all_pbs);
                data_sub(remove,:) = [];
                factor_sub(remove,:)=[];

                clear units
                for u = 1 : length(takeit)
                    units(:,u) = data_sub(factor_sub(:,2)==takeit(u));
                end
                factor_sub = factor_sub(factor_sub(:,2)==takeit(u),1);

                chPb=mod(factor_sub(:,1),100);
                chJu = floor(factor_sub(:,1)/100);

                %- kfold validation..
                nK = 10;
                %cvPb = cvpartition(chPb,'KFold',nK);
                %cvJu = cvpartition(chJu,'KFold',nK);
                cv = cvpartition(factor_sub,'KFold',nK); %- partition using both juice and proba. but sometimes missing some categories (chose J2 at 30% for ex, when strong J1 pref)
                clear units_pbavg units_juavg angl_PC units_pb1avg units_pb2avg
                x = x + 1;
                %  subplot(4,4,x);

                for k = 1 : nK

                    %- average trials for Pb or Juice // using same cv for both
                    sub = cv.training(k);
                    for p = 1 : length(probas)
                        units_pbavg(p,:) = mean(units(chPb==probas(p) & sub==1,:));
                    end
                    units_pbavg_ctr = units_pbavg  - repmat(mean(units_pbavg),size(units_pbavg,1),1);

                    for p = 1 : 2 %length(juices)
                        units_juavg(p,:) = mean(units(chJu==p & sub==1,:));
                    end
                    units_juavg_ctr = units_juavg  - repmat(mean(units_juavg),size(units_juavg,1),1);

                    units_ctr =  units  - repmat(mean(units),size(units,1),1);

                    PC.nComp = length(takeit) ;
                    [PC.eigenvectors,PC.score,PC.eigenvalues,~,PC.explained,PC.mu] = pca(units_pbavg_ctr,'NumComponents',PC.nComp,'Centered',false);
                    [PC2.eigenvectors,PC2.score,PC2.eigenvalues,~,PC2.explained,PC2.mu] = pca(units_juavg_ctr,'NumComponents',PC.nComp,'Centered',false);

                    angl_PC(k) = rad2deg(subspace(PC.eigenvectors(:,1),PC2.eigenvectors(:,1)));
                    angl_PC_rad(k) = subspace(PC.eigenvectors(:,1),PC2.eigenvectors(:,1));


                    
                                [XX,YY,param_LDA] = pop_init_noRdmTrials(data_centered,factor_sub,time,'perm',1,'minTr',tr2take,'pop',true); % ,'window',[-.5 1.5]

                                if ~isempty(param.window)
                                    XX_avg{1} = mean(cat(3,XX{time>=param.window(1) & time<=param.window(2)}),3);
                                    XX = XX_avg;
                                    YY = YY(1);
                                end






                    %- for example: only use 3 dimensions + no normalization?
%                     col_Pb = cbrewer('seq','Greys',6); col_Pb = col_Pb(2:end,:);
%                     col_Ju = cbrewer('seq','Oranges',6) ; col_Ju = col_Ju([2 4 5],:);
% 
%                     [PC.eigenvectors,PC.score,PC.eigenvalues,~,PC.explained,PC.mu] = pca(units_pbavg_ctr(:,1:3),'NumComponents',3,'Centered',false);
%                     [PC2.eigenvectors,PC2.score,PC2.eigenvalues,~,PC2.explained,PC2.mu] = pca(units_juavg_ctr(:,1:3),'NumComponents',3,'Centered',false);
%                     fact_v = 1;
%                     for cd = 1 : 4
%                         plot3(units_pbavg_ctr(cd,1),units_pbavg_ctr(cd,2),units_pbavg_ctr(cd,3),'.','MarkerSize',15,'Color',col_Pb(cd,:))  ; hold on
%                     end
%                     line([0 PC.eigenvectors(1,1)]*fact_v,[0 PC.eigenvectors(2,1)]*fact_v,[0 PC.eigenvectors(3,1)]*fact_v,...
%                         'Color','k','LineWidth',1);
%                     grid on
% 
%                     plot3(units_juavg_ctr(1,1),units_juavg_ctr(1,2),units_juavg_ctr(1,3),'.','MarkerSize',15,'Color',col_Ju(1,:))
%                     plot3(units_juavg_ctr(2,1),units_juavg_ctr(2,2),units_juavg_ctr(2,3),'.','MarkerSize',15,'Color',col_Ju(2,:))
%                     line([0 PC2.eigenvectors(1,1)]*fact_v,[0 PC2.eigenvectors(2,1)]*fact_v,[0 PC2.eigenvectors(3,1)]*fact_v,...
%                         'Color',col_Ju(end,:),'LineWidth',1);
%                     grid on
%                     lim = max(max(abs([units_pbavg_ctr(:,1:3) ; units_juavg_ctr(:,1:3)])));
%                     lim = lim+(0.1*lim);
%                     xlim([-lim lim]);xlabel('neuron 1')
%                     ylim([-lim lim]);ylabel('neuron 2')
%                     zlim([-lim lim]);zlabel('neuron 3')
%                     title([num2str((angl_PC)) ])

                    %- Same but comparing ChosenProba J1 vs ChosenProba J2
                    %- average trials for Pb j1 and Pb J2 // using same cv for both
                    sub = cv.training(k);
                    for p = 1 : length(probas)
                        units_pb1avg(p,:) = mean(units(chPb==probas(p) & chJu==1 & sub==1,:));
                    end
                    units_pb1avg_ctr = units_pb1avg  - repmat(mean(units_pb1avg),size(units_pb1avg,1),1);

                    sub = cv.training(k);
                    for p = 1 : length(probas)
                        units_pb2avg(p,:) = mean(units(chPb==probas(p) & chJu==2 & sub==1,:));
                    end
                    units_pb2avg_ctr = units_pb2avg  - repmat(mean(units_pb2avg),size(units_pb2avg,1),1);

                    PC.nComp = length(takeit) ;
                    [PC_j1.eigenvectors,PC_j1.score,PC_j1.eigenvalues,~,PC_j1.explained,PC_j1.mu] = pca(units_pb1avg_ctr,'NumComponents',PC.nComp,'Centered',false);
                    [PC_j2.eigenvectors,PC_j2.score,PC_j2.eigenvalues,~,PC_j2.explained,PC_j2.mu] = pca(units_pb2avg_ctr,'NumComponents',PC.nComp,'Centered',false);

                    angl_PC_j12(k) = rad2deg(subspace(PC_j1.eigenvectors(:,1),PC_j2.eigenvectors(:,1)));
                    angl_PC_rad_j12(k) = subspace(PC_j1.eigenvectors(:,1),PC_j2.eigenvectors(:,1));

                end
                % title([num2str(s) ' - ' area2test{ar} ' - ' num2str(mean(angl_PC))])
                % pause ; hold off

                theta(s,ar) = mean(angl_PC);
                theta_j12(s,ar) = mean(angl_PC_j12);
                modeldata = [modeldata ; table(mean(angl_PC),mean(angl_PC_j12),area2test(ar),{name(1)},{name(2:7)},'VariableNames',{'angle' 'angle_j12' 'area' 'mk' 'session'})];

            else
                theta(s,ar) = NaN;
                theta_j12(s,ar) = NaN;
            end
        end
    else
        theta(s,:) = NaN;
        theta_j12(s,:) = NaN;
    end
end

save('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\POTT_Subspace_angle.mat')


%%
clear
load('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\POTT_Subspace_angle.mat')

order = [3 1 2 5 4 7];
colorsArea = cbrewer('qual', 'Set2', 8);
colorsArea = colorsArea(order,:);
colorsArea_sub = cbrewer('qual', 'Pastel2', 8);
colorsArea_sub = colorsArea_sub(order,:);

figure
subplot(1,2,1)
mm = 0;
for i  = 1 : length( area2test)
    X = theta(~isnan(theta(:,i)),i);
    %  X = theta_jr_deg(~isnan(theta_jr_deg(:,i)),i);
    %  X = theta_pr_deg(~isnan(theta_pr_deg(:,i)),i);
    yl=i+mm;
    wdth = .5;
    boxplot_ind(X,yl,wdth,[colorsArea_sub(i,:) ; colorsArea(i,:)])
end
set(gca,'view',[90 -90],'color','none','FontSize',16);
set(gca,'YTick',1:length(area2test),'YTickLabel',area2test,'YTickLabelRotation',45)
ylim([0 length(area2test)+1]);
title('Proba vs Juice subspaces')
xlabel('Angle (degree, 0 = parallel // 90 = orthogonal)')

%- stats
modeldata.angle = double(modeldata.angle);
modeldata.angle_j12 = double(modeldata.angle_j12);
models_form = {'angle ~ 1 + area  + (1|mk) + (1|session)' ; 'angle ~ 1 + area  + (1|mk)'};
%[lme,model_final] = model_comparison(modeldata,models_form,false);
lme = fitglme(modeldata,models_form{1}); model_final = models_form{1};

[~,wald_angNeurons,thr_corr,pval_adj_angNeurons]  = area_posthoc(lme,area2test,'n');
disp(['%%%%%%%%%% Angle of proba/juice with model = ' model_final ' %%%%%%%%%%'])
disp(anova(lme));disp(pval_adj_angNeurons);disp(wald_angNeurons);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%- plot the significance
pval2=[pval_adj_angNeurons , zeros(length(area2test),1)]+ [pval_adj_angNeurons' ; zeros(1,length(area2test))];
for ar1  = 1 : length( area2test)
    Xmax =max(theta(~isnan(theta(:,ar1)),ar1));
    Xmin = min(theta(~isnan(theta(:,ar1)),ar1));
    Xmean = mean(theta(~isnan(theta(:,ar1)),ar1));
    updown=[1.5 2];
    for ar2 = 1 : length(area2test)
        Xmean2 = mean(theta(~isnan(theta(:,ar2)),ar2));
        if ar1~=ar2 & pval2(ar1,ar2)<thr_corr & Xmean<Xmean2

            text(Xmax+(updown(1)),ar1,'*','Color',colorsArea(ar2,:),'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
            updown(1) = updown(1) + 2;
        elseif ar1~=ar2 & pval2(ar1,ar2)<thr_corr & Xmean>Xmean2
            text(Xmin-(updown(2)),ar1,'*','Color',colorsArea(ar2,:),'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
            updown(2) = updown(2) + 2;
        end
    end
end



subplot(1,2,2)
mm = 0;
for i  = 1 : length( area2test)
    X = theta_j12(~isnan(theta_j12(:,i)),i);
    %  X = theta_jr_deg(~isnan(theta_jr_deg(:,i)),i);
    %  X = theta_pr_deg(~isnan(theta_pr_deg(:,i)),i);
    yl=i+mm;
    wdth = .5;
    boxplot_ind(X,yl,wdth,[colorsArea_sub(i,:) ; colorsArea(i,:)])
end
set(gca,'view',[90 -90],'color','none','FontSize',16);
set(gca,'YTick',1:length(area2test),'YTickLabel',area2test,'YTickLabelRotation',45)
ylim([0 length(area2test)+1]);
title('Proba J1 vs Proba J2 subspaces')
xlabel('Angle (degree, 0 = parallel // 90 = orthogonal)')
set(gcf, 'Color', [1 1 1]);


median(modeldata.angle(ismember(modeldata.area,'vlPFC')))
median(modeldata.angle(ismember(modeldata.area,'OFC')))
median(modeldata.angle(ismember(modeldata.area,'LAI')))
median(modeldata.angle(ismember(modeldata.area,'AMG')))
median(modeldata.angle(ismember(modeldata.area,'IFG')))


models_form = {'angle_j12 ~ 1 + area  + (1|mk) + (1|session)' ; 'angle_j12 ~ 1 + area  + (1|mk)'};
%[lme,model_final] = model_comparison(modeldata,models_form,false);
lme = fitglme(modeldata,models_form{1}); model_final = models_form{1};

[~,wald_angpbNeurons,thr_corr,pval_adj_angpbNeurons]  = area_posthoc(lme,area2test,'n');
disp(['%%%%%%%%%% Angle of proba J1/proba J2 with model = ' model_final ' %%%%%%%%%%'])
disp(anova(lme));disp(pval_adj_angpbNeurons);disp(wald_angpbNeurons);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


%- plot the significance
pval2=[pval_adj_angpbNeurons , zeros(length(area2test),1)]+ [pval_adj_angpbNeurons' ; zeros(1,length(area2test))];
for ar1  = 1 : length( area2test)
    Xmax =max(theta_j12(~isnan(theta_j12(:,ar1)),ar1));
    Xmin = min(theta_j12(~isnan(theta_j12(:,ar1)),ar1));
    Xmean = mean(theta_j12(~isnan(theta_j12(:,ar1)),ar1));
    updown=[1.5 2];
    for ar2 = 1 : length(area2test)
        Xmean2 = mean(theta_j12(~isnan(theta_j12(:,ar2)),ar2));
        if ar1~=ar2 & pval2(ar1,ar2)<thr_corr & Xmean<Xmean2

            text(Xmax+(updown(1)),ar1,'*','Color',colorsArea(ar2,:),'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
            updown(1) = updown(1) + 2;
        elseif ar1~=ar2 & pval2(ar1,ar2)<thr_corr & Xmean>Xmean2
            text(Xmin-(updown(2)),ar1,'*','Color',colorsArea(ar2,:),'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
            updown(2) = updown(2) + 2;
        end
    end
end









%%

behav = load('/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/Morbier_behav_norm_ALL_prevJuice_only.mat')
behav2 = load('/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/Mimic_behav_norm_ALL_prevJuice_only.mat')
ALL = [behav.ALL ; behav2.ALL]
% pref_session = [behav.ALL.pref];
clear days
for i = 1 : length(ALL)
    pref_session(i) = ALL(i).pref;
    trend_session(i) = ALL(i).Z_trend_bins;
    days{i,:} = ALL(i).name(end-17:end-11);
end

juice_pref = NaN(length(modeldata.area),1);
for i = 1 : length(modeldata.area)
    sess = ismember(days,[modeldata.mk{i} modeldata.session{i}]);
    juice_pref(i) = pref_session(sess);
end
modeldata = [modeldata , table(juice_pref,'VariableNames',{'pref'})];

figure;
for ar = 1 : length(area2test)
    sub = modeldata(ismember(modeldata.area,area2test(ar)),:)
    subplot(1,length(area2test),ar)
    plot(abs(sub.pref-.5),sub.angle_j12,'o')
end


%%
colorsPr = cbrewer('seq', 'Greys', 100);
colorsArea = cbrewer('qual', 'Paired', length(area2test));

figure;
for ar = 1 : length(area2test)
    sub = modeldata(ismember(modeldata.area,area2test(ar)),:)

    x = abs(sub.pref-.5);
    y=sub.angle_j12;

    [r,p] = corrcoef(x,y)
    allRP(ar,:) = [r(2,1) , p(2,1)];

    xgrid=linspace(0,.5,100);
    ygrid=linspace(0,90,100);
    [x1,y1] = meshgrid(xgrid, ygrid);
    % Perform kernel density estimate
    % [x y] is actual data, xi is the desired grid points to evaluate
    % f is an estimate of the density, ep(:,1) is the X location, ep(:,2) is the y location
    xi = [x1(:) y1(:)];
    [f,ep]=ksdensity([x y],xi); % remove the outputs to see a 3D plot of the distribution
    % format data in matrix for contourf and plot
    X = reshape(ep(:,1),length(xgrid),length(ygrid));
    Y = reshape(ep(:,2),length(xgrid),length(ygrid));
    Z = reshape(f,length(xgrid),length(ygrid));
    subplot(1,length(area2test),ar);
    imagesc(X(1,:),Y(:,1),Z);axis xy
    colormap(colorsPr)
    hold on
    plot(x,y,'.','Color',colorsArea(ar,:))
    text(.1,5,['r=' num2str(round(allRP(ar,1),3))])
    text(.3,5,['p=' num2str(round_pval(allRP(ar,2)))])

    cf = fit(x,y,'poly1'); % fit
    p_handle = plot(cf,'k','predfunc',.95); % Plot fit
    set(p_handle,'Color',colorsArea(ar,:),'LineWidth',2);
    legend off
    xlabel('Absolute Juice Preference');ylabel('Proba J1 vs Proba J2 Subspace angle')
    title(area2test{ar})

end

%% ONLY USE ANGLES WHEN SIGNIFICANT DECODING OF BOTH PARAM


clear
load('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\POTT_Subspace_angle.mat')
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/'
LDA = load([path2go 'res_LDA_stim.mat']);
res_LDA = LDA.res_LDA;

measures = {'I_chosenjuice' 'I_chosenproba'}
area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'AMG'}
order = [3 1 2 5 4 7];
colorsArea = cbrewer('qual', 'Set2', 8);
colorsArea = colorsArea(order,:);
colorsArea_sub = cbrewer('qual', 'Pastel2', 8);
colorsArea_sub = colorsArea_sub(order,:);
% figure;
maxPerf = [];
sessions_all=char();
for m = 1 : length(measures)
    tab(m).modeldata=[];

    for ar = 1 : length(area2test)
        eval(['curr = res_LDA.' measures{m} '.' area2test{ar}])
        ss = size(curr(1).perf);
        perf = [];
        sess=[];
        mk = [];
        for s = 1 : length(curr)
            perf(s,:) = curr(s).perf_pval;

            sess{s,1} = curr(s).lda_sess ;
            mk{s,1} = curr(s).lda_sess(1) ;
        end
        curr_area = repmat(area2test(ar),size(mk));
        if size(perf,2)==1
            tab(m).modeldata = [tab(m).modeldata ; table(perf,curr_area,sess,mk,'VariableNames',{'perf' 'area' 'session' 'mk'})];
        else
            tab(m).modeldata = [tab(m).modeldata ; table(perf(:,1),perf(:,2),curr_area,sess,mk,'VariableNames',{'perf1' 'perf2' 'area' 'session' 'mk'})];
        end
    end
end

sigOrNot=false;
for s = 1 : length(modeldata.session)
    curr_sess = {[modeldata.mk{s} modeldata.session{s} 'a']};
    curr_ar = modeldata.area(s);

    current_1 = find(ismember(tab(1).modeldata.area ,curr_ar{1}) & ismember(tab(1).modeldata.session ,curr_sess{1}));
    current_2 = find(ismember(tab(2).modeldata.area ,curr_ar{1}) & ismember(tab(2).modeldata.session ,curr_sess{1}));
    sigOrNot(s,1) = tab(1).modeldata.perf(current_1)<0.05 & tab(2).modeldata.perf(current_2)<0.05 ;
end


modeldata_sig = modeldata(sigOrNot,:);


figure
subplot(1,2,1)
mm = 0;
for ar  = 1 : length( area2test)
    X = modeldata_sig.angle(ismember(modeldata_sig.area,area2test{ar}))
    yl=ar+mm;
    wdth = .5;
    boxplot_ind(X,yl,wdth,[colorsArea_sub(ar,:) ; colorsArea(ar,:)])
end
set(gca,'view',[90 -90],'color','none','FontSize',16);
set(gca,'YTick',1:length(area2test),'YTickLabel',area2test,'YTickLabelRotation',45)
ylim([0 length(area2test)+1]);
title('Proba vs Juice subspaces')
xlabel('Angle (degree, 0 = parallel // 90 = orthogonal)')

%- stats
modeldata_sig.angle = double(modeldata_sig.angle);
modeldata_sig.angle_j12 = double(modeldata_sig.angle_j12);
models_form = {'angle ~ 1 + area  + (1|mk) + (1|session)' ; 'angle ~ 1 + area  + (1|mk)'};
%[lme,model_final] = model_comparison(modeldata,models_form,false);
lme = fitglme(modeldata_sig,models_form{1}); model_final = models_form{1};

[~,wald_angNeurons,thr_corr,pval_adj_angNeurons]  = area_posthoc(lme,area2test,'y');
disp(['%%%%%%%%%% Angle of proba/juice with model = ' model_final ' %%%%%%%%%%'])
disp(anova(lme));disp(pval_adj_angNeurons);disp(wald_angNeurons);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%- plot the significance
pval2=[pval_adj_angNeurons , zeros(length(area2test),1)]+ [pval_adj_angNeurons' ; zeros(1,length(area2test))];
for ar1  = 1 : length( area2test)
    Xmax =max(modeldata_sig.angle(ismember(modeldata_sig.area,area2test{ar1})));
    Xmin = min(modeldata_sig.angle(ismember(modeldata_sig.area,area2test{ar1})));
    Xmean = mean(modeldata_sig.angle(ismember(modeldata_sig.area,area2test{ar1})));
    updown=[1.5 2];
    for ar2 = 1 : length(area2test)
        Xmean2 = mean(modeldata_sig.angle(ismember(modeldata_sig.area,area2test{ar2})));
        if ar1~=ar2 & pval2(ar1,ar2)<thr_corr & Xmean<Xmean2

            text(Xmax+(updown(1)),ar1,'*','Color',colorsArea(ar2,:),'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
            updown(1) = updown(1) + 2;
        elseif ar1~=ar2 & pval2(ar1,ar2)<thr_corr & Xmean>Xmean2
            text(Xmin-(updown(2)),ar1,'*','Color',colorsArea(ar2,:),'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
            updown(2) = updown(2) + 2;
        end
    end
end


median(modeldata_sig.angle(ismember(modeldata_sig.area,'vlPFC')))
median(modeldata_sig.angle(ismember(modeldata_sig.area,'OFC')))
median(modeldata_sig.angle(ismember(modeldata_sig.area,'LAI')))
median(modeldata_sig.angle(ismember(modeldata_sig.area,'AMG')))
median(modeldata_sig.angle(ismember(modeldata_sig.area,'IFG')))








