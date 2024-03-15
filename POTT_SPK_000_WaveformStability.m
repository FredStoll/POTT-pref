%% POTT Pref - Extract 'stable' neurons using waveforms.
%-
%- first average across all waveforms and find the peak/trough location, which is kept constant!
%-
%- Neurons can be found 'unstable' but still highly isolated.. Super super
%- conservative method here..!
%- 
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.06
%- Related to: Stoll & Rudebeck, Neuron, 2024

clear
mk2take = {'M' 'X1' 'X2'}; %- X1 and X2 refers to the 2 drives for Monkey X!

for mm = 1 : length(mk2take)

    clearvars -except mm mk2take
    if strcmp(mk2take{mm},'M')
        path2go = 'F:\POTT_MORBIER\SORTED NEURONS\SPK_Files\'
    elseif strcmp(mk2take{mm},'X1')
        path2go = 'F:\POTT_MIMIC\SORTED NEURONS\SPK_Files\'
    elseif strcmp(mk2take{mm},'X2')
        path2go = 'F:\POTT_MIMIC_02\SORTED NEURONS\SPK_Files\'
    end

    list = dir([path2go '*a_Ch*_Clus*.mat']);
    alignment = [30 34];
    ref = 10;
    nb_bins = 50;

    avg_ratio = NaN(length(list),nb_bins);

    parfor n = 1 : length(list)
        disp(n)
        spk = load([path2go list(n).name]);
        nb_spk = 1 : length(spk.waveforms(:,1));
        bins_spk = discretize(nb_spk,quantile(nb_spk,[0:1/nb_bins:1]));

        wf_avg = mean(spk.waveforms);
        minORmax = [min(wf_avg(alignment(1):alignment(2))) max(wf_avg(alignment(1):alignment(2)))];
        trough = minORmax(abs(minORmax)==max(abs(minORmax)));
        trough_loc = alignment(1)+find(wf_avg(alignment(1):alignment(2))==trough)-1;

        minORmax = [min(wf_avg(trough_loc+5:end)) max(wf_avg(trough_loc+5:end))];
        peak = minORmax(sign(minORmax)~=sign(trough));
        if length(peak)>1
            peak = peak(abs(peak)==max(abs(peak)));
        end
        peak_loc = trough_loc+5+find(wf_avg(trough_loc+5:end)==peak)-1;

        if length(spk.waveforms(:,1))>100
            for i = 1 : nb_bins
                wf_avg = mean(spk.waveforms(bins_spk==i,:));
                avg_ratio(n,i) = wf_avg(peak_loc)/wf_avg(trough_loc);
            end

        end
    end

    avg_ratio_norm = (avg_ratio - mean(avg_ratio(:,1:ref),2)) ./   repmat(std(avg_ratio(:,1:ref)')',1,nb_bins);
    figure;plot(avg_ratio_norm')

    %- extract when exceed given threshold
    thr = 2.5;
    avg_ratio_sig = false(size(avg_ratio_norm));
    for n = 1 : length(avg_ratio_norm)
        [idx,idxs] = findenough(avg_ratio_norm(n,:),thr,3,'>');
        if ~isempty(idx)
            avg_ratio_sig(n,idxs)=true;
        end
        [idx,idxs] =findenough(avg_ratio_norm(n,:),-thr,3,'<');
        if ~isempty(idx)
            avg_ratio_sig(n,idxs)=true;
        end
    end
    figure;plot(mean(avg_ratio_sig))

    spk_names = {list(:).name}';

    save(['C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\' mk2take{mm} '_waveforms_ratio.mat'],'spk_names','avg_ratio','avg_ratio_norm','avg_ratio_sig','thr','nb_bins','ref')
end

%% combine all 3 files across monkeys! 

% clear
% load('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\X1_waveforms_ratio.mat','spk_names','avg_ratio','avg_ratio_norm','avg_ratio_sig','thr','nb_bins','ref')
% x = load('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\X2_waveforms_ratio.mat','spk_names','avg_ratio','avg_ratio_norm','avg_ratio_sig','thr','nb_bins','ref')
% 
% spk_names = [spk_names ; x.spk_names];
% avg_ratio = [avg_ratio ; x.avg_ratio];
% avg_ratio_norm = [avg_ratio_norm ; x.avg_ratio_norm];
% avg_ratio_sig = [avg_ratio_sig ; x.avg_ratio_sig];
% 
% save('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\X_waveforms_ratio.mat','spk_names','avg_ratio','avg_ratio_norm','avg_ratio_sig','thr','nb_bins','ref')
% 
% clear
% load('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\M_waveforms_ratio.mat','spk_names','avg_ratio','avg_ratio_norm','avg_ratio_sig','thr','nb_bins','ref')
% x = load('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\X_waveforms_ratio.mat','spk_names','avg_ratio','avg_ratio_norm','avg_ratio_sig','thr','nb_bins','ref')
% 
% spk_names = [spk_names ; x.spk_names];
% avg_ratio = [avg_ratio ; x.avg_ratio];
% avg_ratio_norm = [avg_ratio_norm ; x.avg_ratio_norm];
% avg_ratio_sig = [avg_ratio_sig ; x.avg_ratio_sig];
% 
% save('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\POTT_waveforms_ratio.mat','spk_names','avg_ratio','avg_ratio_norm','avg_ratio_sig','thr','nb_bins','ref')

%% some stats! 

clear
load('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\POTT_waveforms_ratio.mat','spk_names','avg_ratio','avg_ratio_norm','avg_ratio_sig','thr','nb_bins','ref')

for n = 1 : length(spk_names)
    [Z_trend_wf(n),P_wf(n)]=Mann_Kendall(avg_ratio_norm(n,:),0.01);
end
sum(P_wf<0.01)



