%% POTT Pref - Extract unique neuron for Control analyses
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.10
%- Related to: Stoll & Rudebeck, Neuron, 2024

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data-final/neurons/'
list = dir([path2go 'M*a_SPKpool.mat']);
%list = dir([path2go 'X*a_SPKpool.mat']);

minSpacing = 0.075; %- XX mm spacing minimum or neuron not considered! 
days = [];
for ss = 1 : length(list)
    days = [days;datenum(list(ss).name(2:7),'mmddyy') , ss];
end
date_order = sortrows(days);
list = list(date_order(:,2));

warning off
elecs_all = NaN(157,1);
keep_neurons = []; trash_neurons = [];
dist_from_prev=[];
for s = 1 : length(list)
    if strcmp(list(s).name(1),'X') & s == 64 %- transition to 2nd drive so reset!
        elecs_all = NaN(157,1);
    end
    name = list(s).name;
    load([path2go name]); disp(s)

   % ch = neurons_info(:,[1 3]);
   % [ch_unique,bb] = unique(ch(:,1));
   % elecs_all(ch_unique) = ch(bb,2);

    for el = 1 : 157
        if sum(ismember(neurons_info(:,1),el))~=0 
            ch_depth = neurons_info(ismember(neurons_info(:,1),el),[1 3]);
            if abs(abs(elecs_all(el))-abs(ch_depth(1,2)))>=minSpacing | isnan(elecs_all(el))
                keep_neurons = [keep_neurons; neurons(ismember(neurons_info(:,1),el),1:end-4)] ;
                dist_from_prev = [dist_from_prev; abs(ch_depth(1,2))-elecs_all(el)];
                elecs_all(el) = abs(ch_depth(1,2));
            else
                trash_neurons = [trash_neurons; neurons(ismember(neurons_info(:,1),el),1:end-4)] ;
            end

        end
    end
end

size(keep_neurons,1)/(size(keep_neurons,1)+size(trash_neurons,1))
save([path2go 'Morbier_diff_units_only.mat'],'keep_neurons','minSpacing')
%save([path2go 'Mimic_diff_units_only.mat'],'keep_neurons','minSpacing')

%% Number of turns per day/monkeys

clear 
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data-final/neurons/'
list = dir([path2go 'M*a_SPKpool.mat']);
list = dir([path2go 'X*a_SPKpool.mat']);
for i = 1 : length(list)
    days{i} = list(i).name(2:7);
end

d1 = load([path2go 'Mimic_ElecLocations_Drive2_postHisto.mat'])
d2 = load([path2go 'Mimic_ElecLocations_postHisto.mat'])
d3 = load([path2go 'Morbier_ElecLocations_postHisto.mat'])
LOC = [d1.LOC ; d2.LOC ; d3.LOC]
TIME = [d1.TIME ; d2.TIME ; d3.TIME]

for i = 1 : length(list)
    KEEP_LOC{i} = LOC{strcmp(TIME,days(i))};
end


for i = 1 : length(KEEP_LOC) %- skip the first days
    nMoved(i) = sum(KEEP_LOC{i}(:)~=0);
    dumm = KEEP_LOC{i}(KEEP_LOC{i}(:)~=0)
    nMoved_turns(i) = median(dumm(dumm>0));
    nMoved_turns_mean(i) = mean(dumm(dumm>0));
end

[median(nMoved) 1000*(nanmedian(nMoved_turns)/8) 1000*(nanmean(nMoved_turns)/8)]


