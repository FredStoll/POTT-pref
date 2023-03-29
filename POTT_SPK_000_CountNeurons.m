%- Count number of neurons across analyses (Table S1)

clear
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/';
list_mk1 = dir([path2go 'M*_SPKpool.mat'])
list_mk2 = dir([path2go 'X*_SPKpool.mat'])

all_neurons_area_mk1 =[];
for i = 1 : length(list_mk1)
    load([path2go list_mk1(i).name],'neurons_area')
    all_neurons_area_mk1 = [all_neurons_area_mk1 ; neurons_area ];
end
all_neurons_area_mk2 =[];
for i = 1 : length(list_mk2)
    load([path2go list_mk2(i).name],'neurons_area')
    all_neurons_area_mk2 = [all_neurons_area_mk2 ; neurons_area ];
end

areas = utils_POTT_areas;
area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'AMG' };

for ar = 1 : length(area2test)
    eval(['nbUnits(ar,1) = sum(ismember(all_neurons_area_mk1,areas.' area2test{ar} '));'])
    eval(['nbUnits(ar,2) = sum(ismember(all_neurons_area_mk2,areas.' area2test{ar} '));'])
end
nbUnits(:,3) = nbUnits(:,1)+nbUnits(:,2)