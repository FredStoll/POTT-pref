%% Merge neurons to create pseudo populations

clear
% mk = 'BOTH';
% if strcmp(mk,'BOTH')
%         mk = 'MIMIC1'; path2go1 = utils_POTT_SPKfolder(mk); %- can take 'MORBIER' or 'MIMIC'
%         mk = 'MIMIC2'; path2go2 = utils_POTT_SPKfolder(mk); %- can take 'MORBIER' or 'MIMIC'
%         mk = 'MORBIER'; path2go3 = utils_POTT_SPKfolder(mk); %- can take 'MORBIER' or 'MIMIC'
%         list1 = dir([path2go1 'X*a_SPKpool.mat']);
%         list2 = dir([path2go2 'X*a_SPKpool.mat']);
%         list3 = dir([path2go3 'M*a_SPKpool.mat']);
%         path2go_all = [repmat({path2go1},size(list1)) ; repmat({path2go2},size(list2)) ; repmat({path2go3},size(list3))];
%         list = [list1 ;list2;list3];
% end
areas = utils_POTT_areas;
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/'; %- path where SPKpool files are!
%path2go = 'F:\POTT_POOL\subset\'; %- path where SPKpool files are!
list = dir([path2go 'X*a_SPKpool.mat']);
% area2test = {'vlPFCa' 'vlPFCp' 'OFCa' 'OFCp' 'dlPFCv' 'dlPFCd' 'SMA' 'PMd' 'Cd' 'PUT' 'AMG' 'LAI' 'IFG'};
area2test = {'vlPFC' 'OFC' 'IFG' 'AMG' 'LAI' };
%area2test = {'12r' '12m' '12o' '12l' 'a11ml' '13l' '13m' 'LAI'};
area2test = {'vlPFC' };

lastID = zeros(1,length(area2test));

for s = 1 : length(list)

    current = load([path2go list(s).name]);
    load([path2go list(s).name(1:9) 'AreaInfo.mat'])
    disp(s)
    
    for a = 1 : length(area2test)
        
            eval([' neurons2take = find(ismember(area_histology,areas.'  area2test{a} '));']);
        
        if s==1
            SPK(a).INS = [];BEHAV(a).INS = [];
            INFO(a).neurons_info = [];
        end
        
        if ~isempty(neurons2take)
            
            for n = 1 : length(neurons2take)
                lastID(a) = lastID(a)+1;
                INFO(a).neurons_info  = [INFO(a).neurons_info ; [lastID(a) , current.neurons_info(neurons2take(n),:)] ];
                
                if ~isempty(current.SPK_INS)
                    spkINS = current.SPK_INS(current.Tr_Clust_INS(:,2)==neurons2take(n),:);
                    SPK(a).INS = [SPK(a).INS ; spkINS];
                    tr_dummy = ones(1,length(current.TrialType_INS(1,:)));
                    BEHAV(a).INS = [BEHAV(a).INS , [lastID(a)*tr_dummy ; s*tr_dummy ; current.TrialType_INS   ]   ] ;
                end

            end
        end
    end
end


%- extract matrices for each area2test!
for a = 1 : length(area2test)
    SPK_INS = SPK(a).INS;
    
    BEHAV_INS = BEHAV(a).INS;
    
    neurons_info = INFO(a).neurons_info;
    length(neurons_info)
    save([path2go 'POOL_' num2str(length(list)) 'sessions_' area2test{a} '.mat'],'SPK_INS','BEHAV_INS','neurons_info','-v7.3')
end




%% keep only neurons far apart


clear
% mk = 'BOTH';
% if strcmp(mk,'BOTH')
%         mk = 'MIMIC1'; path2go1 = utils_POTT_SPKfolder(mk); %- can take 'MORBIER' or 'MIMIC'
%         mk = 'MIMIC2'; path2go2 = utils_POTT_SPKfolder(mk); %- can take 'MORBIER' or 'MIMIC'
%         mk = 'MORBIER'; path2go3 = utils_POTT_SPKfolder(mk); %- can take 'MORBIER' or 'MIMIC'
%         list1 = dir([path2go1 'X*a_SPKpool.mat']);
%         list2 = dir([path2go2 'X*a_SPKpool.mat']);
%         list3 = dir([path2go3 'M*a_SPKpool.mat']);
%         path2go_all = [repmat({path2go1},size(list1)) ; repmat({path2go2},size(list2)) ; repmat({path2go3},size(list3))];
%         list = [list1 ;list2;list3];
% end
areas = utils_POTT_areas;
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/'; %- path where SPKpool files are!
%path2go = 'F:\POTT_POOL\subset\'; %- path where SPKpool files are!
list = dir([path2go 'X*a_SPKpool.mat']);
% area2test = {'vlPFCa' 'vlPFCp' 'OFCa' 'OFCp' 'dlPFCv' 'dlPFCd' 'SMA' 'PMd' 'Cd' 'PUT' 'AMG' 'LAI' 'IFG'};
area2test = {'vlPFC' 'OFC' 'IFG' 'AMG' 'LAI' };
%area2test = {'12r' '12m' '12o' '12l' 'a11ml' '13l' '13m' 'LAI'};
%area2test = {'IFG' 'AMG' 'LAI'};

M_keep = load([path2go 'Morbier_diff_units_only.mat'])
X_keep = load([path2go 'Mimic_diff_units_only.mat'])
no_duplicates = [M_keep.keep_neurons ; X_keep.keep_neurons ];


lastID = zeros(1,length(area2test));

for s = 1 : length(list)

    current = load([path2go list(s).name]);
    load([path2go list(s).name(1:9) 'AreaInfo.mat'])
    disp(s)
    
     keepme = false;
    for i = 1 : length(current.neurons_area)
        if ~isempty(find(ismember(no_duplicates,{current.neurons(i,1:end-4)})))
            keepme(i) = true;
        else
            keepme(i) = false;
        end
    end
    keepme = keepme';

    for a = 1 : length(area2test)
        
            eval([' neurons2take = find(ismember(area_histology,areas.'  area2test{a} ') & keepme );']);
        
        if s==1
            SPK(a).INS = [];BEHAV(a).INS = [];
            INFO(a).neurons_info = [];
        end
        
        if ~isempty(neurons2take)
            
            for n = 1 : length(neurons2take)
                lastID(a) = lastID(a)+1;
                INFO(a).neurons_info  = [INFO(a).neurons_info ; [lastID(a) , current.neurons_info(neurons2take(n),:)] ];
                
                if ~isempty(current.SPK_INS)
                    spkINS = current.SPK_INS(current.Tr_Clust_INS(:,2)==neurons2take(n),:);
                    SPK(a).INS = [SPK(a).INS ; spkINS];
                    tr_dummy = ones(1,length(current.TrialType_INS(1,:)));
                    BEHAV(a).INS = [BEHAV(a).INS , [lastID(a)*tr_dummy ; s*tr_dummy ; current.TrialType_INS   ]   ] ;
                end

            end
        end
    end
end


%- extract matrices for each area2test!
for a = 1 : length(area2test)
    SPK_INS = SPK(a).INS;
    
    BEHAV_INS = BEHAV(a).INS;
    
    neurons_info = INFO(a).neurons_info;
    length(neurons_info)
    save([path2go 'POOLsub_' num2str(length(list)) 'sessions_' area2test{a} '.mat'],'SPK_INS','BEHAV_INS','neurons_info','-v7.3')
end


