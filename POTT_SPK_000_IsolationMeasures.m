%% POTT Pref - Isolation metrics - TABLE S1
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.03
%- Related to: Stoll & Rudebeck, Neuron, 2024

clear
spkpath = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/';
cd(spkpath)

listneurons = dir('*SPKpool.mat');
spk_area = [];
isolation =[];
noise_overlap =[];
peak_snr =[];
peak_noise =[];
fr =[];
monkey =[];
name =[];

for ne = 1 : length(listneurons)
    disp(ne)
    
    spike = load(listneurons(ne).name,'neurons_info','neurons_info_header','neurons');
    load([listneurons(ne).name(1:9) 'AreaInfo.mat']);
    
    isolation = [isolation ; spike.neurons_info(:,5)] ;
    noise_overlap = [noise_overlap ; spike.neurons_info(:,6)] ;
    peak_snr = [peak_snr ; spike.neurons_info(:,9)] ;
    peak_noise = [peak_noise ; spike.neurons_info(:,8)] ;
    fr = [fr ; spike.neurons_info(:,4)] ;
    spk_area = [spk_area ; area_histology] ;
    monkey = [monkey ; repmat({listneurons(ne).name(1)},length(area_histology),1)] ;
    name = [name ; spike.neurons];
end
% save([spkpath 'POTT_Isolation.mat'],'isolation','monkey',"name","spk_area","fr","peak_noise","peak_snr","noise_overlap")

load([spkpath 'POTT_Isolation.mat'],'isolation','monkey',"name","spk_area","fr","peak_noise","peak_snr","noise_overlap")

median([isolation , noise_overlap , peak_snr , peak_noise , fr])
%figure;plot(peak_snr,fr,'o')

areas = utils_POTT_areas;
area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'AMG' };

%- COLOR ASSIGNMENT
order = [3 1 2 5 7 4];
colorsArea = cbrewer('qual', 'Set2', 8);
colorsArea = colorsArea(order,:);
colorsArea_sub = cbrewer('qual', 'Pastel2', 8);
colorsArea_sub = colorsArea_sub(order,:);
mk={'M' 'X'}
figure
for ar = 1 : length(area2test)
    for m = 1 : 2
    eval(['takeit = ismember(spk_area,areas.' area2test{ar} ') & ismember(monkey,mk(m));'])
          nb_unit(ar,m) = sum(takeit);     
    subplot(2,5,1+5*(m-1));set(gca,'FontSize',22)
    [a,b]=hist(isolation(takeit),[round(min(isolation),2):.01:round(max(isolation),2)  ]);
    plot(b,a/sum(a),'Color',colorsArea(ar,:),'lineWidth',2);title('ISOLATION');hold on

    subplot(2,5,2+5*(m-1));set(gca,'FontSize',22)
    [a,b]=hist(noise_overlap(takeit),[round(min(noise_overlap),2):.01:round(max(noise_overlap),2)  ]);
    plot(b,a/sum(a),'Color',colorsArea(ar,:),'lineWidth',2);title('NOISE OVERLAP');hold on

    subplot(2,5,3+5*(m-1));set(gca,'FontSize',22)
    [a,b]=hist(peak_snr(takeit),[round(min(peak_snr),2):1:round(max(peak_snr),2)  ]);
    plot(b,a/sum(a),'Color',colorsArea(ar,:),'lineWidth',2);title('PEAK SNR');hold on

    subplot(2,5,4+5*(m-1));set(gca,'FontSize',22)
    [a,b]=hist(peak_noise(takeit),[round(min(peak_noise),2):.25:round(max(peak_noise),2)  ]);
    plot(b,a/sum(a),'Color',colorsArea(ar,:),'lineWidth',2);title('PEAK NOISE');hold on

    subplot(2,5,5+5*(m-1));set(gca,'FontSize',22)
    [a,b]=hist(fr(takeit),[round(min(fr),2):1:round(max(fr),2)  ]);
    plot(b,a/sum(a),'Color',colorsArea(ar,:),'lineWidth',2);title('FR');hold on
    xlim([0 40]);

    end
end

tableS1=cell(length(area2test)*2,5);
mk={'M' 'X'};
x=0;
for ar = 1 : length(area2test)
    for m = 1 : 2
        eval(['takeit = ismember(spk_area,areas.' area2test{ar} ') & ismember(monkey,mk(m));'])
        nb_unit(ar,m) = sum(takeit);
        x=x+1;
        tableS1(x,1:5) = [{[num2str(round(prctile(isolation(takeit),50),3)) ' [' num2str(round(prctile(isolation(takeit),25),3)) ' ' num2str(round(prctile(isolation(takeit),75),3)) ']' ]} ,...
            {[num2str(round(prctile(noise_overlap(takeit),50),3)) ' [' num2str(round(prctile(noise_overlap(takeit),25),3)) ' ' num2str(round(prctile(noise_overlap(takeit),75),3)) ']' ]} ,...
            {[num2str(round(prctile(peak_snr(takeit),50),2)) ' [' num2str(round(prctile(peak_snr(takeit),25),2)) ' ' num2str(round(prctile(peak_snr(takeit),75),2)) ']' ]} ,...
            {[num2str(round(prctile(peak_noise(takeit),50),2)) ' [' num2str(round(prctile(peak_noise(takeit),25),2)) ' ' num2str(round(prctile(peak_noise(takeit),75),2)) ']' ]} ,...
            {[num2str(round(prctile(fr(takeit),50),2)) ' [' num2str(round(prctile(fr(takeit),25),2)) ' ' num2str(round(prctile(fr(takeit),75),2)) ']' ]}];
    end
end



