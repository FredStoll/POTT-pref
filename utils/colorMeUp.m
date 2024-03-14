function [colorsArea,colorsArea_sub] = colorMeUp(dataset)

if strcmp(dataset,'POTT')
    order = [3 1 2 5 4 7];
    colorsArea = cbrewer('qual', 'Set2', 8);
    colorsArea = colorsArea(order,:);
    colorsArea_sub = cbrewer('qual', 'Pastel2', 8);
    colorsArea_sub = colorsArea_sub(order,:);
end