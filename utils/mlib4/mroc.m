function [auroc] = mroc(x,y)
% function [auroc] = mroc(x,y)
%
% computes the area under the receiver operating characteristic curve
% x and y must be vectors hold integers!
%
% results validated a) by simulation with normcdf b) by comparison
% with web-based calculators and c) with code in the MES Toolbox
% 
% also see http://www.mathworks.com/matlabcentral/fileexchange/32398
% 
% Maik C. Stüttgen, July 2012

dispcurve = 0;
start_val = min([min(x) min(y)])-1;    % overall smallest value -1
end_val   = max([max(x) max(y)])-1;    % overall largest value +1
c         = start_val:end_val;         % list of criterion values
a         = zeros(numel(c),2);         % initialize HR and FA array

for i = 1:size(c,2);                   % go through all criterion values
  a(i,1) = sum(y>c(i))/length(y);      % p(y>c)
  a(i,2) = sum(x>c(i))/length(x);      % p(x>c)
end

a     = [a;0,0;1,0;1,1];               % close the curve
auroc = polyarea(a(:,1),a(:,2));

if dispcurve
    figure
  title(['AUROC=' num2str(auroc,'%1.2f')]),hold on
  scatter(a(:,1),a(:,2),'.')
  plot(a(:,1),a(:,2))
  xlabel('p(y>c)'),ylabel('p(x>c)')
  axis([0 1 0 1])
end