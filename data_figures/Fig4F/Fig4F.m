%% Code that generates Figure 4F
%  Aranda-Diaz et al (2024) Assembly of stool-derived bacterial communities follows "early-bird" resource utilization dynamics

clear all;
close all;
clc;

%% Update the file path and figure name accordingly

% Load data and variables

% Diet (SD MD) x Days (0 1 5 14, pre peak res post) x
% Mouse (1 2) x media (bhi [tyg gam ycfa]) x
% replicates (1 2 3) x taxonomic level (6-Genus 7-Species 8-ASV)
load(fullfile('../data/231024_ODvsRichness/datacompiled_171200.mat'));
taxlevel = 5;
load(fullfile('../data/231024_ODvsRichness/alphaBHI.mat'));
load(fullfile('../data/231024_ODvsRichness/finalOD.mat'));

richs = [];
ods = [];
diets = [];
times = [];
for i1 = 1:2
    for i2 = 1:4
        for i3 = 1:2
            for i4 = 1:1
                for i5 = 1:3
                    temp = datacompiled.relative_abundances...%_mapped...
                        {i1,i2,i3,i4,i5,taxlevel};
                    temp = sum(temp>1e-3,1);
                    richs(end+1) = ...
                        median(alphaBHI(i1,i2,i3,i4,i5,end-3:end));
%                         median(temp(end-3:end));
%                         median(alphaBHI(i1,i2,i3,i4,i5,end-3:end));
                    ods(end+1) = ...
                        median(fODBHI(i1,i2,i3,i4,i5,end-3:end));
                    
                    diets(end+1) = i1;
                    times(end+1) = i2;
                end
            end
        end
    end
end

%---Fit for sparsity
od_max = max(ods);
ods = ods;
Sarrs = 0.6:0.01:0.99;
errs = zeros(length(Sarrs),1);
for j = 1:length(Sarrs)
    ods_pred = (1 - Sarrs(j).^richs).*od_max;
    errs(j) = sum((ods - ods_pred).^2);
end
[~,ind] = min(errs);
Sbest = Sarrs(ind)

%---Plot
cmap = [0,130,200;230,25,75;60,180,75;...
    255,225,25;245,130,48;145,30,180;...
    70,240,240;240,50,230;210,245,60;...
    250,190,212;0,128,128;220,190,255;170,110,40;...
    255,250,200;128,0,0;170,255,195;...
    128,128,0;255,215,180;0,0,128;...
    128,128,128;0,0,0]./255;

figthis = figure;
hold on;
ax = gca;
ax.ActivePositionProperty = 'position';

symbols = {'o','^'};
aad_colors = [239,145,189;146,143,196;148,105,52;203,222,151]./255;
for i1 = 1:2
    for i2 = 1:4
        inds = diets==i1 & times==i2;
        xthis = richs(inds);
        ythis = ods(inds);
        plot(xthis,ythis,symbols{i1},'markersize',8,...
            'markeredgecolor','k',...
            'markerfacecolor',aad_colors(i2,:));
    end
end

xtemp = linspace(0,60,1e3);
ytemp = (1 - Sbest.^xtemp).*od_max;
plot(xtemp,ytemp,'-','linewidth',2,...
    'color',cmap(2,:));

xlabel('Richness');
ylabel('Normalized OD');
ylim([0.5,1.4]);
% xlim([0,15]);

set(gca,'linewidth',0.5,'ticklength',[0.015,0.00]);
set(gca,'layer','top');
set(gcf,'units','points','InnerPosition',[100,100,200,160]);
set(gca,'fontsize',10);
xl = get(gca,'XLabel');
set(xl,'FontSize',12);
yl = get(gca,'YLabel');
set(yl,'FontSize',12);
set(gca,'FontName','Arial')
set(gca,'xminortick','off','yminortick','off');
set(gca,'tickdir','out');
set(figthis,'Renderer','painter');

% saveas(figthis,'Fig4F.svg','svg');

