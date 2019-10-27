% set default data
tissue1_time = [14.8 18.51 14.88 13.49];
tissue1_force = [0.02582 0.03101 0.02524 0.02735];
tissue1_deform = [0.63456 0.69179 0.61456 0.66868];
tissue1_success = [1 1 1 1];

tissue2_time = [11.0220 8.1510 19.6020 16.5000];
tissue2_force = [0.3021 0.2980 0.2968 0.2083];
tissue2_deform = [2.0080 2.8019 1.4797 1.3184];
tissue2_success = [0 0 0 1];

tissue3_time = [13.2000   14.3880   16.3350   15.0480];
tissue3_force = [0.02582 0.0467    0.0363    0.0346];
tissue3_deform = [1.0928    1.2214    1.1903    1.0184];
tissue3_success = [1 1 1 1];

tissue4_time = [ 20.4270   16.9290   14.0580   19.1730];
tissue4_force = [ 0.2226    0.2982    0.2082    0.2994];
tissue4_deform = [1.3042    1.5494    1.1929    1.4613];
tissue4_success = [1 0 1 0];

tissue5_time = [13.4970   12.6720   15.6420   15.3780];
tissue5_force = [ 0.0386    0.0363    0.0350    0.0365];
tissue5_deform = [ 1.0525    0.9733    0.9703    0.9856];
tissue5_success = [1 1 1 1];

tissuer_time = [17.8860 16.9950 16.9290 13.49];
tissuer_force = [0.0455 0.09941 0.1382 0.02735];
tissuer_deform = [0.8139 1.1414 1.5494 0.66868];
tissuer_success = [1 1 0 1];
load('this_trial_data.mat');
% load new data
if exist('stat_data.mat','file')
    load('stat_data.mat');
    tissue1_count = tissue1_count+[3,3,1,1];
    tissue2_count = tissue2_count+[1,2,0,2];
    tissue3_count = tissue3_count+[3,3,1,1];
    tissue4_count = tissue4_count+[2,3,0,1];
    tissue5_count = tissue5_count+[3,3,1,1];
    tissuer_count = tissuer_count+[2,2,0,1];
    
    % time  max_force  max_deform success
    if ~isempty(tissue1_stat_data)
        tissue1_time = [tissue1_time tissue1_stat_data(:,1)'];
        tissue1_force = [tissue1_force tissue1_stat_data(:,2)'];
        tissue1_deform = [tissue1_deform tissue1_stat_data(:,3)'];
        tissue1_success = [tissue1_success tissue1_stat_data(:,4)'];
    end
    if ~isempty(tissue2_stat_data)
        tissue2_time = [tissue2_time tissue2_stat_data(:,1)'];
        tissue2_force = [tissue2_force tissue2_stat_data(:,2)'];
        tissue2_deform = [tissue2_deform tissue2_stat_data(:,3)'];
        tissue2_success = [tissue2_success tissue2_stat_data(:,4)'];
    end
    if ~isempty(tissue3_stat_data)
        tissue3_time = [tissue3_time tissue3_stat_data(:,1)'];
        tissue3_force = [tissue3_force tissue3_stat_data(:,2)'];
        tissue3_deform = [tissue3_deform tissue3_stat_data(:,3)'];
        tissue3_success = [tissue3_success tissue3_stat_data(:,4)'];
    end
    if ~isempty(tissue4_stat_data)
        tissue4_time = [tissue4_time tissue4_stat_data(:,1)'];
        tissue4_force = [tissue4_force tissue4_stat_data(:,2)'];
        tissue4_deform = [tissue4_deform tissue4_stat_data(:,3)'];
        tissue4_success = [tissue4_success tissue4_stat_data(:,4)'];
    end
    if ~isempty(tissue5_stat_data)
        tissue5_time = [tissue5_time tissue5_stat_data(:,1)'];
        tissue5_force = [tissue5_force tissue5_stat_data(:,2)'];
        tissue5_deform = [tissue5_deform tissue5_stat_data(:,3)'];
        tissue5_success = [tissue5_success tissue5_stat_data(:,4)'];
    end
    if ~isempty(tissuer_stat_data)
        tissuer_time = [tissuer_time tissuer_stat_data(:,1)'];
        tissuer_force = [tissuer_force tissuer_stat_data(:,2)'];
        tissuer_deform = [tissuer_deform tissuer_stat_data(:,3)'];
        tissuer_success = [tissuer_success tissuer_stat_data(:,4)'];
    end       
end


% collect all data
data_time = [tissue1_time';tissue2_time';tissue3_time';tissue4_time';tissue5_time';tissuer_time'];
data_force = [tissue1_force';tissue2_force';tissue3_force';tissue4_force';tissue5_force';tissuer_force'];
data_deform = [tissue1_deform';tissue2_deform';tissue3_deform';tissue4_deform';tissue5_deform';tissuer_deform'];

% data_success = [sum(tissue1_success)/length(tissue1_success) 1-sum(tissue1_success)/length(tissue1_success);...
%                 sum(tissue2_success)/length(tissue2_success) 1-sum(tissue2_success)/length(tissue2_success);...
%                 sum(tissue3_success)/length(tissue3_success) 1-sum(tissue3_success)/length(tissue3_success);...
%                 sum(tissue4_success)/length(tissue4_success) 1-sum(tissue4_success)/length(tissue4_success);...
%                 sum(tissue5_success)/length(tissue5_success) 1-sum(tissue5_success)/length(tissue5_success);...
%                 sum(tissuer_success)/length(tissuer_success) 1-sum(tissuer_success)/length(tissuer_success)]*100;
data_success = [tissue1_count(1) tissue1_count(2)-tissue1_count(1);...
                tissue1_count(3) tissue1_count(4)-tissue1_count(3);...
                NaN NaN;
                tissue2_count(1) tissue2_count(2)-tissue2_count(1);...
                tissue2_count(3) tissue2_count(4)-tissue2_count(3);...
                NaN NaN;
                tissue3_count(1) tissue3_count(2)-tissue3_count(1);...
                tissue3_count(3) tissue3_count(4)-tissue3_count(3);...
                NaN NaN;
                tissue4_count(1) tissue4_count(2)-tissue4_count(1);...
                tissue4_count(3) tissue4_count(4)-tissue4_count(3);...
                NaN NaN;
                tissue5_count(1) tissue5_count(2)-tissue5_count(1);...
                tissue5_count(3) tissue5_count(4)-tissue5_count(3);...
                NaN NaN;
                tissuer_count(1) tissuer_count(2)-tissuer_count(1);...
                tissuer_count(3) tissuer_count(4)-tissuer_count(3)];
            
% plot the statistical result
pos = [1,2,3,4,5,6];
color = ['c', 'y', 'g', 'k'];
group = [repmat({'tissue1'}, length(tissue1_time), 1); repmat({'tissue2'}, length(tissue2_time), 1); ...
         repmat({'tissue3'}, length(tissue3_time), 1); repmat({'tissue4'}, length(tissue4_time), 1); ...
         repmat({'tissue5'}, length(tissue5_time), 1); repmat({'random'}, length(tissuer_time), 1)];
FIGURE = figure('position', [50, 50, 1200, 580]);
set(FIGURE,'name','The Statistical Result','numbertitle','off');

subplot(3,2,1);
boxplot(data_time,group,'positions',pos);
hold on;
plot([0,7],[TimerArrayAll(length(TimerArrayAll)) TimerArrayAll(length(TimerArrayAll))],'--r','linewidth',2);
hold off;
set(gca,'xtick',pos,'xticklabel',{'tissue1','tissue2','tissue3','tissue4','tissue5','random'});
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(1),'FaceAlpha',.5);
end
if tissue_option~=1
    patch(get(h(8-tissue_option),'XData'),get(h(8-tissue_option),'YData'),color(4));
    set(h(8-tissue_option),'color',color(1),'linewidth',2);
else
    patch(get(h(1),'XData'),get(h(1),'YData'),color(4));
    set(h(1),'color',color(1),'linewidth',2);
end
title('Operation Time');
ylabel('time (sec)');


subplot(3,2,3)
boxplot(data_force,group,'positions',pos);
hold on;
plot([0,7],[max(ForceArrayAll) max(ForceArrayAll)],'--r','linewidth',2);
hold off;
set(gca,'xtick',pos,'xticklabel',{'tissue1','tissue2','tissue3','tissue4','tissue5','random'});
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(2),'FaceAlpha',.5);
end
if tissue_option~=1
    patch(get(h(8-tissue_option),'XData'),get(h(8-tissue_option),'YData'),color(4));
    set(h(8-tissue_option),'color',color(2),'linewidth',2);
else
    patch(get(h(1),'XData'),get(h(1),'YData'),color(4));
    set(h(1),'color',color(2),'linewidth',2);
end

title('Maximum Applied Force');
ylabel('Force (N)');

subplot(3,2,5)
boxplot(data_deform,group,'positions',pos);
hold on;
plot([0,7],[max(DeformMaxAll) max(DeformMaxAll)],'--r','linewidth',2);
hold off;
set(gca,'xtick',pos,'xticklabel',{'tissue1','tissue2','tissue3','tissue4','tissue5','random'});
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(3),'FaceAlpha',.5);
end
if tissue_option~=1
    patch(get(h(8-tissue_option),'XData'),get(h(8-tissue_option),'YData'),color(4));
    set(h(8-tissue_option),'color',color(3),'linewidth',2);
    TISSUE_STRING = ['Tissue Geometry :  Default Tissue' num2str(tissue_option-1)];
else
    patch(get(h(1),'XData'),get(h(1),'YData'),color(4));
    set(h(1),'color',color(3),'linewidth',2);
    TISSUE_STRING = 'Tissue Geometry :  Random';
end

title('Maximum Tissue Deform');
ylabel('Deformation (cm)');

subplot(3,2,[2 4])
plot(TimerArrayAll,XArrayAll*10,'-',TimerArrayAll,YArrayAll*10,'-',TimerArrayAll,ForceArrayAll*1000,'-',...
    TimerArrayAll,CollisionArrayAll,'-','LineWidth',2);
hold on;
plot([0 30],[0,0],'k-',[0 0],[0,300],'k-');
hold off;
title('Surgical Simulation Log');
legend('x trajectory (mm)','y trajectory (mm)','max force (mN)','collisions (count)','location','northeastoutside');
xlabel('time (sec)');
vec_pos = get(get(gca, 'XLabel'), 'Position');
set(get(gca, 'XLabel'), 'Position', vec_pos + [20.0 17.0 0]);
if COMPLETE_STATUS
    COMPLETE_STRING = 'Completion Status:  success!';
else
    COMPLETE_STRING = 'Completion Status:  Fail.';
end
switch tissue_option
    case 1
        my_tissue_count = tissuer_count;
    case 2
        my_tissue_count = tissue1_count;
    case 3
        my_tissue_count = tissue2_count;
    case 4
        my_tissue_count = tissue3_count;
    case 5
        my_tissue_count = tissue4_count;
    case 6
        my_tissue_count = tissue5_count;
end
if Force_See
    FORCE_STRING = 'Force Visualization: Enabled.';
else
    FORCE_STRING = 'Force Visualization: Disabled.';
end
descr1 = {'Summary of Surgical Simulation';
         '----------------------------------------------------';
         TISSUE_STRING; 
         ['Time Spent :          ' num2str(max(TimerArrayAll)) ' seconds'];
         ['Max Force :           ' num2str(max(ForceArrayAll)) ' N'];
         ['Max Collisions :      ' num2str(max(CollisionArrayAll)) ' points'];
         COMPLETE_STRING};

descr2 = {'Success Rate';
         '---------------------------------';
         TISSUE_STRING; 
         FORCE_STRING;
         ['With Visual Force: ' num2str(my_tissue_count(1)) '/' num2str(my_tissue_count(2)) '   (' num2str(round(my_tissue_count(1)*10000/my_tissue_count(2))/100) '%)'];
         ['No Visual Force:   ' num2str(my_tissue_count(3)) '/' num2str(my_tissue_count(4)) '   (' num2str(round(my_tissue_count(3)*10000/my_tissue_count(4))/100) '%)']};
text(32,150*max(max(CollisionArrayAll),300)/300,descr1)
text(40,-160*max(max(CollisionArrayAll),300)/300,descr2)

subplot(3,2,6)
bar(data_success,'stacked');
legend('succ.','fail','location','northeastoutside');
set(gca,'xtick',[1.5 4.5 7.5 10.5 13.5 16.5],'xticklabel',{'tissue1','tissue2','tissue3','tissue4','tissue5','random'});
title('Performance (left / right = with / without visual force)');
ylabel('count (# of times)');
% data=randn(10,5);
% data(:,3)=data(:,3)+2;
% 
% boxplot(data) 


% x = [tissue1_time tissue1_force tissue1_deform tissue2_time tissue2_force tissue2_deform...
%      tissue3_time tissue3_force tissue3_deform tissue4_time tissue4_force tissue4_deform...
%      tissue5_time tissue5_force tissue5_deform tissuer_time tissuer_force tissuer_deform];
% 
% group = [ 1*ones(1,length(tissue1_time)), 2*ones(1,length(tissue1_time)), 3*ones(1,length(tissue1_time)),...
%           4*ones(1,length(tissue2_time)), 5*ones(1,length(tissue2_time)), 6*ones(1,length(tissue2_time)),...
%           7*ones(1,length(tissue3_time)), 8*ones(1,length(tissue3_time)), 9*ones(1,length(tissue3_time)),...
%          10*ones(1,length(tissue4_time)),11*ones(1,length(tissue4_time)),12*ones(1,length(tissue4_time)),...
%          13*ones(1,length(tissue5_time)),14*ones(1,length(tissue5_time)),15*ones(1,length(tissue5_time)),...
%          16*ones(1,length(tissuer_time)),17*ones(1,length(tissuer_time)),18*ones(1,length(tissuer_time))];
% positions = [1 1.2 1.4 2 2.2 2.4 3 3.2 3.4 4 4.2 4.4 5 5.2 5.4 6 6.2 6.4];
% boxplot(time, group, 'positions', positions);
% set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9)) mean(positions(10:12))...
%                  mean(positions(13:15)) mean(positions(16:18))]);
% set(gca,'xticklabel',{'tissue1','tissue2','tissue3','tissue4','tissue5','random'})
% 
% color = ['c', 'y', 'g', 'c', 'y', 'g', 'c', 'y', 'g', 'c', 'y', 'g', 'c', 'y', 'g', 'c', 'y', 'g'];
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%    patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
% end
% 
% c = get(gca, 'Children');
% 
% hleg1 = legend(c(1:3), 'operation time', 'max force','max deform');
clear variables