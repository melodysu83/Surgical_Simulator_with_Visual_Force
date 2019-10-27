% tissue_terrainx = linspace(0,10,1001); % total points 1001 (step size 0.01)
% x = 10*log(1:11)/log(11); % 11 points in the range of 0-10
% tissue_terrainy2=[];
% % make sure upper tissue and lower tissue dont overlap
% while isempty(tissue_terrainy2) || ~isempty(find(tissue_terrainy2 - tissue_terrainy1 < 0.1, 1))...
%         || ~isempty(find(tissue_terrainy1 < 0, 1)) ||  ~isempty(find(tissue_terrainy2 > 10, 1))
%     tmp1 = 0.5*randi(3,1,11)+linspace(1,5,11);
%     tmp2 = 0.5*randi(3,1,11)+linspace(3.5,7.5,11);
%     y1 = min(tmp1,tmp2);
%     y2 = max(tmp1,tmp2);
%     terrain_poly1 = polyfit(x,y1,5);
%     terrain_poly2 = polyfit(x,y2,5); 
%     tissue_terrainy1 = polyval(terrain_poly1,tissue_terrainx);
%     tissue_terrainy2 = polyval(terrain_poly2,tissue_terrainx);
% end
% 
% target_index = [randi(400,1,1) randi(2,1,1)];
% default_tissue5_1 = tissue_terrainy1; % change here: 1-5
% default_tissue5_2 = tissue_terrainy2; % change here: 1-5
% default_target5 = target_index;       % change here: 1-5
% plot(tissue_terrainx,tissue_terrainy1,'r',tissue_terrainx,tissue_terrainy2,'m')

% load to mat file
save('default_tissue_data','default_tissue1_1','default_tissue1_2','default_tissue2_1','default_tissue2_2',...
                           'default_tissue3_1','default_tissue3_2','default_tissue4_1','default_tissue4_2',...
                           'default_tissue5_1','default_tissue5_2','default_target1','default_target2','default_target3','default_target4','default_target5')