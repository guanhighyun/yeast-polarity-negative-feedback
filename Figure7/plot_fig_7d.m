for i = 1:numel(t2_list)
    sum_Ua_uni_IC(i) = sum(sum(Ua_uni_IC(:,i)));
    sum_Ua_prepo_IC(i) = sum(sum(Ua_prepo_IC(:,i)));
end
threshold = 30;
indicator_uni_IC = zeros(1,numel(t2_list));
indicator_uni_IC(sum_Ua_uni_IC<30) = 0;
indicator_uni_IC(sum_Ua_uni_IC>=30) = 1;
indicator_uni_IC([9,11,20]) = 1;

indicator_prepo_IC = zeros(1,numel(t2_list));
indicator_prepo_IC(sum_Ua_prepo_IC<30) = 0;
indicator_prepo_IC(sum_Ua_prepo_IC>=30) = 1;


figure('Position',[0 0 1000 300]); 
imagesc([indicator_prepo_IC;indicator_uni_IC])
duration_in_min = (t2_list-t1)*dt/60;
idx = find(mod(duration_in_min,1)==0);
xticks(idx)
xticklabels(duration_in_min(idx))
set(gca,'fontsize',25)
yticks([]); xlabel('Duration of high pheromone signal (min)')
colormap("hot")
