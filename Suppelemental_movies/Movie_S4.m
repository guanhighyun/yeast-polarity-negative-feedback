filename = 'FigureData/Sample_file_gradually_reduce_positive_feedback.xyz';
nframes = 551; % maximum frames for the movie
moviedir = 'Movie_S4'; mkdir(moviedir);
moviename = 'Movie_S4';

[t,positions]=read_molPos3(filename,nframes);
t(1)=0;
movObj = VideoWriter(sprintf('%s/%s',moviedir,moviename),'MPEG-4');
movObj.FrameRate = 10;
movObj.Quality = 50;
open(movObj)
figure('units','pixels','position',[0 0 1000 500])

min_lambda3 = 95;
max_lambda3 = 500;
tstart = 600;
tstop = 5000;
duration = tstop-tstart;
slope = (min_lambda3-max_lambda3)./duration;
intercept = max_lambda3 - tstart.*slope;

tstart_index = tstart/10+1;
for i=tstart_index-10:nframes % 475
    if isnan(t(i))
       break;
    end
    cla;
        
 
    if ~isempty(positions.Cdc42T{i}) && ~isempty(positions.BemGEF42{i})
        all_active_x = [positions.Cdc42T{i}(:,1);positions.BemGEF42{i}(:,1)];
        all_active_y = [positions.Cdc42T{i}(:,2);positions.BemGEF42{i}(:,2)];
    elseif isempty(positions.Cdc42T{i}) && ~isempty(positions.BemGEF42{i})
        all_active_x = positions.BemGEF42{i}(:,1);
        all_active_y = positions.BemGEF42{i}(:,2);
    elseif ~isempty(positions.Cdc42T{i}) && isempty(positions.BemGEF42{i})
        all_active_x = positions.Cdc42T{i}(:,1);
        all_active_y = positions.Cdc42T{i}(:,2); 
    else
        all_active_x = nan;
        all_active_y = nan;
        
    end
   
    subplot(1,2,1);
    plot(all_active_x,all_active_y,'k.')
    axis([0 8.8623 0 8.8623]);
    xticks([]); yticks([])
    title('Cdc42-GTP distribution')
    xlabel(sprintf('%i min',round(t(i)/60)),'fontsize',25);
    set(gca,'fontsize',20)
    axis square
    
    if t(i) <= 600
        k3 = [500*ones(1,i)];
    elseif (t(i) >= 600) && (t(i) <= 5000)
        k3_before_600 = [500*ones(1,61)];
        k3_after_600 = slope.*(610:10:t(i)) + intercept;
        k3 = [k3_before_600,k3_after_600];
    else
        k3_before_600 = [500*ones(1,61)];
        k3_after_600 = slope.*(610:10:5000) + intercept;
        k3_after_5000 = [95*ones(1,i-501)];
        k3 = [k3_before_600,k3_after_600,k3_after_5000];
    end 
    
    subplot(1,2,2);
    plot((0:10:(i-1)*10)/60,k3,'k','linewidth',5);
    title('Strength of positive feedback')
    ylim([50,500]);
    xlim([500/60,(nframes-1)*10/60]);
    xlabel('Time (min)')
    ylabel('\lambda_3(s^{-1})')
    set(gca,'fontsize',20)
    axis square
    writeVideo(movObj,getframe(gcf));
    drawnow;
end
close(movObj);