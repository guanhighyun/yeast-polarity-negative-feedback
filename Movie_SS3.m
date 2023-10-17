filename = 'FigureData/Sample_file_suddenly_reduce_positive_feedback.xyz';
nframes = 201; % maximum frames for the movie
moviedir = 'Movie_S3'; mkdir(moviedir);
moviename = 'Movie_S3';

[t,positions]=read_molPos3(filename,nframes);
movObj = VideoWriter(sprintf('%s/%s',moviedir,moviename),'MPEG-4');
movObj.FrameRate = 10;
movObj.Quality = 50;
open(movObj)
figure('units','pixels','position',[0 0 1000 500])

% Skip the first 500 seconds of callibration
for i=51:nframes
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
    xlabel(sprintf('%i min',floor(t(i)/60)),'fontsize',25);
    set(gca,'fontsize',20)
    axis square
    
    % drop lambda3 at 600 seconds
    if t(i) >= 600
        lambda3 = [500*ones(1,61),95*ones(1,i-61)];
    else
        lambda3 = [500*ones(1,i)];
    end 
    
    subplot(1,2,2);
    plot((0:10:(i-1)*10)/60,lambda3,'k','linewidth',5);
    title('Strength of positive feedback')
    ylim([50,500]);
    xlim([500/60,2000/60]);
    xlabel('Time (min)')
    ylabel('\lambda_3(s^{-1})')
    set(gca,'fontsize',20)
    axis square
    writeVideo(movObj,getframe(gcf));
    drawnow;
    
end
close(movObj);