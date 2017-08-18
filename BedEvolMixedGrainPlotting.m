function [ Grain ] = BedEvolMixedGrainPlotting( varargin )

    %   UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here

    clc
    close all

    load('stationing.mat','Station');    
    load('Velocity.mat','Vsave');
    load('transportrates.mat','Qbsave');
    load('bedelevation.mat','etasave');
    load('slope.mat','Ssave');
    load('grainsize.mat','Dgsave');
    load('grainsizedistribution.mat','Fisave');
    load('sandfraction.mat','Fssave');
    load('standarddeviation.mat','Sigmasave');
    load('depth.mat','ysave');
    load('time.mat','time');
    %load('numberofsaves.mat','NoS');
    load('watersurface.mat','WSEsave');
    load('rmsmean.mat','RMSSize');
    load('d90save.mat','D90save');
    
    [N] = Station(:,1);            
    
    Grain = [0.2,0.4,0.7,1.4,2.8,5.7,11.3,22.6,45.3];
    
    linecolors = winter(11);
    linecolors1 = gray(11);
    linecolors2 = autumn(11);
    N_Reverse = sort(N,'descend') .* 1000;
    
    % Shift topography and water surface to match the flume coordinates
    etasave_shift = (etasave .* 1000) - 800;
    WSEsave_shift = (WSEsave .* 1000) - 800;
    
    % Make a zero-crossing plot with the simulated bed profile
    % Fit a linear function to the trended bed profile
    [curve,~,~] = fit(N_Reverse,etasave_shift(end,:)','poly1');
    % Grab the curve coefficient values
    CC = coeffvalues(curve);
    % Compute the zero crossing fit curve
    ZeroCrossFunc = (N_Reverse .* CC(1)) + CC(2);
    % Plot the results
    hh = figure(2);
    plot(N_Reverse,etasave_shift(end,:),'k-','LineWidth',1.5)
    hold on
    plot(N_Reverse,ZeroCrossFunc,'k--','LineWidth',0.5)
    hold on
    %
    % Do some calculations to find the exact fill area. Do this call a
    % function call Fill_Between_Curves.
    [output] = Fill_Between_Curves(N_Reverse,etasave_shift(end,:),ZeroCrossFunc);
    % Use the output from the Fill_Between_Curves function to make the
    % filled area plot.
    % cmap = viridis(100);
    color = [115/255, 203/255, 255/255];
    xout = output(1,:);
    topout = output(2,:);
    botout = output(3,:);
    fh = fill([xout fliplr(xout)],[botout fliplr(topout)],color);
    hold on
    fh.EdgeColor = 'none';
    uistack(fh,'bottom')
    hold off
    axis([0,16000,150,550])
    ax = gca;
    set(gca,'FontName','Sans','FontSize',16);
    ax.XTick = ([0 2000 4000 6000 8000 10000 12000 14000 16000]);
    ax.XTickLabel = ({'0','2000','4000','6000','8000','10000','12000','14000','16000'});
    ax.XTickLabelRotation = 45;
    ax.YTick = ([200 300 400 500]);
    ax.YTickLabel = ({'200','300','400','500'});
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridLineStyle = '--';
    ax.GridColor = [0.35 0.35 0.35];
    % Set the figure print size
    x0 = 10;
    y0 = 10;
    width = 25;
    height = 5;
    set(gcf, 'Units','centimeters','Position',[x0 y0 width height])
    % set(gcf,'PaperType','B2','PaperPositionMode','auto','PaperOrientation','portrait') 
    % Ask the user to navigate to the directory which contains the two processed DEMs.
    [Directory] = uigetdir('C:\Users\schartrand\Dropbox\Shawn_PhDExperiments', ...
    'Navigate to the directory to save the zero crossing graphic');
    ZeroCrossFilename1 = ['PRE1' '_ZeroCrossing.tif'];
    ZeroCrossFilename2 = ['PRE1' 'ZeroCrossing.svg'];
    ZeroCrossFilename3 = ['PRE1' 'ZeroCrossing.pdf'];
    ZeroCrossFullFilename1 = fullfile( Directory,ZeroCrossFilename1);
    ZeroCrossFullFilename2 = fullfile( Directory,ZeroCrossFilename2);
    ZeroCrossFullFilename3 = fullfile( Directory,ZeroCrossFilename3);
    print(hh,'-dtiff',ZeroCrossFullFilename1,'-r600');
    print(hh,'-dsvg',ZeroCrossFullFilename2,'-r600');
    print(hh,'-dpdf',ZeroCrossFullFilename3,'-r600');
    clear hh
    clear output
    
    h = figure(10);
    plot(N_Reverse,WSEsave_shift(1,:),'--','color',linecolors1(6,:),'LineWidth',1.0),xlabel('Channel Station (mm)'),ylabel('Bed Elevation (mm)');
    hold on
    plot(N_Reverse,etasave_shift(1,:),'color',linecolors1(6,:),'LineWidth',1.0);
    hold on
    plot(N_Reverse,WSEsave_shift(end,:),'color',linecolors(8,:),'LineWidth',2.5);
    hold on
    plot(N_Reverse,etasave_shift(end,:),'color',linecolors2(4,:),'LineWidth',2.5);
    hold on
    legend({'Initial Bed Surface','Initial Water Surface','Final Bed Surface','Final Water Surface','Evolving Bed Surface'},'Location','northwest','FontSize',8)
    % Format stuff
    xlim([0 16000]);
    ylim([200 500]);
    a = get(gca,'XTickLabel');
    b = get(gca,'YTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',14)
    set(gca,'YTickLabel',b,'FontName','Times','fontsize',14)
    % Set-up a filename 
    % Specify the experiment and filename
    prompt = {'Exp','Filename'};
    dlg_title = 'Input';
    num_lines = 1;
    def = {'PRE1_','Simulated_Profile'};
    answers = inputdlg(prompt,dlg_title,num_lines,def);
    exp = answers{1};
    filename = answers{2};
    % Specify the filename
    ChannelResFilename = [exp filename '.png'];
    % Export a file and specify dimensions, etc.
    x0 = 10;
    y0 = 10;
    width = 20;
    height = 10;
    set(gcf, 'Units','centimeters','Position',[x0 y0 width height])
    set(gcf,'PaperType','B5','PaperPositionMode','auto') 
    print(h,'-dpng',ChannelResFilename,'-r600');
    clear h
    
    figure(1)        
    subplot(2,1,1);
    plot(N,etasave(1,:),'-k','Linewidth',2.5),xlabel('channel station (meters)'),ylabel('bed elevation (meters)'),title('Time Evolution of Bed Elevation');
    hold on
    plot(N,WSEsave(1,:),'-b','Linewidth',1.5);
    hold on
    plot(N,etasave(end,:),'-r','Linewidth',2.5);
    hold on
    plot(N,WSEsave(end,:),'--c','Linewidth',2.5);
    hold on
%     plot(N,etasave(2:end-1,:),'-y','Linewidth',0.5);
    legend('Initial Bed Surface','Initial Water Surface','Final Bed Surface','Final Water Surface','Evolving Bed Surface')
    
    subplot(2,1,2);
    plot(N,Qbsave(1,:),'-k'),xlabel('channel station (meters)'),ylabel('bedload transport (sq. m / s)'),title('Time Evolution of Bedload Transport');
    hold on
    plot(N,Qbsave(5,:),'b');
    hold on
    plot(N,Qbsave(100,:),'-r');
    hold on
    plot(N,Qbsave(end,:),'-g','Linewidth',2.5);
    legend('Time 1','Time 5','Time 100','Last Time')
    
    figure(2)
    subplot(3,1,1);
    plot(N,D90save(1,:)*1000,'-k'),xlabel('channel station (meters)'),ylabel('geometric mean grain size (mm)'),title('Time Evolution of Geometric Mean Grain Size');
    hold on
    plot(N,D90save(5,:)*1000,'b');
    hold on
    plot(N,D90save(50,:)*1000,'-r');
    hold on
    plot(N,D90save(end,:)*1000,'-g','Linewidth',2.5);
    legend('Time 1','Time 50','Time 500','Last Time')
    
    subplot(3,1,2);
    plot(N,ysave(1,:),'-k'),xlabel('channel station (meters)'),ylabel('water depth (meters)'),title('Time Evolution of Water Depth');
    hold on
    plot(N,ysave(5,:),'b');
    hold on
    plot(N,ysave(50,:),'-r');
    hold on
    plot(N,ysave(end,:),'-g','Linewidth',2.5);
    legend('Time 1','Time 50','Time 500','Last Time')    

    subplot(3,1,3);
    plot(N,Vsave(1,:),'-k'),xlabel('time (hours)'),ylabel('average velocity (meters per second)'),title('Time Evolution of Velocity');
    hold on
    plot(N,Vsave(5,:),'b');
    hold on
    plot(N,Vsave(50,:),'-r');
    hold on
    plot(N,Vsave(end,:),'-g','Linewidth',2.5);
    legend('Time 1','Time 50','Time 500','Last Time')
    
%     subplot(1,1,1);
%     semilogx(Runtime,Qwseries(:,1),'-b','Linewidth',1.5),axis([0 10 0 0.10]),xlabel('time (hours)'),ylabel('flow rate (sq. meters per second)'),title('Experimental Unit Hydrograph - Run 1');
%     hold on
%     semilogx(Runtime(:),Qwavgseries(:),'--c','Linewidth',1.5)
%     legend('Unit Hydrograph','Equivalent Flow');
%     
%     subplot(2,1,2);
%     plot(N,Frsave(2,:),'-k'),xlabel('time (hours)'),ylabel('Geometric Grain Size (millimeters)'),title('Time Evolution of Surface Geometric Grain Size');
%     hold on
%     plot(N,Frsave(4,:),'-b');
%     hold on
%     plot(N,Frsave(6,:),'-r');
%     hold on
%     plot(N,Frsave(10,:),'-g');
%     legend('Node 2','Node 24','Node 36','Node 60')
%     hold off
    
%     subplot(2,1,2);
%     plot(N,Dgsave1(1,:),'-k'),xlabel('channel station (meters)'),ylabel('D_s_g (millimeters)'),title('Time Evolution of D_s_g');
%     hold on
%     plot(N,Dgsave(10,:),'-b');
%     hold on
%     plot(N,Dgsave(50,:),'-r');
%     hold on
%     plot(N,Dgsave(100,:),'-g');
%     legend('Save 1','Save 10','Save 100','Save 1000')
%     hold off
%     
%     subplot(1,1,1);
%     plot(time,Fssave(:,1),'-k','Linewidth',2),axis([0 100 20 100]),xlabel('time (hours)'),ylabel('Sand Fraction of Surface and Load (percent)'),title('Time Evolution of Sand Fractions');
%     hold on
%     plot(time,Fssave(:,2),'-b','Linewidth',2);
%     hold on
%     plot(time,Fssave(:,4),'-r','Linewidth',2);
%     hold on
%     plot(time,Fssave(:,5),'-g','Linewidth',2);
%     hold on
%     plot(time,Qbsand(:,1),'-.k');
%     hold on
%     plot(time,Qbsand(:,2),'-.b');
%     hold on
%     plot(time,Qbsand(:,4),'-.r');
%     hold on
%     plot(time,Qbsand(:,5),'-.g');
%     hold on
%     plot(feedtime(1:11),feedsand(1:11),'ok','Linewidth',2);
%     legend('Surface Node 1','Surface Node 2','Surface Node 4','Surface Node 8','Load Node 1','Load Node 2','Load Node 4','Load Node 8','Feed Sand Content')
%     hold off

figure(3)    
subplot(2,1,1);
    plot(time,Dgsave(:,1).*1000,'-b'),xlabel('time (hours)'),ylabel('Geometric Grain Size (millimeters)'),title('Time Evolution of Surface Geometric Grain Size');
    legend('Node 8')
    
    subplot(2,1,2);
    plot(time,Sigmasave(:,1).*1000,'-r'),xlabel('time (hours)'),ylabel('Geometric Standard Deviation (millimeters)'),title('Time Evolution of Surface Standard Deviation');
    legend('Node 8')
    
% %     subplot(3,1,3);
% %     plot(time,Fssave(:,1),'-g','Linewidth',2),axis([0 1000 0 20]),xlabel('time (hours)'),ylabel('Sand Fraction of Surface (percent)'),title('Time Evolution of Surface Sand Fraction');
% %     legend('Node 8')
% % % %
%     subplot(2,1,1);
%     plot(Grain,GSDratio1(1,1:9),'-k','Linewidth',1.5),xlabel('Geometric Grain Size of Grain Classes (mm)'),ylabel('Ratio of Fk to fk'),title('Ratio of Grain Size Composition Between the Surface and Subsurface at t = 96 hours');
%     hold on
%     plot(Grain,GSDratio1(2,1:9),'-b','Linewidth',1.5);
%     hold on
%     plot(Grain,GSDratio1(3,1:9),'-r','Linewidth',1.5);
%     hold on
%     plot(Grain,GSDratio1(4,1:9),'-g','Linewidth',1.5);
%     hold on
%     plot(Grain,GSDratio1(5,1:9),'-c','Linewidth',1.5);
%     legend('Surface Node 1','Surface Node 2','Surface Node 4','Surface Node 6','Surface Node 8')
% 
%     subplot(2,1,2);
%     plot(Grain,GSDratio2(1,1:9),'-k','Linewidth',1.5),xlabel('Geometric Grain Size of Grain Classes (mm)'),ylabel('Ratio of Fk to Pbk'),title('Ratio of Grain Size Composition Between the Surface and the Load at t = 96 hours');
%     hold on
%     plot(Grain,GSDratio2(2,1:9),'-b','Linewidth',1.5);
%     hold on
%     plot(Grain,GSDratio2(3,1:9),'-r','Linewidth',1.5);
%     hold on
%     plot(Grain,GSDratio2(4,1:9),'-g','Linewidth',1.5);
%     hold on
%     plot(Grain,GSDratio2(5,1:9),'-c','Linewidth',1.5);
%     legend('Surface Node 1','Surface Node 2','Surface Node 4','Surface Node 6','Surface Node 8')
%
   
% subplot(1,1,1);
%     semilogy(time,fractrans(:,1),'-k','Linewidth',3),axis([0 100 10.5^-6 10^-4]),xlabel('time (hours)'),ylabel('Fractional Transport Rate (m^2/second)'),title('Time Evolution of Fractional Transport Rate at Surface Node 8');
%     hold on
%     semilogy(time,fractrans(:,2),'-b','Linewidth',3);
%     hold on
%     semilogy(time,fractrans(:,3),'-r','Linewidth',3);
%     hold on
%     semilogy(time,fractrans(:,4),'-g','Linewidth',3);
%     hold on
%     semilogy(time,fractrans(:,5),'-c','Linewidth',3);
%     hold on
%     semilogy(time,fractrans(:,6),'-m','Linewidth',3);
%     hold on
%     semilogy(time,fractrans(:,7),'-.k','Linewidth',1.5);
%     hold on
%     semilogy(time,fractrans(:,8),'-.b','Linewidth',1.5);
%     hold on
%     semilogy(time,fractrans(:,9),'-.r','Linewidth',1.5);
%     hold on
%     semilogy(time,fractrans(:,10),'-.g','Linewidth',1.5);
%     hold on
%     semilogy(time,fractrans(:,11),'-.c','Linewidth',1.5);
%     legend('0.2 mm','0.4 mm','0.7 mm','1.4 mm','2.8 mm','5.7 mm','11.3 mm','22.6 mm','45.3 mm','90.5 mm','181 mm')
%     hold off

%     subplot(4,2,1);
%     semilogx(GSG(:,1),GSG(:,3),'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',6),xlabel('grain size (millimeters)'),ylabel('cummulaitve probability distribution (percent)'),title('Grain Size Distribution');
%     hold on
%     semilogx(GSG(:,1),CummDistsave(j,:),'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',6),xlabel('grain size (millimeters)'),ylabel('cummulaitve probability distribution (percent)');
%     legend('Initial','Final')

%cmap = colormap(autumn(NoS));
%             cmap1 = colormap(winter(NoS));
            
            
%             subplot(2,1,1)
%             plot(N,WSEsave(index,:)-WSEsave1(1,:),'color',cmap1(index,:),'Linewidth',1.5);
%             M(index,:) = strcat('Location (NoS) = ',num2str(NoS(1,:)));
%             hold on
%             plot(N,n1save1(1,:),'-y','Linewidth',1.5);
%             hold on
%             plot(N,WSEsave1(1,:),'-r','Linewidth',1.5);
%             axis([0 15 1 1.40])
%             set(gca,'XTick',[0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30])
%             xlabel('channel station (m)'),ylabel('bed elevation and wse (m)');
%             legend('Initial WSE');
%             title(sprintf('PRFreeEvolve2CoarsePatches - Simulation run time = %0.2f (hours)',T/3600));
            %title(sprintf('percent complete = %0.2f',index/NoS));
            
%            
%             subplot(2,2,2)
%             plot(N,Frsave1(1,:),'-b','Linewidth',1.5);
% %             axis([0 15 1.0 1.40])
%             set(gca,'XTick',[0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30])
%             xlabel('channel station (m)'),ylabel('geometric mean D_i (mm)');
%             legend('Evolving D_g_m');            
%             title(sprintf('computation run time (seconds) = %0.2f',TotalTime));
            
%             hold on
            %subplot(2,1,1)
%             figure (3)
%             plot(N,etasave(index,:),'color',cmap(index,:),'Linewidth',1.5);
%             M(index,:) = strcat('Location (NoS) = ',num2str(NoS(1,:)));
%             hold on
%             plot(N,WSEsave(index,:),'color',cmap1(index,:),'Linewidth',1.5);
%             M(index,:) = strcat('Location (NoS) = ',num2str(NoS(1,:)));
%             set(gca,'XTick',[0 10 20 30 40 50 60 70 80 90 100])
%             xlabel('channel station (m)'),ylabel('elevation (m)');
%             legend('Evolving Bed','Evolving Water Surface');
%             title(sprintf('computation run time = %0.2f (seconds)',TotalTime));
            %title(sprintf('simulation run time (hours) = %0.2f',T/3600));
            
%            

%             hold on
%             subplot(2,2,4)
%             plot(N,Frsave(index,:),'color',cmap1(index,:),'Linewidth',1.5);
%             M(index,:) = strcat('Location (NoS) = ',num2str(NoS(1,:)));
%             hold on
% %             axis([0 15 1.0 1.40])
%             set(gca,'XTick',[0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30])
%             xlabel('channel station (m)'),ylabel('geometric mean D_i (mm)');
%             legend('Evolving D_g_m');
            
            %drawnow

end

