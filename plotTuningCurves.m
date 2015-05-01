function plotTuningCurves(master,slave,layer,muNum)
%layer = 'Layer 6';
muName = [layer '_MU#_' num2str(muNum)];
figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.05], [0.03 0.03], [0.05 0.02]);
    for i = size(master,1):-1:1
        subplot(9,2,i)
        plot(master(i,:)./max(master(i,:)))
        hold on
        plot(slave(i,:)./max(slave(i,:)),'g')
        ylabel([num2str((i*25)-25) ' ms'])
        set(gca,'XTickLabel',{'4','5.7','8','11.3','16','22.6','32','45.3'});
        subplot(9,2,1)
        title(muName,'Interpreter','none')
    end
    set(gcf,'position',[0,-1000,1200,1920]);
    export_fig(['C:\Users\polley_lab\Documents\MATLAB\cortex_genericlaser+tuningmoddelay\Images\Norm_tuning_curves_C\' muName ])
    close
end