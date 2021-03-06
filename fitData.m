function[sdMaster,sdSlave,rsqMaster,rsqSlave] = fitData(masterData,slaveData)
    for j = 1:length(masterData)
        if j == 1
            layer = 'Layer 23';
        elseif j == 2
            layer = 'Layer 4';
        elseif j == 3
            layer = 'Layer 5';
        elseif j == 4
            layer = 'Layer 6';
        end
        for k = 1:length(masterData{j})
            master = masterData{j}{k};
            slave = slaveData{j}{k};
            if isempty(master) && isempty(slave)
                continue
            end
            muNum = k;
            muName = [layer '_MU#_' num2str(muNum)];
            figure
            subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.05], [0.03 0.03], [0.05 0.02]);
                for i = size(master,1):-1:1
                    try
                        yTmp = master(i,:)./max(master(i,:));
                        yIntrp = yTmp;
                    %    yIntrp = interp1(1:8,yTmp,1:.1:8,'spline');
                    %    s = fitoptions('Lower',[0,0,0],...
                     %                  'Upper',[Inf,Inf,Inf]);
                     %   f = fittype('gauss1','options',s);      
                        [fitobject, gof] = fit((1:length(yIntrp))',yIntrp','gauss1');
                        sdMaster{j}(i,k) = fitobject.a1;
                        rsqMaster{j}(i,k) = gof.rsquare;      
                        subplot(9,2,i)
                        %plot(yTmp,'k*')
                        plot(yIntrp)
                        hold on
                        plot(fitobject,'k')
                    catch
                        sdMaster{j}(i,k) = nan;
                        rsqMaster{j}(i,k) = nan; 
                        [muName 'master' num2str((i*25)-25) ' ms']
                    end
                    try
                        yTmp = slave(i,:)./max(slave(i,:));
                        yIntrp = yTmp;
                        %yIntrp = interp1(1:8,yTmp,1:.1:8,'spline');
                        [fitobject, gof] = fit((1:length(yIntrp))',yIntrp','gauss1');
                        sdSlave{j}(i,k) = fitobject.a1;
                        rsqSlave{j}(i,k) = gof.rsquare;      
                        subplot(9,2,i)
                        %plot(yTmp,'k*')
                        plot(yIntrp,'g')
                        hold on
                        plot(fitobject,'m')
                        ylabel([num2str((i*25)-25) ' ms'])
                        set(gca,'XTickLabel',{'4','5.7','8','11.3','16','22.6','32','45.3'});
                        subplot(9,2,1)
                        title(muName,'Interpreter','none')
                    catch
                        sdSlave{j}(i,k) = nan;
                        rsqSlave{j}(i,k) = nan;
                        [muName 'slave' num2str((i*25)-25) ' ms']
                    end
                end
                set(gcf,'position',[0,100,1200,1920]);
                export_fig(['C:\Users\polley_lab\Documents\MATLAB\cortex_genericlaser+tuningmoddelay\Images\Norm_tuning_curves_C\' muName ])
                close
        end
    end
end


