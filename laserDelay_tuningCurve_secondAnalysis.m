function [masterDataExp,slaveDataExp] = laserDelay_tuningCurve_secondAnalysis(indices,toneResponses,startingBins,doSave)
fileList = descendingFileList('modLaserdelay');
for k = indices;
    path = ['C:\Users\asa\Documents\MATLAB\' fileList{k}];
    load(path);
    
    
    % Create a vector (trials) containing the indexes of all the "real" trials
    chanCounter = 1;
    for chan = 1:6:96;
        if find(toneResponses{k} == chanCounter) > 0
            isIncluded = 1;
        else
            isIncluded = 0;
        end
        for i=1:size(SCL,1)
            real_trials(i)=~isempty(SCL(i).innerIndex);
        end
        trials=find(real_trials==1);
        
        
        % Create a vector containing a list of all outer loop stimulus elements
        for i=trials
            outer_stim_seq_vector(i)=SCL(i).outerIndex;
        end
        
        % Define start and end points
        start_points=[stimChans{1, 1}.Gate.Delay+2,stimChans{1, 4}.Gate.Delay+2];
        end_points=[(stimChans{1, 1}.Gate.Width+stimChans{1, 1}.Gate.Delay)-2,(stimChans{1, 1}.Gate.Width+stimChans{1, 4}.Gate.Delay)-2];
        
        
        
        bin_ms=end_points(1)-start_points(1);
        
        % Loop through each unique outer loop element
        for f=1:length(unique( outer_stim_seq_vector))
            
            % Find indexes of all outer loop values equivalent to the current
            % unique outer loop value.
            curr_outer=find( outer_stim_seq_vector==f);
            % If none exist, then make curr_outer the size of
            % outer_stim_seq_vector?
            if isempty(curr_outer)
                curr_outer=1:length(outer_stim_seq_vector);
            end
            
            % For each outer loop value of the current type, find all corresponding
            % inner indeices.
            for i=1:length(curr_outer)
                stim_seq_vector(i)=SCL(curr_outer(i)).innerIndex;
            end
            
            % Sort inner indices (value), but also return vector of their original
            % index (stim_seq_vector_sorted
            [value,stim_seq_vector_sorted] = sort( stim_seq_vector);
            
            % Loop through unique inner values
            values_inner=unique(value);
            for j=1:length(values_inner)
                % Find indices of sorted inner values that are equal to the current
                % unique inner value
                curr_value_inner=find(value==values_inner(j));
                
                raster_t = [];
                
                rep=[curr_value_inner];
                % Loop for the number of times that each inner value appears.
                for i=1:length(rep)
                    % Find indeices of the spikes  on a specified channel during
                    % the whole of trial rep(i)
                    ch_choosen=find(SCL(curr_outer(stim_seq_vector_sorted(rep(i)))).ch==chan);
                    
                    if isempty(ch_choosen)
                        %display (['channel ' num2str(chan) ' was not defined in trial ' num2str(stim_seq_vector_sorted(i)) ])
                        fr(i,1)= 0;
                        fr_spntns(i,1)=0;
                    else
                        % Find the times of the spikes on a specified channel during
                        % the whole of trial rep(i)
                        fr_t=SCL(curr_outer(stim_seq_vector_sorted(rep(i)))).t(ch_choosen)- SCL(curr_outer(stim_seq_vector_sorted(rep(i)))).t(1);
                        % Calculate spike rate over the time window specieifed by
                        % start and end points
                        if isIncluded
                            newStart = ((startingBins(k,chanCounter)-601)*5)+1000;
                            fr(i,1)=length(find(and(fr_t>newStart,fr_t<=newStart+20)))*1000/bin_ms;
                        elseif ~isIncluded
                            fr(i,1)=length(find(and(fr_t>start_points(1),fr_t<=end_points(1))))*1000/bin_ms;
                        end
                        % Calculate spontaneous spike rate
                        fr_spntns(i,1)=length(find(and(fr_t>start_points(1)+200,fr_t<=end_points(1)+200)))*1000/bin_ms;
                        clear  fr_t;
                    end
                end
                
                % Calculate mean firing rates for the F outer loop value, J inner
                % loop value.
                fra_values(j,f)=nanmean(fr);
                fra_SEM(j,f)=std(fr)/sqrt(length(fr));
                %fra_values(j,f)=sum(fr);
                fra_spntns(j,f)=nanmean(fr_spntns);
                fra_values_se(j,f)=nanmean(fr)/sqrt(length(fr));
                
                clear fr fr_t
            end
            
        end
        
        % Plot it
        
        % subplot(8,2,chanCounter)
        % imagesc(1:length(outerSeq.master.values),-innerSeq.master.values,(fra_values))
        % set(gca,'XTick',[1 11 21 31])
        % set(gca,'XTickLabel',{outerSeq.master.values([1 11 21 31])/1000})
        % % figure;imagesc(outerSeq.master.values,innerSeq.master.values,(fra_values))
        % % set(gca,'xscale','log')
        % set(gca,'YDir','normal')
        % colormap(hot)
        
        % Cells containing data across all 16 channels
        masterData{chanCounter} = fra_values;
        masterSEMData{chanCounter} = fra_SEM;
        chanCounter = chanCounter +1;
    end
    
    %% slave
    chanCounter = 1;
    for chan = 1:6:96;
        if find(toneResponses{k} == chanCounter) > 0
            isIncluded = 1;
        else
            isIncluded = 0;
        end
        
        for i=trials
            outer_stim_seq_vector_slave(i)=SCL(i).outerSlave(1);
        end
        
        for f=1:length(unique( outer_stim_seq_vector_slave))
            
            curr_outer_slave=find( outer_stim_seq_vector_slave==f);
            if isempty(curr_outer_slave)
                curr_outer_slave=1:length(outer_stim_seq_vector_slave);
            end
            
            
            for i=1:length(curr_outer_slave)
                stim_seq_vector_slave(i)=SCL(curr_outer_slave(i)).innerIndex(1);
            end
            
            [value_slave,stim_seq_vector_sorted_slave] = sort( stim_seq_vector_slave);
            
            values_inner_slave=unique(value_slave);
            for j=1:length(values_inner_slave)
                curr_value_inner_slave=find(value_slave==values_inner_slave(j));
                
                raster_t = [];
                
                rep=[curr_value_inner_slave];
                for i=1:length(rep)
                    ch_choosen_slave=find(SCL(curr_outer_slave(stim_seq_vector_sorted_slave(rep(i)))).ch==chan);
                    
                    if isempty(ch_choosen_slave)
                        %  display (['channel ' num2str(chan) ' was not defined in trial ' num2str(stim_seq_vector_sorted_slave(i)) ])
                        fr_slave(i,1)= 0;
                        fr_spntns_slave(i,1)=0;
                    else
                        fr_t_slave=SCL(curr_outer_slave(stim_seq_vector_sorted_slave(rep(i)))).t(ch_choosen_slave)- SCL(curr_outer_slave(stim_seq_vector_sorted_slave(rep(i)))).t(1);
                        if isIncluded
                            newStart = ((startingBins(k,chanCounter)-601)*5)+3000;
                            fr_slave(i,1)=length(find(and(fr_t_slave>newStart,fr_t_slave<=newStart+20)))*1000/bin_ms;
                        elseif ~isIncluded
                            fr_slave(i,1)=length(find(and(fr_t_slave>start_points(2),fr_t_slave<=end_points(2))))*1000/bin_ms;
                        end
                        fr_spntns_slave(i,1)=length(find(and(fr_t_slave>start_points(2)+200,fr_t_slave<=end_points(2)+200)))*1000/bin_ms;
                        clear tmpRaster_slave fr_t_slave;
                    end
                end
                
                fra_values_slave(j,f)=nanmean(fr_slave);
                fra_SEM_slave(j,f)=std(fr_slave)/sqrt(length(fr_slave));
                fra_spntns_slave(j,f)=nanmean(fr_spntns_slave);
                fra_values_se_slave(j,f)=nanmean(fr_slave)/sqrt(length(fr_slave));
                
                %         rasters(j,f) = {mean(raster_t_slave')};
                clear fr fr_t
            end
            
        end
        
        % figure;imagesc(1:length(outerSeq.slave(1).values),innerSeq.slave(1).values,(fra_values_slave))
        % set(gca,'XTick',[1 11 21 31])
        % set(gca,'XTickLabel',{outerSeq.slave(1).values([1 11 21 31])/1000})
        % % figure;imagesc(outerSeq.master.values,innerSeq.master.values,(fra_values))
        % % set(gca,'xscale','log')
        % set(gca,'YDir','normal')
        % colormap(hot)
        
        slaveData{chanCounter} = fra_values_slave;
        slaveSEMData{chanCounter} = fra_SEM_slave;
        chanCounter = chanCounter +1;
    end
    masterDataExp(k,:) = masterData;
    slaveDataExp(k,:) = slaveData;
    
    
    
    
    for i = 1:length(masterData)
        diffData{i} = (masterData{i}-slaveData{i})./slaveData{i};
        %     diffData{i}(isnan(diffData{i})) = 0;
        %     diffData{i}(~isfinite(diffData{i})) = 0;
    end
    
    minRate = [min(min(masterData{1})),min(min(slaveData{1}))];
    maxRate = [max(max(masterData{1})),max(max(slaveData{1}))];
    minDiff = min(min(diffData{1}));
    maxDiff = max(max(diffData{1}));
    
    for i = 1:length(masterData)
        if minRate(1)>min(min(masterData{i}))
            minRate(1) = min(min(masterData{i}));
        end
        if minRate(2)>min(min(slaveData{i}))
            minRate(2) = min(min(slaveData{i}));
        end
        if minDiff(1)>min(min(diffData{i}))
            minDiff(1) = min(min(diffData{i}));
        end
        if maxRate(1)<max(max(masterData{i}))
            maxRate(1) = max(max(masterData{i}));
        end
        if maxRate(2)<max(max(slaveData{i}))
            maxRate(2) = max(max(slaveData{i}));
        end
        if maxDiff(1)<max(max(diffData{i}))
            maxDiff(1) = max(max(diffData{i}));
        end
    end
    
    figure
    j =1;
    rgbVals = 0:0.125:1;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.03], [0.05 0.02], [0.05 0.02]);
    
    % Plot difference functions alongside heatmaps for each channel
    for i = [7,13,4,12,5,15,2,16,1,14,3,11,6]
        clims = [min(min([masterData{i},slaveData{i}])),max(max([masterData{i},slaveData{i}]))];
        subplot(13,3,j)
        imagesc(length(outerSeq.master.values),innerSeq.master.values,slaveData{i},clims)
        if find(toneResponses{k} == i) > 0
            ylabel([num2str(i) ' - included'])
        else
            ylabel([num2str(i)])
        end
        set(gca,'box','off')
        colormap(hot)
        axis xy
        set(gca,'XTickLabel',{'4','5.7','8','11.3','16','22.6','32','45.3'});
        j = j+1;
        
        subplot(13,3,j)
        imagesc(length(outerSeq.master.values),innerSeq.master.values,masterData{i},clims)
        set(gca,'box','off')
        colormap(hot)
        axis xy
        set(gca,'XTickLabel',{'4','5.7','8','11.3','16','22.6','32','45.3'});
        j = j+1;
        
        subplot(13,3,j)
        h = 1;
        for line = 2:2:length(diffData{i})
            plot(diffData{i}(line,:),'color', [rgbVals(h),0,1-rgbVals(h)],'linewidth',2)
            hold on
            h = h+1;
        end
        plot(1:8,zeros(8,1),'k--')
        set(gca,'box','off')
        % axis([outerSeq.master.values(1)./1000 outerSeq.master.values(8)./1000 min(min(diffData{i})) max(max(diffData{i}))])
        axis([1 8 -1 8])
        set(gca,'XTickLabel',{'4','5.7','8','11.3','16','22.6','32','45.3'});
        j = j+1;
        
    end
    subplot(13,3,1)
    title('Tone')
    subplot(13,3,2)
    title('Tone + Laser')
    subplot(13,3,3)
    title('Difference')
    subplot(13,3,37)
    xlabel('Frequency (KHz)')
    subplot(13,3,38)
    xlabel('Frequency (KHz)')
    subplot(13,3,39)
    xlabel('Frequency (KHz)')
    
    if doSave
        if ~isdir(['C:\Users\asa\Documents\MATLAB\cortex_genericlaser+tuningmoddelay\Images\' fileList{k}(73:98) '\'])
            mkdir(['C:\Users\asa\Documents\MATLAB\cortex_genericlaser+tuningmoddelay\Images\' fileList{k}(73:98) '\']);
        end
        set(gcf,'position',[0,-1000,1200,1920]);
        export_fig(['C:\Users\asa\Documents\MATLAB\cortex_genericlaser+tuningmoddelay\Images\'  fileList{k}(73:end-4) '-2ndRndHeatmap.png'])
        close
    end
end

