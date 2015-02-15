function[masterData, slaveData, psthData,repPsthData,repPsthStdData] = laserDelay_tuningCurve_firstAnalysis(path,speed)

load(path);

if strcmp(speed,'fast')
    % Create a vector (trials) containing the indexes of all the "real" trials
    chanCounter = 1;
    for chan = 1:6:96;
        
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
                        fr(i,1)=length(find(and(fr_t>start_points(1),fr_t<=end_points(1))))*1000/bin_ms;
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
elseif strcmp(speed,'smoothPSTH')
    % Create a vector (trials) containing the indexes of all the "real" trials
    chanCounter = 1;
    for chan = 1:6:96;
        
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
                        fr(i,1)=length(find(and(fr_t>start_points(1),fr_t<=end_points(1))))*1000/bin_ms;
                        % Calculate spontaneous spike rate
                        fr_spntns(i,1)=length(find(and(fr_t>start_points(1)+200,fr_t<=end_points(1)+200)))*1000/bin_ms;
                        fr_t(find(fr_t <= 0)) = 1/1000;
                        fr_t(find(fr_t>3500)) = 3500;
                        tmpRaster = zeros(dacInt.SampleRate * (str2num(dacInt.TotalDuration)/10000),1);
                        tmpRaster(ceil(fr_t*(dacInt.SampleRate/10000))) = 1;
                        psth_t(:,i) = downsample(quickPSTH(tmpRaster,70),10);
                        clear  fr_t tmpRaster;
                    end
                end
                
                % Calculate mean firing rates for the F outer loop value, J inner
                % loop value.
                fra_values(j,f)=nanmean(fr);
                fra_SEM(j,f)=std(fr)/sqrt(length(fr));
                %fra_values(j,f)=sum(fr);
                fra_spntns(j,f)=nanmean(fr_spntns);
                fra_values_se(j,f)=nanmean(fr)/sqrt(length(fr));
                
                % Construct PSTHs.
                % Calculate and store the repetiton averaged PSTH for each
                % delay/frequency combo.
                PSTHs(j,f,:) = (mean(psth_t'));
                % Store PSTHs for all 15 repetitions for each delay/frequency combo
                repPSTHs(j,f) = {psth_t'};
                repStdPSTHs(j,f) = std(psth_t(:));
                clear fr fr_t raster_t psth_t
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
        psthData{chanCounter} = PSTHs;
        repPsthData{chanCounter} = repPSTHs;
        repPsthStdData{chanCounter} = repStdPSTHs;
        chanCounter = chanCounter +1;
    end
    
elseif strcmp(speed,'PSTH')
    % Create a vector (trials) containing the indexes of all the "real" trials
    chanCounter = 1;
    for chan = 1:6:96;
        
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
                        fr(i,1)=length(find(and(fr_t>start_points(1),fr_t<=end_points(1))))*1000/bin_ms;
                        % Calculate spontaneous spike rate
                        fr_spntns(i,1)=length(find(and(fr_t>start_points(1)+200,fr_t<=end_points(1)+200)))*1000/bin_ms;
                        fr_t(find(fr_t <= 0)) = 1/1000;
                        fr_t(find(fr_t>3500)) = 3500;
                        tmpRaster = zeros(dacInt.SampleRate * (str2num(dacInt.TotalDuration)/10000),1);
                        tmpRaster(ceil(fr_t*(dacInt.SampleRate/10000))) = 1;
                        psth_t(:,i) = histc(find(tmpRaster>0),1:50:35000)*200;
                        clear  fr_t tmpRaster;
                    end
                end
                
                % Calculate mean firing rates for the F outer loop value, J inner
                % loop value.
                fra_values(j,f)=nanmean(fr);
                fra_SEM(j,f)=std(fr)/sqrt(length(fr));
                %fra_values(j,f)=sum(fr);
                fra_spntns(j,f)=nanmean(fr_spntns);
                fra_values_se(j,f)=nanmean(fr)/sqrt(length(fr));
                
                % Construct PSTHs.
                % Calculate and store the repetiton averaged PSTH for each
                % delay/frequency combo.
                % Store PSTHs for all 15 repetitions for each delay/frequency combo
                PSTHs(j,f,:) = (mean(psth_t'));
                repPSTHs(j,f) = {psth_t'};
                repStdPSTHs(j,f) = std(psth_t(:));
                clear fr fr_t raster_t psth_t
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
        psthData{chanCounter} = PSTHs;
        repPsthData{chanCounter} = repPSTHs;
        repPsthStdData{chanCounter} = repStdPSTHs;
        chanCounter = chanCounter +1;
    end
    
end
%% slave
chanCounter = 1;
for chan = 1:6:96;
    
    
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
                    fr_slave(i,1)=length(find(and(fr_t_slave>start_points(2),fr_t_slave<=end_points(2))))*1000/bin_ms;
                    fr_spntns_slave(i,1)=length(find(and(fr_t_slave>start_points(2)+200,fr_t_slave<=end_points(2)+200)))*1000/bin_ms;
                    %                 tmpRaster_slave = zeros(dacInt.SampleRate * (str2num(dacInt.TotalDuration)/1000),1);
                    %                 tmpRaster_slave(round(fr_t_slave*(dacInt.SampleRate/1000))) = 1;
                    %                 raster_t_slave(:,i) = tmpRaster_slave;
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


