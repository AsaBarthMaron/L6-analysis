% clear all
% load('Z:\File Transfer\jr\recordings\122414\pvcrexai32_8-1-3-FRA')
% 
% clc

% Create a vector (trials) containing the indexes of all the "real" trials
for i=1:size(SCL,1)
    real_trials(i)=~isempty(SCL(i).innerIndex);
end
trials=find(real_trials==1);

electrode=9;
unit=1;

chan=(electrode-1)*6+unit;

% Create a vector containing a list of all outer loop stimulus elements
for i=trials
    outer_stim_seq_vector(i)=SCL(i).outerIndex;
end

% Define start and end points
start_points=[stimChans{1, 1}.Gate.Delay];
end_points=[stimChans{1, 1}.Gate.Width];


start_points=start_points+8;
end_points=end_points-10;
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
            end
        end
        
        % Calculate mean firing rates for the F outer loop value, J inner
        % loop value
        fra_values(j,f)=nanmean(fr);
        %fra_values(j,f)=sum(fr);
        fra_spntns(j,f)=nanmean(fr_spntns);
        fra_values_se(j,f)=nanmean(fr)/sqrt(length(fr));
        
        clear fr fr_t
    end
    
end

% Plot it
figure;imagesc(1:length(outerSeq.master.values),innerSeq.master.values,(fra_values))
set(gca,'XTick',[1 11 21 31])
set(gca,'XTickLabel',{outerSeq.master.values([1 11 21 31])/1000})
% figure;imagesc(outerSeq.master.values,innerSeq.master.values,(fra_values))
% set(gca,'xscale','log')
set(gca,'YDir','normal')
colormap(hot)



%% slave

for i=trials
    outer_stim_seq_vector_slave(i)=SCL(i).outerSlave(1);
end

for f=1:length(unique( outer_stim_seq_vector_slave))
    
    curr_outer_slave=find( outer_stim_seq_vector_slave==f);
    if isempty(curr_outer_slave)
        curr_outer_slave=1:length(outer_stim_seq_vector_slave);
    end
    
    
    for i=1:length(curr_outer_slave)
        stim_seq_vector_slave(i)=SCL(curr_outer_slave(i)).innerSlave(1);
    end
    
    [value_slave,stim_seq_vector_sorted_slave] = sort( stim_seq_vector_slave);
    
    values_inner_slave=unique(value_slave);
    for j=1:length(values_inner_slave)
        curr_value_inner_slave=find(value_slave==values_inner_slave(j));
        
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
            end
        end
        
        fra_values_slave(j,f)=nanmean(fr_slave);
        fra_spntns_slave(j,f)=nanmean(fr_spntns_slave);
        fra_values_se_slave(j,f)=nanmean(fr_slave)/sqrt(length(fr_slave));
        
        clear fr fr_t
    end
    
end

figure;imagesc(1:length(outerSeq.slave(1).values),innerSeq.slave(1).values,(fra_values_slave))
set(gca,'XTick',[1 11 21 31])
set(gca,'XTickLabel',{outerSeq.slave(1).values([1 11 21 31])/1000})
% figure;imagesc(outerSeq.master.values,innerSeq.master.values,(fra_values))
% set(gca,'xscale','log')
set(gca,'YDir','normal')
colormap(hot)
