for i = 1:4
    for j=1:size(sdMaster{i},1)
        for k = 1:size(sdMaster{i},2)
            sdMaster{i}(j,k) = 1/((sdMaster{i}(j,k)^2)*2*pi);
            sdSlave{i}(j,k) = 1/((sdSlave{i}(j,k)^2)*2*pi);
            if (rsqMaster{i}(j,k)<0.8) | (rsqSlave{i}(j,k)<0.8)
                sdMaster{i}(j,k) = nan;
                sdSlave{i}(j,k) = nan;
            end
        end
    end
            
    sdMaster{i}(sdMaster{i} == 0) = nan;
    sdSlave{i}(sdSlave{i} == 0) = nan;
    gSelectivity(:,i) = nanmean((sdMaster{i}-sdSlave{i})./sdSlave{i},2);
end
    
figure
plot(0:25:400,fliplr(gSelectivity(:,1)')*100)
hold on
plot(0:25:400,fliplr(gSelectivity(:,2)')*100)
plot(0:25:400,fliplr(gSelectivity(:,3)')*100)
plot(0:25:400,fliplr(gSelectivity(:,4)')*100)