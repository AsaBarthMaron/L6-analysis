function quantifyGain(indices,toneResponses,startingBins)
    
    [masterDataExp,slaveDataExp] = laserDelay_tuningCurve_secondAnalysis(indices,toneResponses,startingBins,0)
    
    for k = indices
        for i = toneResponses{k}
        	gainQuant{k,i} = (mean(slaveDataExp)-mean(masterDataExp{}))/(mean(slaveDataExp)+mean(masterDataExp{}))
        end
    end
   