function plotTuningCurves(master,slave)

subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.05], [0.05 0.06], [0.05 0.02]);
for i = 1:size(master,1)
    subplot(9,2,1)
    subplot(master(1,:))
    subplot(slave(1,:),'g')
end
end