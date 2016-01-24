function [] = plot2D(x,Y,color,width)
subplot(2,1,1);hold on;grid on;
for k=1:size(Y,3)
    plot(x',Y(:,1,k),color,'lineWidth',width);
end
for i=3:4
    subplot(2,2,i);hold on;grid on;
    for k=1:size(Y,3)
        plot(x',Y(:,i,k),color,'lineWidth',width);
    end
end