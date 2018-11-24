
figure
for j=1:1
    figure;
   scrsz = get(0,'ScreenSize'); figure('Position',[1 1 2000 2000]);
    

    for i=1:9
            
        subplot(3,3,i);
        plot(td(:,i));hold on;plot(asd(:,i),'r');   
        plot(td(:,i)-asd(:,i),'g')
                    

         
    end
    % export_fig (savefile,'-append', '-transparent')
end
%========================================================
