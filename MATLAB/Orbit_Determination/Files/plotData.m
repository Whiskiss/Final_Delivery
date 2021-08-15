function plotData( dataAll, a, b )
%PLOTDATA Summary of this function goes here
%   Detailed explanation goes here
    if size(dataAll,2)==1
        for idx=2:size(dataAll,1) 
            figure;
            for jdx=1:size(dataAll{idx,1},1)
                subplot(size(dataAll{idx,1},1),1,jdx);
                plot(dataAll{idx,1}(jdx,a:min(b,size(dataAll{idx,1},2))));

                grid on;
                if jdx == 1
                    title(idx);
                    legend 1 2;
                end
            end
        end
    else if size(dataAll,2)>1
        figure;
        subplot(2,1,1);
        plot(dataAll{1,1}(:,7));
        title('El');
        hold on;
        plot(dataAll{1,2}(:,7));
        grid on;   
        subplot(2,1,2);
        plot(dataAll{1,1}(:,8));
        title('Az');
        hold on;
        plot(dataAll{1,2}(:,8));
        grid on;  

        for idx=2:size(dataAll,1) 
            figure;
            for jdx=1:size(dataAll{idx,1},1)
                subplot(size(dataAll{idx,1},1),1,jdx);
                plot(dataAll{idx,1}(jdx,a:min(b,size(dataAll{idx,1},2))));
                hold on;
                plot(dataAll{idx,2}(jdx,a:min(b,size(dataAll{idx,2},2))));
%                 if idx == 4
%                     temp = dataAll{1,1}(a:min(b,size(dataAll{idx,1},2)),jdx+8);
%                     if jdx > 2
%                         temp = mod(deg2rad(temp),2*pi);
%                     end
%                     plot(temp);              
%                 end
                grid on;
                if jdx == 1
                    title(idx);
                    legend 1 2;
                end
            end
        end
	end
end

