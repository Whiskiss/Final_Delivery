elt = 1;
for iData = 1:size(GlobalRec,2)
    a = [];
    for jData = 1:size(GlobalRec,1)
    %a
        ilen = length([zeros([1, jData-1]) GlobalRec{jData,iData}{4,1}(elt,:)]);
        temp = [];
        for kData = 1:size(a,1)
           temp =  [ temp; [a(kData,:) zeros(1,ilen - size(a(kData,:),2))]];
        end
        a = temp;    
            
        a = [a; [zeros([1, jData-1]) GlobalRec{jData,iData}{4,1}(elt,:)]];
        %GlobalRec{jData,iData}{4,1}(1,:)
    
    %dummy = 
    %GlobalRec(size(GlobalRec,1)+1,iData)
    end
    figure; imagesc(a);
    diff = [];
    for jData=1:size(a,2)
        idx = find(a(:,jData)~=0);
        if idx>=1
            diff(jData) = max(a(idx,jData)) - min(a(idx,jData)) ;
        else
            diff(jData) = 0;
        end
    end
    figure; plot(diff);
    
end

%%
% close all;
% for jData = 20:min(30,size(GlobalRec,1))
%     figure;
%     for iData = 1:size(GlobalRec,2)
%         plot(GlobalRec{jData,iData}{4,1}(1,:));
%         hold on;
%         legend('data','dataMan','Location', 'Best');
%     end    
%     title(jData);
% end

%%
close all;
% value = {};
for iData = 1:size(GlobalRec,2)
    value{iData} = [];
    for jData = 1:size(GlobalRec,1)
        value{iData} = [value{iData}; GlobalRec{jData,iData}{4,1}(1,:)];

    end 
    figure;
    imagesc(value{iData}(1:end,:));
    figure; mesh(value{iData});
    xlabel('x: time');
    ylabel('y: time translation');
    zlabel('z: sma');
end


%%
close all;
%%
figure;
idx = 2;
for i = 1:size(value{idx},1)
    plot3([1:length(value{idx}(i,:))], ...
            value{idx}(i,:), ...
            i*ones(1,length(value{idx}(i,:))), 'r');
    hold on;
end
grid on;
xlabel('x: time');
ylabel('y: sma');
zlabel('z: time translation');


%%
figure;
idx = 2;
for i = 1:size(value{idx},2)
    plot3([1:length(value{idx}(:,i))], ...
            value{idx}(:,i), ...
            i*ones(1,length(value{idx}(:,i))));
    hold on;
end
grid on;
xlabel('x: time');
ylabel('y: sma');
zlabel('z: time translation');



%% Plot3 all elements on same figure
figure;
idx = 2;
elt = 2;
for i = 1:size(GlobalRec,1)
    plot3([1:length(GlobalRec{i,idx}{4,1}(elt,:))], ...
            GlobalRec{i,idx}{4,1}(elt,:), ...
            i*ones(1,length(GlobalRec{i,idx}{4,1}(elt,:))), 'g');
    hold on;
end
grid on;
xlabel('x: time');
ylabel('y: sma');
zlabel('z: time translation');

idx = 1;
for i = 1:size(GlobalRec,1)
    plot3([1:length(GlobalRec{i,idx}{4,1}(elt,:))], ...
            GlobalRec{i,idx}{4,1}(elt,:), ...
            i*ones(1,length(GlobalRec{i,idx}{4,1}(elt,:))), 'r');
    hold on;
end
grid on;
xlabel('x: time');
ylabel('y: sma');
zlabel('z: time translation');

%% plot only last elements of each iteration
idx = 1;
elt = 1;
noMan = [];
Man = [];
hMan = figure;
%hNoMan = figure;

for j = 1:size(GlobalRec{1,1}{4,1},2)
 noMan = [];
 Man = [];
%     for i = 1:size(GlobalRec,1)
%         noMan = [noMan; i, GlobalRec{i,1}{4,1}(elt,1), ...
%                         i, GlobalRec{i,1}{4,1}(elt,25), ...
%                         i, GlobalRec{i,1}{4,1}(elt,end)];
%         Man   = [Man  ; i, GlobalRec{i,2}{4,1}(elt,1), ...
%                         i, GlobalRec{i,2}{4,1}(elt,25), ...
%                         i, GlobalRec{i,2}{4,1}(elt,end)];
%     end
    for i = 1:size(GlobalRec,1)
        noMan = [noMan; i, GlobalRec{i,1}{4,1}(elt,j)];
        Man   = [Man  ; i, GlobalRec{i,2}{4,1}(elt,j)];
    end 
    figure(hMan);
    plot(Man(:,1), Man(:,2), 'r');
%     figure(hNoMan);
%     plot(noMan(:,1), noMan(:,2), 'g');
    pause(0.5);
end
figure;
plot(noMan(:,1), noMan(:,2), 'r');
hold on;
plot(noMan(:,3), noMan(:,4), 'g');
plot(noMan(:,5), noMan(:,6), 'b');
title('NoMan');

figure;
plot(Man(:,1), Man(:,2), 'r');
hold on;
plot(Man(:,3), Man(:,4), 'g');
plot(Man(:,5), Man(:,6), 'b');
title('Man');






