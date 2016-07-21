
allData = rand(10,2);
for i = 1:5
    dataToAppend = ones(10,2);
     allData = vertcat(allData, dataToAppend);
end


allData=ones(1,10);
for i = 1:5
    A = rand(1,10);
    allData = vertcat(allData,A);
end
allData=allData(2:end,:);




    
