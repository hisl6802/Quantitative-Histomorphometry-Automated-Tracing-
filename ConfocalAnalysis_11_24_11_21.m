clear; clc; close all force

inputs = input('How many analyses would like to perform? ');clc;
files = cell(inputs,2);
for a = 1:inputs
    %prompt the user inputs times to allow them to run as many analyses as
    %possible in a single run of code
    [file path] = uigetfile('*.png');
    fileImg = fullfile(path,file);
    files{a,1} = fileImg;

    %open the BP image and select the rectangle for analysis, then save
    %coordinates for usage in later analysis
    BP = imread(files{a,1});
    imshow(BP);
    rect = getrect();
    files{a,3} = round(rect(1));
    files{a,4} = round(rect(2));
    width = round(rect(3));
    height = round(rect(4));
    files{a,5} = files{a,3} + width;
    files{a,6} = files{a,4} + height;
    %% getting the scale for each analyzed image
    %points selected. 
    [scaleFile path] = uigetfile('*.png');
    scaleFile = string(path) + string(scaleFile);
    scaleFile = char(scaleFile);

    scaleImg = imread(scaleFile);
    imshow(scaleImg);
    rectScale = getrect();
    %determine the x,y, width,height for scale bar rectangle.
    xmin = round(rectScale(1));
    ymin = round(rectScale(2));
    width = round(rectScale(3));
    height = round(rectScale(4));
    xmax = xmin + width;
    ymax = ymin + height;

    BW_scale = zeros(size(scaleImg,1),size(scaleImg,2));

    for i = 1:size(scaleImg,1)
        for j = 1:size(scaleImg,2)
            if scaleImg(i,j,1) > 250 && scaleImg(i,j,2) > 250 && scaleImg(i,j,3) > 250
               if j < xmax && j > xmin % Need to automate the boundaries of the image. 
                   if i < ymax && i > ymin
                       BW_scale(i,j) = 1;
                   end
               end
            end
        end
    end

    BW_scale = imbinarize(BW_scale);
    imshow(BW_scale);
    s_scale = regionprops(BW_scale,'all');


    %find the scale factor for each set of images (i.e., 43L_1_BP and Labels)
    scale = 200/s_scale(1).BoundingBox(1,3);
    files{a,7} = scale;

    file = fileImg(1:end-19);
    file = string(file) + 'Labels.png';
    file = char(file);
    files{a,2} = file;
    close all force

end 
%create a cell array to store the data
outputParam = cell(4,inputs+1);
outputParam{1,1} = 'Parameters';
outputParam{2,1} = 'BS';
outputParam{3,1} = 'Calcein Surface';
outputParam{4,1} = 'Alizarin Surface';
outputParam{5,1} = 'Cal_to_BS';
outputParam{6,1} = 'Alz_to_BS';
outputParam{7,1} = 'Scale';

for a = 1:size(files,1)
    BP_total = 0;
    LabTotal = 0;
    clc; close all force

    %% Reading in files for analysis 
    BP = imread(files{a,1});
    Labels = imread(files{a,2});

    %% Getting region of interest for image that will be analyzed
    xmin = files{a,3};
    ymin = files{a,4};
    xmax = files{a,5};
    ymax = files{a,6};


    %% Bone Perimeter Analysis
    BW_BP = zeros(size(BP,1),size(BP,2));
    for i = 1:size(BP,1)
        for j = 1:size(BP,2)
            if BP(i,j,1) > 253 && BP(i,j,2) > 253 && BP(i,j,3) > 253
               if j < xmax && j > xmin % Need to automate the boundaries of the image. 
                   if i < ymax-10 && i > ymin+10
                       BW_BP(i,j) = 1;
                   end
               end
            end
        end
    end


    BW_BP = imbinarize(BW_BP);
    %get all the properties for the regions of interest
    % I think I will need the pixel list
    s_BP = regionprops(BW_BP,'PixelList');
    
    for p = 1:length(s_BP)
        
        if size(s_BP(p).PixelList,1) < 2500
            continue;
        end
        BW_Cur = zeros(size(BP,1),size(BP,2));
        for j = 1:size(s_BP(p).PixelList)
            BW_Cur(s_BP(p).PixelList(j,2),s_BP(p).PixelList(j,1)) = 1; 
        end
        
        %Binzarize the current region of interest
        BW_Cur_orig= imbinarize(BW_Cur);
        BW_Cur = imbinarize(BW_Cur);
        
        %% get a single regions properties and erode the images thick line.
        %get a region props output for single region
        SE = strel("disk",8);
        BW_Cur = imerode(BW_Cur,SE);
        %BW_Cur = imbinarize(BW_Cur);
        s_cur = regionprops(BW_Cur,'PixelList');
        radius = 8;
        while length(s_cur) > 1
            radius = radius -1;
            SE = strel("disk",radius);
            BW_Cur = imerode(BW_Cur_orig,SE);
            s_cur = regionprops(BW_Cur,'PixelList');
        end
        xx = min(s_cur(1).PixelList(:,1)):1:max(s_cur(1).PixelList(:,1));
        value = csapi(s_cur(1).PixelList(:,1),s_cur(1).PixelList(:,2),xx);

        %% find the X values in the region and there medians or means
        for k=1:length(xx)
            count = 0;
            %find the range x range in the pixel list
            curVals = find(s_cur(1).PixelList(:,1)==xx(k));
            curVals = s_cur(1).PixelList(curVals,:);
            
            if curVals(size(curVals,1),2)-curVals(1,2) > (size(curVals,1)-1)
                initialIndex = 1;
                for j=2:size(curVals,1)
                    %look for dicontinuities in y.
                    if curVals(j,2)-curVals(j-1,2) > 1
                        %found discontinuity
                        count = count + 1;
                        endIndex = j-1;
                        curCheck = curVals(initialIndex:endIndex,:);
                        if rem(size(curCheck,1),2) == 1
                            if size(curCheck,1) == 1
                                newVal = curCheck;
                            else
                                newVal = median(curCheck);
                            end

                            if count == 1
                                curFind = find(xx == curCheck(1,1));
                                value(curFind) = newVal(1,2);
                            else
                                xx(size(xx,2)+1) = newVal(1,1);
                                value(size(value,2)+1) = newVal(1,2);
                            end
                        else
                            newVal(1,1) = curCheck(size(curCheck,1)/2,1);
                            newVal(1,2) = curCheck(size(curCheck,1)/2,2);
                            
                            if count == 1
                                curFind = find(xx == curCheck(1,1));
                                value(curFind) = newVal(1,2);
                            else
                                xx(size(xx,2)+1) = newVal(1,1);
                                value(size(value,2)+1) = newVal(1,2);
                            end
                        end
                        initialIndex = endIndex + 1;
                    elseif j == size(curVals,1)
                        %found discontinuity
                        count = count + 1;
                        endIndex = j;
                        curCheck = curVals(initialIndex:endIndex,:);
                        if rem(size(curCheck,1),2) == 1
                            newVal = median(curCheck);
                            if count == 1
                                curFind = find(xx == curCheck(1,1));
                                value(curFind) = newVal(1,2);
                            else
                                xx(size(xx,2)+1) = newVal(1,1);
                                value(size(value,2)+1) = newVal(1,2);
                            end
                        else
                            newVal(1,1) = curCheck(size(curCheck,1)/2,1);
                            newVal(1,2) = curCheck(size(curCheck,1)/2,2);
                            
                            if count == 1
                                curFind = find(xx == curCheck(1,1));
                                value(curFind) = newVal(1,2);
                            else
                                xx(size(xx,2)+1) = newVal(1,1);
                                value(size(value,2)+1) = newVal(1,2);
                            end
                        end
                        initialIndex = 1;
                    end 
                end
            end
        end
        
        %% flip the object and find the unique Y values and there medians or means
        uniqueYVals = unique(s_cur(1).PixelList(:,2));
        
        for i = 1:length(uniqueYVals)
            %find the locations of the unique Y values in the pixel list
            curYCheck = find(s_cur(1).PixelList(:,2)==uniqueYVals(i));
            
            %convert the found list to an array of the (x,y) coordinates
            %for that unique y-value.
            curYChecks = s_cur(1).PixelList(curYCheck,:);
            
            %find any values in value that match the current Y value (+/-
            %0.5)
            matchY = find(value == uniqueYVals(i));
            
            %check for a single line in the Y direction.
            if curYChecks(size(curYChecks,1),1)-curYChecks(1,1) > (size(curYChecks,1)-1)
                initialIndex = 1;
                for j=2:size(curYChecks,1)
                    %look for dicontinuities in x.
                    if curYChecks(j,1)-curYChecks(j-1,1) > 1
                        %found discontinuity
                        endIndex = j-1;
                        curCheck = curYChecks(initialIndex:endIndex,:);
                        
                        if rem(size(curCheck,1),2) == 1
                            if size(curCheck,1) == 1
                                newVal = curCheck;
                            else
                                newVal = median(curCheck);
                            end
                            xx(size(xx,2)+1) = newVal(1,1);
                            value(size(value,2)+1) = newVal(1,2);
                        else
                            newVal(1,1) = curCheck(size(curCheck,1)/2,1);
                            newVal(1,2) = curCheck(size(curCheck,1)/2,2);
                            xx(size(xx,2)+1) = newVal(1,1);
                            value(size(value,2)+1) = newVal(1,2);
                        end
                        initialIndex = endIndex + 1;
                    elseif j == size(curYChecks,1)
                        %found discontinuity
                        count = count + 1;
                        endIndex = j;
                        curCheck = curYChecks(initialIndex:endIndex,:);
                        if rem(size(curCheck,1),2) == 1
                            newVal = median(curCheck);
                            xx(size(xx,2)+1) = newVal(1,1);
                            value(size(value,2)+1) = newVal(1,2);
                        else
                            newVal(1,1) = curCheck(size(curCheck,1)/2,1);
                            newVal(1,2) = curCheck(size(curCheck,1)/2,2);
                            xx(size(xx,2)+1) = newVal(1,1);
                            value(size(value,2)+1) = newVal(1,2);
                        end
                        initialIndex = 1;
                    end 
                end
            else
                if rem(size(curYChecks,1),2) == 1
                    if size(curYChecks,1) == 1
                        newVal = curYChecks;
                    else
                        newVal = median(curYChecks);
                    end
                    xx(size(xx,2)+1) = newVal(1,1);
                    value(size(value,2)+1) = newVal(1,2);
                else
                    newVal(1,1) = curYChecks(size(curYChecks,1)/2,1);
                    newVal(1,2) = curYChecks(size(curYChecks,1)/2,2);
                    xx(size(xx,2)+1) = newVal(1,1);
                    value(size(value,2)+1) = newVal(1,2);
                end
            end
        end
        clear curCord 
        
        %% pair down your fit to a single set of points (no duplicates).
        %find the unique xx locations.
        uniqueXX = unique(xx);
        curCord = [];
        for i = 1:length(uniqueXX)
            %find the location of these unique values
            locsXX = find(xx == uniqueXX(i));
            if length(locsXX) == 1
                % add this to the output curCord immediately
                curCord(size(curCord,1)+1,1) = uniqueXX(i);
                curCord(size(curCord,1),2) = value(locsXX);
            elseif length(locsXX) > 1
                YY = value(locsXX);
                uniqueYY = unique(YY);
                if length(unique(YY)) == 1
                    curCord(size(curCord,1)+1,1) = uniqueXX(i);
                    curCord(size(curCord,1),2) = value(locsXX(1));
                else
                    for j =1:length(uniqueYY)
                        curCord(size(curCord,1)+1,1) = uniqueXX(i);
                        curCord(size(curCord,1),2) = uniqueYY(j);
                    end
                end
            else
                disp('Not all values will be in curCord');
            end
        end
        

        if p == 1
            imshow(BW_BP);
        end
        %% distance analysis
        BPdist = 0;
        startCord = curCord(1,:);
        %lineTrace = startCord;
        dists = pdist(curCord);
        dists = squareform(dists);
        for i = 1:size(dists,1)
            for j = 1:size(dists,2)
                if dists(i,j) > 15
                    dists(i,j) = 1e6;
                end
            end
        end
            
        %take the first coordinate and connect it to the next point.
        curDist = dists(:,1);
        nextStep = min(curDist(curDist>0));
        BPdist = BPdist+nextStep;
        %set the value in the distance matrix to 1e6 to remove from future
        %calculations.
        
        nextLoc = find(curDist == nextStep);
        if length(nextLoc) > 1
            %pick the second biggest value with the same distance.
            nextLoc = nextLoc(2);
        end
        
        %update the distances so they cannot connect the points which just
        %connected. 
        dists(nextLoc,1) = 1e6;
        dists(1,nextLoc) = 1e6;
        nextStep = curCord(nextLoc,:);
       
        LineShow = [startCord; nextStep];
        %hold on;plot(LineShow(:,1),LineShow(:,2),'r','LineWidth',2);
        curDist = dists(:,nextLoc);
        dists(nextLoc,:) = 1e6;
        dists(:,nextLoc) = 1e6;
        
        nextStep = min(curDist(curDist>0));
        BPdist = BPdist + nextStep;
        nextLoc1 = find(curDist == nextStep);
        if length(nextLoc1) > 1
            nextLoc1 = nextLoc1(2);
        end

        while nextStep < 1e6
            dists(nextLoc1,nextLoc) = 1e6;
            dists(nextLoc,nextLoc1) = 1e6;
            nextLoc = find(curDist == nextStep);
            if length(nextLoc) > 1
                nextLoc = nextLoc(2);
            end
            nextStep = curCord(nextLoc,:);
            LineShow = [LineShow; nextStep];
            curDist = dists(:,nextLoc);
            %update any nextLoc cooridinates to be all really large so no
            %reversing of the coordinates can occur. 
            dists(nextLoc,:) = 1e6;
            dists(:,nextLoc) = 1e6;
            nextStep = min(curDist(curDist>0));
            if nextStep < 1e6
                BPdist = BPdist+nextStep;
            end
            nextLoc1 = find(curDist == nextStep);
            if length(nextLoc1) > 1
                nextLoc1 = nextLoc1(2);
            end
        end
        secondCord = LineShow(2,:);
        LineShowOrig = LineShow;
        yy = smooth(LineShow(:,2));
        LineShow(:,2) = yy;
        if p == 1
            subplot(1,2,1);
            imshow(BW_BP);
            subplot(1,2,2);
            imshow(BW_BP);hold on;plot(LineShow(:,1),LineShow(:,2),'g','LineWidth',2.5);
        else
            subplot(1,2,2);
            hold on;plot(LineShow(:,1),LineShow(:,2),'g','LineWidth',2.5);hold on;
        end
        
        if size(LineShow(:,1),1)/size(curDist(:,1),1) < 0.9
            %update the list so as to remove points included in the initial
            %line analysis that didn't capture all the object.
            curLineLocs = [];
            for i = 1:size(LineShow,1)
                %get the first point from the line
                curLineCheck = LineShowOrig(i,:);
                for j = 1:size(curCord,1)
                    %get a coordinate to check points against.
                    curCordCheck = curCord(j,:);
                    check = curCordCheck == curLineCheck;
                    if check == [true true]
                        curLineLocs = [curLineLocs; j];
                    end
                end
            end
            curCordUpdated = [];
            
            for i =2:size(curCord,1)
                curCheck = find(curLineLocs == i);
                if isempty(curCheck)
                    curCordUpdated = [curCordUpdated; curCord(i,:)];
                end
            end
            curCordUpdated = [startCord; curCordUpdated];    
            LineShow = [];
            dists = pdist(curCordUpdated);
            dists = squareform(dists);
            for i = 1:size(dists,1)
                for j = 1:size(dists,2)
                    if dists(i,j) > 15
                        dists(i,j) = 1e6;
                    end
                end
            end
            
            %take the first coordinate and connect it to the next point.
            curDist = dists(:,1);
            nextStep = min(curDist(curDist>0));
            BPdist = BPdist+nextStep;
            nextLoc = find(curDist == nextStep);
            if length(nextLoc) > 1
                %pick the second biggest value with the same distance.
                nextLoc = nextLoc(1);
            end
            
            %update the distances so they cannot connect the points which just
            %connected. 
            dists(nextLoc,1) = 1e6;
            dists(1,nextLoc) = 1e6;
            nextStep = curCordUpdated(nextLoc,:);

            LineShow = [startCord; nextStep];
            curDist = dists(:,nextLoc);
            dists(nextLoc,:) = 1e6;
            dists(:,nextLoc) = 1e6;

            nextStep = min(curDist(curDist>0));
            BPdist = BPdist+nextStep;
            nextLoc1 = find(curDist == nextStep);
            if length(nextLoc1) > 1
                nextLoc1 = nextLoc1(1);
            end

            while nextStep < 1e6
                dists(nextLoc1,nextLoc) = 1e6;
                dists(nextLoc,nextLoc1) = 1e6;
                nextLoc = find(curDist == nextStep);
                if length(nextLoc) > 1
                    nextLoc = nextLoc(2);
                end
                nextStep = curCordUpdated(nextLoc,:);
                LineShow = [LineShow; nextStep];
                curDist = dists(:,nextLoc);
                %update any nextLoc cooridinates to be all really large so no
                %reversing of the coordinates can occur. 
                dists(nextLoc,:) = 1e6;
                dists(:,nextLoc) = 1e6;
                nextStep = min(curDist(curDist>0));
                if nextStep <1e6
                    BPdist=BPdist+nextStep;
                end
                nextLoc1 = find(curDist == nextStep);
                if length(nextLoc1) > 1
                    nextLoc1 = nextLoc1(2);
                end
            end
            yy = smooth(LineShow(:,2));
            LineShow(:,2) = yy;
            if p == 1
                hold on;plot(LineShow(:,1),LineShow(:,2),'r','LineWidth',2.5);
            else
                hold on;plot(LineShow(:,1),LineShow(:,2),'r','LineWidth',2.5);hold on;
            end
        end
        BP_total = BP_total + BPdist;
        outputParam{2,a+1} = BP_total;
    end
    close all force

    %% Calcein Labels Analysis 
    BW_Lab = zeros(size(Labels,1),size(Labels,2));
    for i = 1:size(Labels,1)
        for j = 1:size(Labels,2)
            %looking for the white labels
            if Labels(i,j,1) > 253 && Labels(i,j,2) > 253 && Labels(i,j,3) > 253
               if j < xmax && j > xmin % Need to automate the boundaries of the image. 
                   if i < ymax-10 && i > ymin+10
                       BW_Lab(i,j) = 1;
                   end
               end
            %looking for the blue labels 
            elseif Labels(i,j,3) > 240 && Labels(i,j,2) > 125 && Labels(i,j,1) < 120
                if j <xmax && j >xmin
                    if i < ymax-10 && i > ymin+10
                        BW_Lab(i,j) = 1;
                    end
                end
            end
        end
    end

    BW_Lab = imbinarize(BW_Lab);
    imshow(BW_Lab)
    %get all the properties for the regions of interest
    s_Lab = regionprops(BW_Lab,'all');
    TotalLabDist = 0;

    for p =1:length(s_Lab)
        if length(s_Lab(p).PixelList) < 25
            continue
        end
        %binarize image with a single object.
        clear BW_Cur
        x_start = s_Lab(p).BoundingBox(1,1)-.5;
        x_stop = x_start + 1 + s_Lab(p).BoundingBox(1,3);
        y_start = s_Lab(p).BoundingBox(1,2)+.5;
        y_stop = y_start + 1 + s_Lab(p).BoundingBox(1,4); 
        %Isolate the object for easy selection
        BW_Cur = BW_Lab(y_start:y_stop,x_start:x_stop); 
        close all force
        imshow(BW_Cur,'InitialMagnification',750)
        
        %number of drawings needed.
        numDraws = input('How many separate sets of points are needed? ');clc;
        if numDraws > 3
            numDrawsCheck = input('Are you wanting to draws more than 3 sets of points (Y or N)?','s');clc;
            if numDrawsCheck == 'N'
                numDraws = input('How many separate sets of points are needed? ');clc;
            end
        end
                
        for i = 1:numDraws
            %grab the appropriate number of dots for using geodesic distances then press enter to exit.
            [x y] = ginput;
            clear points
            points(:,1) = x;
            points(:,2) = y;
            if length(points) > 2
                %create an array of all solution points
                clear full_solution
                full_solution =[];
                %put the Cordinates in an array
                Cords = [];
                for i = 1:length(points)
                    %create a list of values which are white pixels. 
                    clear curCord
                    curCord = [];
                    curCord(1,1) = points(i,1);
                    curCord(1,2) = points(i,2);
                    for i=1:size(BW_Cur,1)
                        for j =1:size(BW_Cur,2)
                            if BW_Cur(i,j)
                                curCord(size(curCord,1)+1,1) = j;
                                curCord(size(curCord,1),2) = i;
                            end
                        end
                    end
                    curDist = pdist(curCord);
                    curDist = squareform(curDist);
                    curDist = curDist(:,1);
                    nextStep = min(curDist(curDist>0));
                    loc = find(curDist == nextStep);
                    if length(loc) > 1
                        loc = loc(1);
                    end

                    Cords(size(Cords,1)+1,:) = curCord(loc,:);   

                end

                for i = 2:length(Cords)
                    D1 = bwdistgeodesic(BW_Cur, Cords(i-1,1), Cords(i-1,2), 'quasi-euclidean');
                    D2 = bwdistgeodesic(BW_Cur, Cords(i,1), Cords(i,2), 'quasi-euclidean');
                    D = D1 + D2;
                    D = round(D * 8) / 8;

                    D(isnan(D)) = inf;
                    clear paths solution_path solution_points
                    paths = imregionalmin(D);

                    solution_path = bwmorph(paths, 'thin', inf);
                    solution_points = [];
                    for i = 1:size(solution_path,1)
                        for j = 1:size(solution_path,2)
                            if solution_path(i,j) == 1
                                solution_points(size(solution_points,1)+1,1) = i;
                                solution_points(size(solution_points,1),2) = j;
                            end
                        end
                    end
                    %put solution points into an array for analysis of distance
                    full_solution =[full_solution; solution_points];


                end
                %start = full_solution(1,:);
                %start = Cords(1,:);start = flip(start);
                start = points(1,:);start = flip(start);
                full_solution = unique(full_solution(:,1:2),'rows');
                full_solution = [start; full_solution];
                %% Distance analysis Bone Labels
                clear curCord
                Labdist = 0;

                %start = full_solution(1,:);
                locStart = ismember(full_solution,start,'rows');
                for i =1:length(locStart)
                    if locStart(i) == 1
                        startCord = full_solution(i,:);
                    end
                end
                %startCord = full_solution(1,:);
                %lineTrace = startCord;
                dists = pdist(full_solution);
                dists = squareform(dists);
                for i = 1:size(dists,1)
                    for j = 1:size(dists,2)
                        if dists(i,j) > 5
                            dists(i,j) = 1e6;
                        end
                    end
                end

                %take the first coordinate and connect it to the next point.
                curDist = dists(:,1);
                nextStep = min(curDist(curDist>0));
                if nextStep < 1e6
                    Labdist = Labdist+nextStep;
                end
                %set the value in the distance matrix to 1e6 to remove from future
                %calculations.

                nextLoc = find(curDist == nextStep);
                if length(nextLoc) > 1
                    %pick the second biggest value with the same distance.
                    nextLoc = nextLoc(2);
                end

                %update the distances so they cannot connect the points which just
                %connected. 
                dists(nextLoc,1) = 1e6;
                dists(1,nextLoc) = 1e6;
                nextStep = full_solution(nextLoc,:);

                LineShow = [startCord; nextStep];
                %hold on;plot(LineShow(:,1),LineShow(:,2),'r','LineWidth',2);
                curDist = dists(:,nextLoc);
                dists(nextLoc,:) = 1e6;
                dists(:,nextLoc) = 1e6;

                nextStep = min(curDist(curDist>0));
                Labdist = Labdist + nextStep;
                nextLoc1 = find(curDist == nextStep);
                if length(nextLoc1) > 1
                    nextLoc1 = nextLoc1(2);
                end

                while nextStep < 1e6
                    dists(nextLoc1,nextLoc) = 1e6;
                    dists(nextLoc,nextLoc1) = 1e6;
                    nextLoc = find(curDist == nextStep);
                    if length(nextLoc) > 1
                        nextLoc = nextLoc(2);
                    end
                    nextStep = full_solution(nextLoc,:);
                    LineShow = [LineShow; nextStep];
                    curDist = dists(:,nextLoc);
                    %update any nextLoc cooridinates to be all really large so no
                    %reversing of the coordinates can occur. 
                    dists(nextLoc,:) = 1e6;
                    dists(:,nextLoc) = 1e6;
                    nextStep = min(curDist(curDist>0));
                    if nextStep < 1e6
                        Labdist = Labdist+nextStep;
                    end
                    nextLoc1 = find(curDist == nextStep);
                    if length(nextLoc1) > 1
                        nextLoc1 = nextLoc1(2);
                    end
                end
                close all force
                imshow(BW_Cur);hold on;plot(LineShow(:,2),LineShow(:,1),'g','LineWidth',2);
                
                %disp(Labdist)
                TotalLabDist = TotalLabDist + Labdist;
                %disp(TotalLabDist)

            else
                %create a list of values which are white pixels. 
                clear curCord
                curCord = [];
                curCord(1,1) = points(1,1);
                curCord(1,2) = points(1,2);
                for i=1:size(BW_Cur,1)
                    for j =1:size(BW_Cur,2)
                        if BW_Cur(i,j)
                            curCord(size(curCord,1)+1,1) = j;
                            curCord(size(curCord,1),2) = i;
                        end
                    end
                end
                curDist = pdist(curCord);
                curDist = squareform(curDist);
                curDist = curDist(:,1);
                nextStep = min(curDist(curDist>0));
                loc = find(curDist == nextStep);
                startCord = curCord(loc,:);

                clear curCord
                curCord = [];
                curCord(1,1) = points(2,1);
                curCord(1,2) = points(2,2);
                for i=1:size(BW_Cur,1)
                    for j =1:size(BW_Cur,2)
                        if BW_Cur(i,j)
                            curCord(size(curCord,1)+1,1) = j;
                            curCord(size(curCord,1),2) = i;
                        end
                    end
                end
                curDist = pdist(curCord);
                curDist = squareform(curDist);
                curDist = curDist(:,1);
                nextStep = min(curDist(curDist>0));
                loc = find(curDist == nextStep);
                endCord = curCord(loc,:);


                D1 = bwdistgeodesic(BW_Cur, startCord(:,1), startCord(:,2), 'quasi-euclidean');
                D2 = bwdistgeodesic(BW_Cur, endCord(:,1), endCord(:,2), 'quasi-euclidean');
                D = D1 + D2;
                D = round(D * 8) / 8;

                D(isnan(D)) = inf;
                paths = imregionalmin(D);

                solution_path = bwmorph(paths, 'thin', inf);
                solution_points = [];
                for i =1:size(solution_path,1)
                    for j = 1:size(solution_path,2)
                        if solution_path(i,j) == 1
                            solution_points(size(solution_points,1)+1,1) = i;
                            solution_points(size(solution_points,1),2) = j;
                        end
                    end
                end

                clear full_solution
                full_solution = solution_points;

                full_solution = unique(full_solution(:,1:2),'rows');
                %% Distance analysis Bone Labels
                clear curCord
                Labdist = 0;
                start = startCord;
                start = full_solution(1,:);
                locStart = ismember(full_solution,start,'rows');
                for i =1:length(locStart)
                    if locStart(i) == 1
                        startCord = full_solution(i,:);
                    end
                end
                %startCord = full_solution(1,:);
                %lineTrace = startCord;
                dists = pdist(full_solution);
                dists = squareform(dists);
                for i = 1:size(dists,1)
                    for j = 1:size(dists,2)
                        if dists(i,j) > 5
                            dists(i,j) = 1e6;
                        end
                    end
                end

                %take the first coordinate and connect it to the next point.
                curDist = dists(:,1);
                nextStep = min(curDist(curDist>0));
                if nextStep <1e6
                    Labdist = Labdist+nextStep;
                end
                    %set the value in the distance matrix to 1e6 to remove from future
                %calculations.

                nextLoc = find(curDist == nextStep);
                if length(nextLoc) > 1
                    %pick the second biggest value with the same distance.
                    nextLoc = nextLoc(2);
                end

                %update the distances so they cannot connect the points which just
                %connected. 
                dists(nextLoc,1) = 1e6;
                dists(1,nextLoc) = 1e6;
                nextStep = full_solution(nextLoc,:);

                LineShow = [startCord; nextStep];
                %hold on;plot(LineShow(:,1),LineShow(:,2),'r','LineWidth',2);
                curDist = dists(:,nextLoc);
                dists(nextLoc,:) = 1e6;
                dists(:,nextLoc) = 1e6;

                nextStep = min(curDist(curDist>0));
                Labdist = Labdist + nextStep;
                nextLoc1 = find(curDist == nextStep);
                if length(nextLoc1) > 1
                    nextLoc1 = nextLoc1(2);
                end

                while nextStep < 1e6
                    dists(nextLoc1,nextLoc) = 1e6;
                    dists(nextLoc,nextLoc1) = 1e6;
                    nextLoc = find(curDist == nextStep);
                    if length(nextLoc) > 1
                        nextLoc = nextLoc(2);
                    end
                    nextStep = full_solution(nextLoc,:);
                    LineShow = [LineShow; nextStep];
                    curDist = dists(:,nextLoc);
                    %update any nextLoc cooridinates to be all really large so no
                    %reversing of the coordinates can occur. 
                    dists(nextLoc,:) = 1e6;
                    dists(:,nextLoc) = 1e6;
                    nextStep = min(curDist(curDist>0));
                    if nextStep < 1e6
                        Labdist = Labdist+nextStep;
                    end
                    nextLoc1 = find(curDist == nextStep);
                    if length(nextLoc1) > 1
                        nextLoc1 = nextLoc1(2);
                    end
                end
                close all force
                imshow(BW_Cur,'InitialMagnification',750);hold on;plot(LineShow(:,2),LineShow(:,1),'g','LineWidth',2);
                %disp(Labdist)
                TotalLabDist = TotalLabDist + Labdist;
                %disp(TotalLabDist)
            end
        end
    end
    outputParam{3,a+1} = TotalLabDist;
    
    %% Alizarin Labels Analysis 
    clear TotalLabDist Labdist
    LabTotal = 0;
    TotalLabDist = 0;
    BW_Lab = zeros(size(Labels,1),size(Labels,2));
    for i = 1:size(Labels,1)
        for j = 1:size(Labels,2)
            %looking for the yellow labels 
            if Labels(i,j,1) > 150 && Labels(i,j,2) > 150 && Labels(i,j,3) < 50 && Labels(i,j,1) < 175 && Labels(i,j,2) < 175
                if j <xmax && j >xmin
                    if i < ymax-10 && i > ymin+10
                        BW_Lab(i,j) = 1;
                    end
                end
            %looking for the blue labels 
            elseif Labels(i,j,3) > 240 && Labels(i,j,2) > 125 && Labels(i,j,1) < 120
                if j <xmax && j >xmin
                    if i < ymax-10 && i > ymin+10
                        BW_Lab(i,j) = 1;
                    end
                end
            end
        end
    end

    BW_Lab = imbinarize(BW_Lab);
    imshow(BW_Lab)
    %get all the properties for the regions of interest
    s_Lab = regionprops(BW_Lab,'all');
    TotalLabDist = 0;
    for p =1:length(s_Lab)
        if length(s_Lab(p).PixelList) < 25
            continue
        end

        %binarize image with a single object.
        clear BW_Cur
        x_start = s_Lab(p).BoundingBox(1,1)-.5;
        x_stop = x_start + 1 + s_Lab(p).BoundingBox(1,3);
        y_start = s_Lab(p).BoundingBox(1,2)+.5;
        y_stop = y_start + 1 + s_Lab(p).BoundingBox(1,4); 
        %Isolate the object for easy selection
        BW_Cur = BW_Lab(y_start:y_stop,x_start:x_stop); 
        close all force
        imshow(BW_Cur,'InitialMagnification',750)
        
        %number of drawings needed.
        numDraws = input('How many separate sets of points are needed? ');clc;
        if numDraws > 3
            numDrawsCheck = input('Are you wanting to draws more than 3 sets of points (Y or N)?','s');clc;
            if numDrawsCheck == 'N'
                numDraws = input('How many separate sets of points are needed? ');clc;
            end
        end
        for i = 1:numDraws
            %grab the appropriate number of dots for using geodesic distances then press enter to exit.
            [x y] = ginput;
            clear points
            points(:,1) = x;
            points(:,2) = y;
            if length(points) > 2
                %create an array of all solution points
                clear full_solution
                full_solution =[];
                %put the Cordinates in an array
                Cords = [];
                for i = 1:length(points)
                    %create a list of values which are white pixels. 
                    clear curCord
                    curCord = [];
                    curCord(1,1) = points(i,1);
                    curCord(1,2) = points(i,2);
                    for i=1:size(BW_Cur,1)
                        for j =1:size(BW_Cur,2)
                            if BW_Cur(i,j)
                                curCord(size(curCord,1)+1,1) = j;
                                curCord(size(curCord,1),2) = i;
                            end
                        end
                    end
                    curDist = pdist(curCord);
                    curDist = squareform(curDist);
                    curDist = curDist(:,1);
                    nextStep = min(curDist(curDist>0));
                    loc = find(curDist == nextStep);
                    if length(loc) > 1
                        loc = loc(1);
                    end

                    Cords(size(Cords,1)+1,:) = curCord(loc,:);   

                end

                for i = 2:length(Cords)
                    D1 = bwdistgeodesic(BW_Cur, Cords(i-1,1), Cords(i-1,2), 'quasi-euclidean');
                    D2 = bwdistgeodesic(BW_Cur, Cords(i,1), Cords(i,2), 'quasi-euclidean');
                    D = D1 + D2;
                    D = round(D * 8) / 8;

                    D(isnan(D)) = inf;
                    clear paths solution_path solution_points
                    paths = imregionalmin(D);

                    solution_path = bwmorph(paths, 'thin', inf);
                    solution_points = [];
                    for i = 1:size(solution_path,1)
                        for j = 1:size(solution_path,2)
                            if solution_path(i,j) == 1
                                solution_points(size(solution_points,1)+1,1) = i;
                                solution_points(size(solution_points,1),2) = j;
                            end
                        end
                    end
                    %put solution points into an array for analysis of distance
                    full_solution =[full_solution; solution_points];


                end
                %start = full_solution(1,:);
                %start = Cords(1,:);start = flip(start);
                start = points(1,:);start = flip(start);
                full_solution = unique(full_solution(:,1:2),'rows');
                full_solution = [start; full_solution];
                %% Distance analysis Bone Labels
                clear curCord
                Labdist = 0;

                %start = full_solution(1,:);
                locStart = ismember(full_solution,start,'rows');
                for i =1:length(locStart)
                    if locStart(i) == 1
                        startCord = full_solution(i,:);
                    end
                end
                %startCord = full_solution(1,:);
                %lineTrace = startCord;
                dists = pdist(full_solution);
                dists = squareform(dists);
                for i = 1:size(dists,1)
                    for j = 1:size(dists,2)
                        if dists(i,j) > 5
                            dists(i,j) = 1e6;
                        end
                    end
                end

                %take the first coordinate and connect it to the next point.
                curDist = dists(:,1);
                nextStep = min(curDist(curDist>0));
                if nextStep <1e6
                    Labdist = Labdist+nextStep;
                end
                %set the value in the distance matrix to 1e6 to remove from future
                %calculations.

                nextLoc = find(curDist == nextStep);
                if length(nextLoc) > 1
                    %pick the second biggest value with the same distance.
                    nextLoc = nextLoc(2);
                end

                %update the distances so they cannot connect the points which just
                %connected. 
                dists(nextLoc,1) = 1e6;
                dists(1,nextLoc) = 1e6;
                nextStep = full_solution(nextLoc,:);

                LineShow = [startCord; nextStep];
                %hold on;plot(LineShow(:,1),LineShow(:,2),'r','LineWidth',2);
                curDist = dists(:,nextLoc);
                dists(nextLoc,:) = 1e6;
                dists(:,nextLoc) = 1e6;

                nextStep = min(curDist(curDist>0));
                if nextStep <1e6
                    Labdist = Labdist + nextStep;
                end
                nextLoc1 = find(curDist == nextStep);
                if length(nextLoc1) > 1
                    nextLoc1 = nextLoc1(2);
                end

                while nextStep < 1e6
                    dists(nextLoc1,nextLoc) = 1e6;
                    dists(nextLoc,nextLoc1) = 1e6;
                    nextLoc = find(curDist == nextStep);
                    if length(nextLoc) > 1
                        nextLoc = nextLoc(2);
                    end
                    nextStep = full_solution(nextLoc,:);
                    LineShow = [LineShow; nextStep];
                    curDist = dists(:,nextLoc);
                    %update any nextLoc cooridinates to be all really large so no
                    %reversing of the coordinates can occur. 
                    dists(nextLoc,:) = 1e6;
                    dists(:,nextLoc) = 1e6;
                    nextStep = min(curDist(curDist>0));
                    if nextStep < 1e6
                        Labdist = Labdist+nextStep;
                    end
                    nextLoc1 = find(curDist == nextStep);
                    if length(nextLoc1) > 1
                        nextLoc1 = nextLoc1(2);
                    end
                end
                close all force
                imshow(BW_Cur);hold on;plot(LineShow(:,2),LineShow(:,1),'g','LineWidth',2);
                
                %disp(Labdist)
                TotalLabDist = TotalLabDist + Labdist;
                %disp(TotalLabDist)

            else
                %create a list of values which are white pixels. 
                clear curCord
                curCord = [];
                curCord(1,1) = points(1,1);
                curCord(1,2) = points(1,2);
                for i=1:size(BW_Cur,1)
                    for j =1:size(BW_Cur,2)
                        if BW_Cur(i,j)
                            curCord(size(curCord,1)+1,1) = j;
                            curCord(size(curCord,1),2) = i;
                        end
                    end
                end
                curDist = pdist(curCord);
                curDist = squareform(curDist);
                curDist = curDist(:,1);
                nextStep = min(curDist(curDist>0));
                loc = find(curDist == nextStep);
                startCord = curCord(loc,:);

                clear curCord
                curCord = [];
                curCord(1,1) = points(2,1);
                curCord(1,2) = points(2,2);
                for i=1:size(BW_Cur,1)
                    for j =1:size(BW_Cur,2)
                        if BW_Cur(i,j)
                            curCord(size(curCord,1)+1,1) = j;
                            curCord(size(curCord,1),2) = i;
                        end
                    end
                end
                curDist = pdist(curCord);
                curDist = squareform(curDist);
                curDist = curDist(:,1);
                nextStep = min(curDist(curDist>0));
                loc = find(curDist == nextStep);
                endCord = curCord(loc,:);


                D1 = bwdistgeodesic(BW_Cur, startCord(:,1), startCord(:,2), 'quasi-euclidean');
                D2 = bwdistgeodesic(BW_Cur, endCord(:,1), endCord(:,2), 'quasi-euclidean');
                D = D1 + D2;
                D = round(D * 8) / 8;

                D(isnan(D)) = inf;
                paths = imregionalmin(D);

                solution_path = bwmorph(paths, 'thin', inf);
                solution_points = [];
                for i =1:size(solution_path,1)
                    for j = 1:size(solution_path,2)
                        if solution_path(i,j) == 1
                            solution_points(size(solution_points,1)+1,1) = i;
                            solution_points(size(solution_points,1),2) = j;
                        end
                    end
                end

                clear full_solution
                full_solution = solution_points;

                full_solution = unique(full_solution(:,1:2),'rows');
                %% Distance analysis Bone Labels
                clear curCord
                Labdist = 0;
                start = startCord;
                start = full_solution(1,:);
                locStart = ismember(full_solution,start,'rows');
                for i =1:length(locStart)
                    if locStart(i) == 1
                        startCord = full_solution(i,:);
                    end
                end
                %startCord = full_solution(1,:);
                %lineTrace = startCord;
                dists = pdist(full_solution);
                dists = squareform(dists);
                for i = 1:size(dists,1)
                    for j = 1:size(dists,2)
                        if dists(i,j) > 5
                            dists(i,j) = 1e6;
                        end
                    end
                end

                %take the first coordinate and connect it to the next point.
                curDist = dists(:,1);
                nextStep = min(curDist(curDist>0));
                Labdist = Labdist+nextStep;
                %set the value in the distance matrix to 1e6 to remove from future
                %calculations.

                nextLoc = find(curDist == nextStep);
                if length(nextLoc) > 1
                    %pick the second biggest value with the same distance.
                    nextLoc = nextLoc(2);
                end

                %update the distances so they cannot connect the points which just
                %connected. 
                dists(nextLoc,1) = 1e6;
                dists(1,nextLoc) = 1e6;
                nextStep = full_solution(nextLoc,:);

                LineShow = [startCord; nextStep];
                %hold on;plot(LineShow(:,1),LineShow(:,2),'r','LineWidth',2);
                curDist = dists(:,nextLoc);
                dists(nextLoc,:) = 1e6;
                dists(:,nextLoc) = 1e6;

                nextStep = min(curDist(curDist>0));
                Labdist = Labdist + nextStep;
                nextLoc1 = find(curDist == nextStep);
                if length(nextLoc1) > 1
                    nextLoc1 = nextLoc1(2);
                end

                while nextStep < 1e6
                    dists(nextLoc1,nextLoc) = 1e6;
                    dists(nextLoc,nextLoc1) = 1e6;
                    nextLoc = find(curDist == nextStep);
                    if length(nextLoc) > 1
                        nextLoc = nextLoc(2);
                    end
                    nextStep = full_solution(nextLoc,:);
                    LineShow = [LineShow; nextStep];
                    curDist = dists(:,nextLoc);
                    %update any nextLoc cooridinates to be all really large so no
                    %reversing of the coordinates can occur. 
                    dists(nextLoc,:) = 1e6;
                    dists(:,nextLoc) = 1e6;
                    nextStep = min(curDist(curDist>0));
                    if nextStep < 1e6
                        Labdist = Labdist+nextStep;
                    end
                    nextLoc1 = find(curDist == nextStep);
                    if length(nextLoc1) > 1
                        nextLoc1 = nextLoc1(2);
                    end
                end
                close all force
                %imshow(BW_Cur,'InitialMagnification',750);hold on;plot(LineShow(:,2),LineShow(:,1),'g','LineWidth',2);
                %disp(Labdist)
                TotalLabDist = TotalLabDist + Labdist;
                %disp(TotalLabDist)
            end
        end
    end
    outputParam{4,a+1} = TotalLabDist;
    
    outputParam{1,a+1} = files{a,1}(length(files{a,1})-45:length(files{a,1})-20);
    outputParam{5,a+1} = outputParam{3,a+1}/outputParam{2,a+1};
    outputParam{6,a+1} = outputParam{4,a+1}/outputParam{2,a+1};
    outputParam{7,a+1} = scale;
    
end


try
    disp(files{1}(50:length(files{1})))
catch
    disp("Didn't display file name");
end

[filename, pathname] = uiputfile({'*.csv'});
fullname = fullfile(pathname,filename);
T = cell2table(outputParam);
writetable(T,fullname);
