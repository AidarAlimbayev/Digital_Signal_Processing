function ALPR()
clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
imtool close all;	% Close all figure windows created by imtool.
load trained2

assignin('base','net',net);

%workspace;	% Make sure the workspace panel is showing.

% Change the current folder to the folder of this m-file.
% (The "cd" line of code below is from Brett Shoelson of The Mathworks.)
if(~isdeployed)
	cd(fileparts(which(mfilename))); % From Brett
end

	% Continue with the demo.  Do some initialization stuff.
	close all;
	fontSize = 16;
	figure;
	% Maximize the figure. 
	set(gcf, 'Position', get(0, 'ScreenSize')); 

	% Change the current folder to the folder of this m-file.
	% (The line of code below is from Brett Shoelson of The Mathworks.)
	if(~isdeployed)
		cd(fileparts(which(mfilename)));
	end

	% Open an image
		% They want to pick their own.
		% Change default directory to the one containing the standard demo images for the MATLAB Image Processing Toolbox. 
		originalFolder = pwd; 
		folder = 'C:\Program Files\MATLAB\R2010a\toolbox\images\imdemos'; 
		if ~exist(folder, 'dir') 
			folder = pwd; 
		end 
		cd(folder); 
		% Browse for the image file. 
		[baseFileName, folder] = uigetfile('*.*', 'Specify an image file'); 
		fullImageFileName = fullfile(folder, baseFileName); 
		% Set current folder back to the original one. 
		cd(originalFolder);
		selectedImage = 'My own image'; % Need for the if threshold selection statement later.



	% Check to see that the image exists.  (Mainly to check on the demo images.)
	if ~exist(fullImageFileName, 'file')
		message = sprintf('This file does not exist:\n%s', fullImageFileName);
		uiwait(msgbox(message));
		return;
	end

	% Read in image into an array.
	[rgbImage storedColorMap] = imread(fullImageFileName); 
	[rows columns numberOfColorBands] = size(rgbImage); 
	% If it's monochrome (indexed), convert it to color. 
	% Check to see if it's an 8-bit image needed later for scaling).
	if strcmpi(class(rgbImage), 'uint8')
		% Flag for 256 gray levels.
		eightBit = true;
	else
		eightBit = false;
	end
	if numberOfColorBands == 1
		if isempty(storedColorMap)
			% Just a simple gray level image, not indexed with a stored color map.
			% Create a 3D true color image where we copy the monochrome image into all 3 (R, G, & B) color planes.
			rgbImage = cat(3, rgbImage, rgbImage, rgbImage);
		else
			% It's an indexed image.
			rgbImage = ind2rgb(rgbImage, storedColorMap);
			% ind2rgb() will convert it to double and normalize it to the range 0-1.
			% Convert back to uint8 in the range 0-255, if needed.
			if eightBit
				rgbImage = uint8(255 * rgbImage);
			end
		end
	end 
	
	% Display the original image.
	subplot(3, 4, 1);
	imshow(rgbImage);
	drawnow; % Make it display immediately. 
	if numberOfColorBands > 1 
		title('Original Color Image', 'FontSize', fontSize); 
	else 
		caption = sprintf('Original Indexed Image\n(converted to true color with its stored colormap)');
		title(caption, 'FontSize', fontSize);
	end

	% Convert RGB image to HSV
	hsvImage = rgb2hsv(rgbImage);
	% Extract out the H, S, and V images individually
	hImage = hsvImage(:,:,1);
	sImage = hsvImage(:,:,2);
	vImage = hsvImage(:,:,3);
	
	% Display them.
	subplot(3, 4, 2);
	imshow(hImage);
	title('Hue Image', 'FontSize', fontSize);
	subplot(3, 4, 3);
	imshow(sImage);
	title('Saturation Image', 'FontSize', fontSize);
	subplot(3, 4, 4);
	imshow(vImage);
	title('Value Image', 'FontSize', fontSize);

	% Compute and plot the histogram of the "hue" band.
	hHuePlot = subplot(3, 4, 6); 
	[hueCounts, hueBinValues] = imhist(hImage); 
	maxHueBinValue = find(hueCounts > 0, 1, 'last'); 
	maxCountHue = max(hueCounts); 
	bar(hueBinValues, hueCounts, 'r'); 
	grid on; 
	xlabel('Hue Value'); 
	ylabel('Pixel Count'); 
	title('Histogram of Hue Image', 'FontSize', fontSize);

	% Compute and plot the histogram of the "saturation" band.
	hSaturationPlot = subplot(3, 4, 7); 
	[saturationCounts, saturationBinValues] = imhist(sImage); 
	maxSaturationBinValue = find(saturationCounts > 0, 1, 'last'); 
	maxCountSaturation = max(saturationCounts); 
	bar(saturationBinValues, saturationCounts, 'g', 'BarWidth', 0.95); 
	grid on; 
	xlabel('Saturation Value'); 
	ylabel('Pixel Count'); 
	title('Histogram of Saturation Image', 'FontSize', fontSize);

	% Compute and plot the histogram of the "value" band.
	hValuePlot = subplot(3, 4, 8); 
	[valueCounts, valueBinValues] = imhist(vImage); 
	maxValueBinValue = find(valueCounts > 0, 1, 'last'); 
	maxCountValue = max(valueCounts); 
	bar(valueBinValues, valueCounts, 'b'); 
	grid on; 
	xlabel('Value Value'); 
	ylabel('Pixel Count'); 
	title('Histogram of Value Image', 'FontSize', fontSize);

	% Set all axes to be the same width and height.
	% This makes it easier to compare them.
	maxCount = max([maxCountHue,  maxCountSaturation, maxCountValue]); 
	axis([hHuePlot hSaturationPlot hValuePlot], [0 1 0 maxCount]); 

	% Plot all 3 histograms in one plot.
	subplot(3, 4, 5); 
	plot(hueBinValues, hueCounts, 'r', 'LineWidth', 2); 
	grid on; 
	xlabel('Values'); 
	ylabel('Pixel Count'); 
	hold on; 
	plot(saturationBinValues, saturationCounts, 'g', 'LineWidth', 2); 
	plot(valueBinValues, valueCounts, 'b', 'LineWidth', 2); 
	title('Histogram of All Bands', 'FontSize', fontSize); 
	maxGrayLevel = max([maxHueBinValue, maxSaturationBinValue, maxValueBinValue]); 
	% Make x-axis to just the max gray level on the bright end. 
	xlim([0 1]); 

	% Now select thresholds for the 3 color bands.


	% Assign the low and high thresholds for each color band.

		% Take a guess at the values that might work for the user's image.
		hueThresholdLow = 0;
		hueThresholdHigh = graythresh(hImage);
		saturationThresholdLow = 0.5*graythresh(sImage);
		saturationThresholdHigh = 1.0;
		valueThresholdLow = graythresh(vImage);
		valueThresholdHigh = 1.0;
	

	% Show the thresholds as vertical red bars on the histograms.
	PlaceThresholdBars(6, hueThresholdLow, hueThresholdHigh);
	PlaceThresholdBars(7, saturationThresholdLow, saturationThresholdHigh);
	PlaceThresholdBars(8, valueThresholdLow, valueThresholdHigh);

	% Now apply each color band's particular thresholds to the color band
	hueMask = (hImage >= hueThresholdLow) & (hImage <= hueThresholdHigh);
	saturationMask = (sImage >= saturationThresholdLow) & (sImage <= saturationThresholdHigh);
	valueMask = (vImage >= valueThresholdLow) & (vImage <= valueThresholdHigh);

	% Display the thresholded binary images.
	fontSize = 16;
	subplot(3, 4, 10);
	imshow(hueMask, []);
	title('=   Hue Mask', 'FontSize', fontSize);
	subplot(3, 4, 11);
	imshow(saturationMask, []);
	title('&   Saturation Mask', 'FontSize', fontSize);
	subplot(3, 4, 12);
	imshow(valueMask, []);
	title('&   Value Mask', 'FontSize', fontSize);
	% Combine the masks to find where all 3 are "true."
	% Then we will have the mask of only the red parts of the image.
	yellowObjectsMask = uint8(hueMask & saturationMask & valueMask);
	subplot(3, 4, 9);
	imshow(yellowObjectsMask, []);
	caption = sprintf('Mask of Only\nThe Yellow Objects');
	title(caption, 'FontSize', fontSize);

	% Tell user that we're going to filter out small objects.
	smallestAcceptableArea = 100; % Keep areas only if they're bigger than this.

	% Open up a new figure, since the existing one is full.
	figure;  
	% Maximize the figure. 
	set(gcf, 'Position', get(0, 'ScreenSize'));

	% Get rid of small objects.  Note: bwareaopen returns a logical.
	yellowObjectsMask = uint8(bwareaopen(yellowObjectsMask, smallestAcceptableArea));
	subplot(3, 3, 1);
	imshow(yellowObjectsMask, []);
	fontSize = 13;
	caption = sprintf('bwareaopen() removed objects\nsmaller than %d pixels', smallestAcceptableArea);
	title(caption, 'FontSize', fontSize);

	% Smooth the border using a morphological closing operation, imclose().
	structuringElement = strel('disk', 4);
	yellowObjectsMask = imclose(yellowObjectsMask, structuringElement);
	subplot(3, 3, 2);
	imshow(yellowObjectsMask, []);
	fontSize = 16;
	title('Border smoothed', 'FontSize', fontSize);

	% Fill in any holes in the regions, since they are most likely red also.
	yellowObjectsMask = uint8(imfill(yellowObjectsMask, 'holes'));
	subplot(3, 3, 3);
	imshow(yellowObjectsMask, []);
	title('Regions Filled', 'FontSize', fontSize);


	% You can only multiply integers if they are of the same type.
	% (yellowObjectsMask is a logical array.)
	% We need to convert the type of yellowObjectsMask to the same data type as hImage.
	yellowObjectsMask = cast(yellowObjectsMask, class(rgbImage)); 

	% Use the yellow object mask to mask out the yellow-only portions of the rgb image.
	maskedImageR = yellowObjectsMask .* rgbImage(:,:,1);
	maskedImageG = yellowObjectsMask .* rgbImage(:,:,2);
	maskedImageB = yellowObjectsMask .* rgbImage(:,:,3);
	% Show the masked off red image.
	subplot(3, 3, 4);
	imshow(maskedImageR);
	title('Masked Red Image', 'FontSize', fontSize);
	% Show the masked off saturation image.
	subplot(3, 3, 5);
	imshow(maskedImageG);
	title('Masked Green Image', 'FontSize', fontSize);
	% Show the masked off value image.
	subplot(3, 3, 6);
	imshow(maskedImageB);
	title('Masked Blue Image', 'FontSize', fontSize);
	% Concatenate the masked color bands to form the rgb image.
	maskedRGBImage = cat(3, maskedImageR, maskedImageG, maskedImageB);
%    imwrite(maskedRGBImage,'elp2.png')
	% Show the masked off, original image.
	subplot(3, 3, 8);
	imshow(maskedRGBImage);
	fontSize = 13;
	caption = sprintf('Masked Original Image\nShowing Only the Yellow Objects');
	title(caption, 'FontSize', fontSize);
	% Show the original image next to it.
	subplot(3, 3, 7);
	imshow(rgbImage);
	title('The Original Image (Again)', 'FontSize', fontSize);

	% Measure the mean HSV and area of all the detected blobs.
	[meanHSV, boundingBox,areas, numberOfBlobs] = MeasureBlobs(yellowObjectsMask, hImage, sImage, vImage);
	if numberOfBlobs > 0
		fprintf(1, '\n----------------------------------------------\n');
		fprintf(1, 'Blob #, Area in Pixels, Mean H, Mean S, Mean V\n');
		fprintf(1, '----------------------------------------------\n');
		for blobNumber = 1 : numberOfBlobs
			fprintf(1, '#%5d, %14d, %6.2f, %6.2f, %6.2f\n', blobNumber, areas(blobNumber), ...
				meanHSV(blobNumber, 1), meanHSV(blobNumber, 2), meanHSV(blobNumber, 3));
		end
	else
		% Alert user that no yellow blobs were found.
		message = sprintf('No yellow blobs were found in the image:\n%s', fullImageFileName);
		fprintf(1, '\n%s\n', message);
		uiwait(msgbox(message));
    end

%imwrite(maskedRGBImage,'haha.jpg');

%figure,imshow(newim);
im1=rgb2gray(maskedRGBImage);
%im1=medfilt2(im1,[3 3]); %Median filtering the image to remove noise%lossy
BW = edge(im1,'sobel'); %finding edges
[w,h]=size(BW);
CC=bwconncomp(BW)
s=regionprops(CC,'Perimeter','BoundingBox','Area');
perimfiltered=find([s.Perimeter]>0.1*(w+h))
areafiltered1=find([s.Area]<0.005*w*h)
areafiltered2=find([s.Area]>0.0008*w*h)
filtered=intersect(perimfiltered,areafiltered1);
filtered=intersect(filtered,areafiltered2)
F=ismember(labelmatrix(CC),filtered);


%labelmatrix(CC)
%figure,imshow(im);
%figure,imshow(im1);

[m,n]=size(filtered);
a=[];
for i=1:n
I=imcrop(im1, s(filtered(i)).BoundingBox);
%figure,imshow(I);%ff
bw = im2bw(I,0.5);
bw2 = edu_imgcrop(bw);
charvec = edu_imgresize(bw2);
figure,imshow(bw);
selected_net = evalin('base','net');

result = sim(selected_net,charvec);

a=[a vec2ind(result)];
end
a
[m,n]= size(a);
lnumber=[];
for i=1:n
if  a(i)==1
        if i==1
    lnumber=[lnumber 'K'];
    else lnumber=[lnumber '4'];
    end
    
end

if a(i)==2
    lnumber=[lnumber 'P'];
end

if a(i)==3
    lnumber=[lnumber 'T'];
end

if a(i)==4
    lnumber=[lnumber 'B'];
end

if a(i)==5
    lnumber=[lnumber 'U'];
end

if a(i)==6
    lnumber=[lnumber 'M'];
end

if a(i)==7
    lnumber=[lnumber '2'];
end

if a(i)==8
    lnumber=[lnumber '9'];
end

if a(i)==9
    lnumber=[lnumber '5'];
end

if a(i)==10
    lnumber=[lnumber '4'];
end

if a(i)==11
    lnumber=[lnumber 'J'];
end
end
lnumber
fid=fopen('output.txt','w+');
fprintf(fid,'%s\n',lnumber);
fclose(fid);
winopen('output.txt');
return


function [meanHSV,boundingBox, areas, numberOfBlobs] = MeasureBlobs(maskImage, hImage, sImage, vImage)
	[labeledImage numberOfBlobs] = bwlabel(maskImage, 8);     % Label each blob so we can make measurements of it
	if numberOfBlobs == 0
		% Didn't detect any yellow blobs in this image.
		meanHSV = [0 0 0];
		areas = 0;
		return;
	end
	% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
	blobMeasurementsHue = regionprops(labeledImage, hImage, 'area', 'MeanIntensity','BoundingBox');   
	blobMeasurementsSat = regionprops(labeledImage, sImage, 'area', 'MeanIntensity','BoundingBox');   
	blobMeasurementsValue = regionprops(labeledImage, vImage, 'area', 'MeanIntensity','BoundingBox');   
	
	meanHSV = zeros(numberOfBlobs, 3);  % One row for each blob.  One column for each color.
	meanHSV(:,1) = [blobMeasurementsHue.MeanIntensity]';
	meanHSV(:,2) = [blobMeasurementsSat.MeanIntensity]';
	meanHSV(:,3) = [blobMeasurementsValue.MeanIntensity]';
	
  boundingBox= blobMeasurementsHue.BoundingBox
    % BoundingBox(:,1)=[blobMeasurementsHue.BoundingBox]';
   % BoundingBox(:,1)=[blobMeasurementsSat.BoundingBox]';
   % BoundingBox(:,1)=[blobMeasurementsValue.BoundingBox]';
    
	% Now assign the areas.
	areas = zeros(numberOfBlobs, 3);  % One row for each blob.  One column for each color.
	areas(:,1) = [blobMeasurementsHue.Area]';
	areas(:,2) = [blobMeasurementsSat.Area]';
	areas(:,3) = [blobMeasurementsValue.Area]';

	return; % from MeasureBlobs()
    
    function PlaceThresholdBars(plotNumber, lowThresh, highThresh)
	% Show the thresholds as vertical red bars on the histograms.
	subplot(3, 4, plotNumber); 
	hold on;
	maxYValue = ylim;
	maxXValue = xlim;
	hStemLines = stem([lowThresh highThresh], [maxYValue(2) maxYValue(2)], 'r');
	children = get(hStemLines, 'children');
	set(children(2),'visible', 'off');
	% Place a text label on the bar chart showing the threshold.
	fontSizeThresh = 14;
	annotationTextL = sprintf('%d', lowThresh);
	annotationTextH = sprintf('%d', highThresh);
	% For text(), the x and y need to be of the data class "double" so let's cast both to double.
	text(double(lowThresh + 5), double(0.85 * maxYValue(2)), annotationTextL, 'FontSize', fontSizeThresh, 'Color', [0 .5 0], 'FontWeight', 'Bold');
	text(double(highThresh + 5), double(0.85 * maxYValue(2)), annotationTextH, 'FontSize', fontSizeThresh, 'Color', [0 .5 0], 'FontWeight', 'Bold');
	
	% Show the range as arrows.
	% Can't get it to work, with either gca or gcf.
% 	annotation(gca, 'arrow', [lowThresh/maxXValue(2) highThresh/maxXValue(2)],[0.7 0.7]);

	return; % from PlaceThresholdBars()
    
