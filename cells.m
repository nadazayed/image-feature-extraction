lbl='null';
clear area;
clear majorAxis;
clear solidity;
clear filled;
clear minorAxis;
clear circ;
clear convex;
clear ecc;
clear peri;
clear label;
%%clear previous variables.
img_folder = '/Users/nadaaz/NAZ/AAST/LISC Database/Main Dataset/Baso';
cd (img_folder);
file_name = dir(fullfile(img_folder,'*.bmp'));
figure;
totl_img = numel(file_name);
for n = 1:totl_img
    f= fullfile(img_folder,file_name(n).name);
    img = imread(f);
    
    original_image=img;
image = rgb2gray(original_image);
image = uint8(image);

image_adjust = imadjust(image); % Adjust intensity
image_histo = histeq(image); % histogram eq
% Addition to brighten most of details in image except the nuclei
image_add = imadd(image_adjust, image_histo); 
% Subtraction to highlight objects and borders including the nuclei
image_subtract = imsubtract(image_add, image_histo); 
% Addition to remove all other components but the nuclei
image = imadd(image_add, image_subtract); 

B = zeros(size(image)); %h3ml array of zeros b 7gm el image
A = padarray(image,[1 1]); % Pad the image, hy3ml frame black 7wliha

x=[1:3]';
y=[1:3]';       
for i = 1 : size(A, 1) - 2 %size of img 578  722 .. size(A,1) = 578
    for j = 1 : size(A, 2) - 2 %kol el cols 720
       % Vectorized method 
       window = reshape(A(i + x - 1, j + y - 1), [], 1);
       B(i,j) = min(window);
    end
end
B = uint8(B); % Convert to 0 - 255 intensity scale
image = B;
clear i; clear j;

image = imbinarize(image, graythresh(image));
image = imcomplement(image);

se = strel('disk', 7);
image = imopen(image, se); % Morphological opening

image_cc = bwconncomp(image); % Find connected components
% Calculate number of pixels of each connected component
numPixels = cellfun(@numel, image_cc.PixelIdxList); 
numPixels_max = max(numPixels);
for i = 1 : size(numPixels, 2)
    
    if contains(img_folder,'Baso')
        lbl=1;
        if numPixels(i) < 0.38*numPixels_max
            image(image_cc.PixelIdxList{i}) = 0;
        end
        
        
    elseif contains(img_folder,'eosi')
        lbl=2;
        if numPixels(i) < 0.985*numPixels_max
            image(image_cc.PixelIdxList{i}) = 0;
        end
        
    elseif contains(img_folder,'lymp')
        lbl=3;
        if numPixels(i) < 0.95*numPixels_max
            image(image_cc.PixelIdxList{i}) = 0;
        end
        
    elseif contains(img_folder,'mono')
        lbl=4;
        if numPixels(i) < 0.7*numPixels_max
            image(image_cc.PixelIdxList{i}) = 0;
        end
        
    elseif contains(img_folder,'neut')
        lbl=5;
        if numPixels(i) < 0.5*numPixels_max
            image(image_cc.PixelIdxList{i}) = 0;
        end
        
    end %el if
  
%%extract features of selected image.    
stats= regionprops(image,'Area','MajorAxisLength','Solidity','FilledArea','MinorAxisLength','Circularity','ConvexArea','Eccentricity','Perimeter');
area{n}=[stats.Area];
majorAxis{n}=[stats.MajorAxisLength];
solidity{n}=[stats.Solidity];
filled{n}=[stats.FilledArea];
minorAxis{n}=[stats.MinorAxisLength];
circ{n}=[stats.Circularity];
convex{n}=[stats.ConvexArea];
ecc{n}=[stats.Eccentricity];
peri{n}=[stats.Perimeter];
label{n}=lbl;    
end %el for

%result_image = labeloverlay(original_image, image);
imshow(image);

%%save images in current directory.
saveas(gcf,[num2str(n),'.jpg'])
end
%%
%%make a table of selected features.
clear data;

data=table(area',majorAxis',solidity',filled',minorAxis',circ',convex',ecc',peri',label');
%%
%%save minor values as initial value in the array.
minArea=            size(data.Var1{1});
minMajor=           size(data.Var2{1});
minSolidity=        size(data.Var3{1});
minFilled=          size(data.Var4{1});
minMinor=           size(data.Var5{1});
minCirc=            size(data.Var6{1});
minConvex=          size(data.Var7{1});
minEcc=             size(data.Var8{1});
minPeri=            size(data.Var9{1});
%%
%%check for minimum values compared to the first one.
for i = 2:totl_img
    if sum(size(data.Var1{i})) < sum(minArea)
        minArea=size(data.Var1{i});
    end
    
    if sum(size(data.Var2{i})) < sum(minMajor)
        minMajor=size(data.Var2{i}); 
    end
    
    if sum(size(data.Var3{i})) < sum(minSolidity)
        minSolidity=size(data.Var3{i}); 
    end
    
    if sum(size(data.Var4{i})) < sum(minFilled)
        minFilled=size(data.Var4{i}); 
    end
    
    if sum(size(data.Var5{i})) < sum(minMinor)
        minMinor=size(data.Var5{i}); 
    end
    
    if sum(size(data.Var6{i})) < sum(minCirc)
        minCirc=size(data.Var6{i}); 
    end
    
    if sum(size(data.Var7{i})) < sum(minConvex)
        minConvex=size(data.Var7{i}); 
    end
    
    if sum(size(data.Var8{i})) < sum(minEcc)
        minEcc=size(data.Var8{i}); 
    end
    
    if sum(size(data.Var9{i})) < sum(minPeri)
        minPeri=size(data.Var9{i}); 
    end
end
%%
%%remove un-equal cells in size 
for i = 1:totl_img
    if ~isequal(minArea,(size(data.Var1{i})))
        data([i],:) = [];
    end
end
%%
%%remove the rest.
for i = 1:totl_img
    if ~isequal(minArea,(size(data.Var1{i})))
        data([i],:) = [];
    end
end

for i = 1:totl_img
    if ~isequal(minMajor,(size(data.Var2{i})))
        data([i],:) = [];
    end
end

for i = 1:totl_img
    if ~isequal(minSolidity,(size(data.Var3{i}))) 
        data([i],:) = [];
    end
end

for i = 1:totl_img
    if ~isequal(minFilled,(size(data.Var4{i}))) 
        data([i],:) = [];
    end
end

for i = 1:totl_img
    if ~isequal(minMinor,(size(data.Var5{i}))) 
        data([i],:) = [];
    end
end

for i = 1:totl_img
    if ~isequal(minCirc,(size(data.Var6{i}))) 
        data([i],:) = [];
    end
end

for i = 1:totl_img
    if ~isequal(minConvex,(size(data.Var7{i}))) 
        data([i],:) = [];
    end
end

for i = 1:totl_img
    if ~isequal(minEcc,(size(data.Var8{i}))) 
        data([i],:) = [];
    end
end

for i = 1:totl_img
    if ~isequal(minPeri,(size(data.Var9{i}))) 
        data([i],:) = [];
    end
end
%%
for i = 1:totl_img
    if ~isequal(minPeri,(size(data.Var9{i}))) 
        data([i],:) = [];
    end
end
%%
%%convert cells to matrices in the table.
data.Var1 = cell2mat (data.Var1);
data.Var2 = cell2mat (data.Var2);
data.Var3 = cell2mat (data.Var3);
data.Var4 = cell2mat (data.Var4);
data.Var5 = cell2mat (data.Var5);
data.Var6 = cell2mat (data.Var6);
data.Var7 = cell2mat (data.Var7);
data.Var8 = cell2mat (data.Var8);
data.Var9 = cell2mat (data.Var9);
data.Var10 = cell2mat (data.Var10);
%%
%%create custom table for each class.
if contains(img_folder,'Baso')
    baso=data;
elseif contains(img_folder,'lymp')
    lymp=data;
elseif contains(img_folder,'eosi')
    eosi=data;
elseif contains(img_folder,'mono')
    mono=data;
elseif contains(img_folder,'neut')
    neut=data;
end

%%
final=[baso;eosi;lymp;mono;neut];
%%
final.Properties.VariableNames={'Area' 'MajorAxisLength' 'Solidity' 'FilledArea' 'MinorAxisLength' 'Circularity' 'ConvexArea' 'Eccentricity' 'Perimeter' 'Label'};
writetable(final,'dataset2.xlsx','Sheet',1);