
% From colorbrewer.org

function cmap = cbrew(n,cname)

if nargin<1, n=256; end
if nargin<2, cname=''; end
qual=0;

switch lower(cname)
    %%% Sequential Schemes
    case 'bugn'
        cmap =[247, 252, 253; 229, 245, 249; 204, 236, 230; 153, 216, 201; 102, 194, 164; 65, 174, 118; 35, 139, 69; 0, 109, 44; 0, 68, 27];
    case 'bupu'
        cmap =[247, 252, 253; 224, 236, 244; 191, 211, 230; 158, 188, 218; 140, 150, 198; 140, 107, 177; 136, 65, 157; 129, 15, 124; 77, 0, 75];
    case 'gnbu'
        cmap =[247, 252, 240; 224, 243, 219; 204, 235, 197; 168, 221, 181; 123, 204, 196; 78, 179, 211; 43, 140, 190; 8, 104, 172; 8, 64, 129];
    case 'orrd'
        cmap =[255, 247, 236; 254, 232, 200; 253, 212, 158; 253, 187, 132; 252, 141, 89; 239, 101, 72; 215, 48, 31; 179, 0, 0; 127, 0, 0];
    case 'pubu'
        cmap =[255, 247, 251; 236, 231, 242; 208, 209, 230; 166, 189, 219; 116, 169, 207; 54, 144, 192; 5, 112, 176; 4, 90, 141; 2, 56, 88];
    case 'pubugn'
        cmap =[255, 247, 251; 236, 226, 240; 208, 209, 230; 166, 189, 219; 103, 169, 207; 54, 144, 192; 2, 129, 138; 1, 108, 89; 1, 70, 54];
    case 'purd'
        cmap =[247, 244, 249; 231, 225, 239; 212, 185, 218; 201, 148, 199; 223, 101, 176; 231, 41, 138; 206, 18, 86; 152, 0, 67; 103, 0, 31];
    case 'rdpu'
        cmap =[255, 247, 243; 253, 224, 221; 252, 197, 192; 250, 159, 181; 247, 104, 161; 221, 52, 151; 174, 1, 126; 122, 1, 119; 73, 0, 106];
    case 'ylgn'
        cmap =[255, 255, 229; 247, 252, 185; 217, 240, 163; 173, 221, 142; 120, 198, 121; 65, 171, 93; 35, 132, 67; 0, 104, 55; 0, 69, 41];
    case 'ylgnbu'
        cmap =[255, 255, 217; 237, 248, 217; 199, 233, 180; 127, 205, 187; 65, 182, 196; 29, 145, 192; 34, 94, 168; 37, 52, 148; 8, 29, 88];
    case 'ylorbr'
        cmap =[255, 255, 229; 255, 247, 188; 254, 227, 145; 254, 196, 79; 254, 153, 41; 236, 112, 20; 204, 76, 2; 153, 52, 4; 102, 37, 6];
    case 'ylorrd'
        cmap =[255, 255, 204; 255, 237, 160; 254, 217, 118; 254, 178, 76; 253, 141, 60; 252, 78, 42; 227, 26, 28; 189, 0, 38; 128, 0, 38];
    case 'blues'
        cmap =[247, 251, 255; 222, 235, 247; 198, 219, 239; 158, 202, 225; 107, 174, 214; 66, 146, 198; 33, 113, 181; 8, 81, 156; 8, 48, 107];
    case 'greens'
        cmap =[247, 252, 245; 229, 245, 224; 199, 233, 192; 161, 217, 155; 116, 196, 118; 65, 171, 93; 35, 139, 69; 0, 109, 44; 0, 68, 27];
    case 'oranges'
        cmap =[255, 245, 235; 254, 230, 206; 253, 208, 162; 253, 174, 107; 253, 141, 60; 241, 105, 19; 217, 72, 1; 166, 54, 3; 127, 39, 4];
    case 'purples'
        cmap =[252, 251, 253; 239, 237, 245; 218, 218, 235; 188, 189, 220; 158, 154, 200; 128, 125, 186; 106, 81, 163; 84, 39, 143; 63, 0, 125];
    case 'reds'
        cmap =[255, 245, 240; 254, 224, 210; 252, 187, 161; 252, 146, 114; 251, 106, 74; 239, 59, 44; 203, 24, 29; 165, 15, 21; 103, 0, 13];
        
        %%% Qualitative Schemes
    case 'paired'
        cmap = [166, 206, 227; 31, 120, 180; 178, 223, 138; 51, 160, 44; 251, 154, 153; 227, 26, 28; 253, 191, 111; 255, 127, 0; 202, 178, 214; 106, 61, 154; 255, 255, 153];
        qual = 1;
    case 'set3'
        cmap = [141, 211, 199; 255, 255, 179; 190, 186, 218; 251, 128, 114; 128, 177, 211; 253, 180, 98; 179, 222, 105; 252, 205, 229; 217, 217, 217; 188, 128, 189; 204, 235, 197];
        qual = 1;
    case 'pastel1'
        cmap = [251, 180, 174; 179, 205, 227; 204, 235, 197; 222, 203, 228; 254, 217, 166; 255, 255, 204; 229, 216, 189; 253, 218, 236; 242, 242, 242];
        qual = 1;
    case 'set1'
        cmap = [228, 26, 28; 55, 126, 184; 77, 175, 74; 152, 78, 163; 255, 127, 0; 255, 255, 51; 166, 86, 40; 247, 129, 191; 153, 153, 153];
        qual = 1;
    case 'accent'
        cmap = [127, 201, 127; 190, 174, 212; 253, 192, 134; 255, 255, 153; 56, 108, 176; 240, 2, 127; 191, 91, 23; 102, 102, 102];
        qual = 1;
    case 'dark2'
        cmap = [27, 158, 119; 217, 95, 2; 117, 112, 179; 231, 41, 138; 102, 166, 30; 230, 171, 2; 166, 118, 29; 102, 102, 102];
        qual = 1;
    case 'pastel2'
        cmap = [179, 226, 205; 253, 205, 172; 203, 213, 232; 244, 202, 228; 230, 245, 201; 255, 242, 174; 241, 226, 204; 204, 204, 204];
        qual = 1;
    case 'set2'
        cmap = [102, 194, 165; 252, 141, 98; 141, 160, 203; 231, 138, 195; 166, 216, 84; 255, 217, 47; 229, 196, 148; 179, 179, 179];
        qual = 1;
        
        %%% Diverging Schemes
    case 'brbg'
        cmap = [84, 48, 5; 140, 81, 10; 191, 129, 45; 223, 129, 125; 246, 232, 195; 245, 245, 245; 199, 234, 229; 128, 205, 193; 53, 151, 143; 1, 102, 94; 0, 60, 48];
    case 'piyg'
        cmap = [142, 1, 82; 197, 27, 125; 222, 119, 174; 241, 182, 218; 253, 224, 239; 247, 247, 247; 230, 245, 208; 184, 225, 134; 127, 188, 65; 77, 146, 33; 39, 100, 25];
    case 'prgn'
        cmap = [64, 0, 75; 118, 42, 131; 153, 112, 171; 194, 165, 207; 231, 212, 232; 247, 247, 247; 217, 240, 211; 166, 219, 160; 90, 174, 97; 27, 120, 55; 0, 68, 27];
    case 'puor'
        cmap = [127, 59, 8; 179, 88, 6; 224, 130, 20; 253, 184, 99; 254, 224, 182; 247, 247, 247; 216, 218, 235; 178, 171, 210; 128, 115, 172; 84, 39, 136; 45, 0, 75];
    case 'rdbu'
        cmap = [103, 0, 31; 178, 24, 43; 214, 96, 77; 244, 165, 130; 253, 219, 199; 247, 247, 247; 209, 229, 240; 146, 197, 222; 67, 147, 195; 33, 102, 172; 5, 48, 97];
    case 'rdgy'
        cmap = [103, 0, 31; 178, 24, 43; 214, 96, 77; 244, 165, 130; 253, 219, 199; 255, 255, 255; 224, 224, 224; 186, 186, 186; 135, 135, 135; 77, 77, 77; 26, 26, 26];
    case 'rdylbu'
        cmap = [165, 0, 38; 215, 48, 39; 244, 109, 67; 253, 174, 97; 254, 224, 144; 255, 255, 191; 224, 243, 248; 171, 217, 233; 116, 173, 209; 69, 117, 180; 49, 54, 149];
    case 'rdylgn'
        cmap = [165, 0, 38; 215, 48, 39; 244, 109, 67; 253, 174, 97; 254, 224, 139; 255, 255, 191; 217, 239, 139; 166, 217, 106; 102, 189, 99; 26, 152, 80; 0, 104, 55];
    case 'paired_div'
        cmap = [227, 26, 28; 255, 127, 0; 51, 160, 44; 31, 120, 180; 106, 61, 154; 227, 26, 28];
        
        %%% External Schemes...
    case 'seismic'
        cmap = seismic(25)*255;
    case 'rainbow'
        cmap = rainbow(25)*255;        
    case 'hyps'
        cmap = hyps(25)*255;
        
    otherwise %spectral
        cmap = [158, 1, 66; 213, 62, 79; 244, 109, 67; 253, 174, 97; 254, 224, 139; 255, 255, 191; 230, 245, 152; 171, 221, 164; 102, 194, 165; 50, 136, 189; 94, 79, 162];
end
cmap=cmap/255;

if n>0
    if ~qual
        % interpolate...
        x=1:size(cmap,1);
        xi=linspace(1,size(cmap,1),n);
        cmap2=zeros(n,3);
        for i=1:3
            cmap2(:,i)=interp1(x,cmap(:,i),xi,'pchip');
        end
        cmap=cmap2;
    else
        % loop for qualitative maps
        idx=mod([1:n]-1,size(cmap,1))+1;
        cmap=cmap(idx,:);
    end
else
    clist = {'brbg','piyg','prgn','puor','rdbu','rdgy','rdylbu','rdylgn','spectral',...
        'bugn','bupu','gnbu','orrd','pubu','pubugn','purd','rdpu','ylgn','ylgnbu','ylorbr','ylorrd','blues','greens','oranges','purples','reds',...
        'seismic','rainbow','hyps','paired_div'...
        'paired','set3','pastel1','set1','accent','dark2','pastel2','set2'};
    for i=1:length(clist)
        %         clist{i}
        subplot(13,3,i)
        if i>30 % qualitative
            image(permute(cbrew(10,clist{i}),[3 1 2]))
        else
            image(permute(cbrew(128,clist{i}),[3 1 2]))
        end
        
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        ylabel(clist{i})
    end
end




function rgb = seismic(n,sh);

% seismic(n) creates a colormap, ranging from dark blue via white to dark red.
% 
% Nico Sneeuw
% Munich, 31/08/94

if nargin == 0, n = 256; end
if nargin<2, sh=0; end

m = ceil(n/3);
top = ones(m,1);
bot = zeros(m,1);
up = (0:m-1)'/m;
down = flipud(up);

r = [bot; up; 1; top; down];
g = [bot; up; 1; down; bot];
b = [up; top; 1; down; bot];
rgb = [r g b];

% rgb-map has size 4m+1 now. The central part will be extracted.

xlarge = 4*m+1-n;
xblue = round(xlarge/2);
xred = xlarge - xblue;
rgb([1:xblue 4*m-xred+2:4*m+1],:) = [];

c = rgb2hsv(rgb);
c(:,1) = c(:,1)+sh;
c(c(:,1)>1,1) = c(c(:,1)>1,1)-1;
rgb = hsv2rgb(c);

function rgb = rainbow(n);

% RAINBOW(n) creates a colormap, ranging from blue via green to red.
% Similar to 'jet', but without the darkening at the ends.
% 
% Nico Sneeuw
% Munich, 30/08/94

if nargin == 0, n = 256; end

m = fix(n/3);
step = 1/m;
ltop = ones(m+1,1);
stop = ones(m,1);
lbot = zeros(m+1,1);
sbot = zeros(m,1);
lup = (0:step:1)'; 
sup = (step/2:step:1)'; 
ldown = (1:-step:0)';
sdown = (1-step/2:-step:0)';
if n - 3*m == 2
   rgb = [lbot lup ltop; sup stop sdown; ltop ldown lbot];
elseif n-3*m == 1 
   rgb = [sbot sup stop; sup stop sdown; ltop ldown lbot];
else 
   rgb = [sbot sup stop; sup stop sdown; stop sdown sbot];
end

function map = hyps(n)

% HYPS hypsometric (height) colormap
% Linear interpolation between the 'corner colors':
% dark green - green - yellow - golden - brown - wheat
%
% Nico Sneeuw Munich 30/03/95
  
if nargin == 0, n = size(get(gcf,'Colormap'),1); end


map = [0 .4 0 ; % dark green
       .5 .7 .3;	 % light green
       1 .9 0;	 % yellow
       .8 .6 .1;	 % golden rod
       .5 .2 .1;	 % brown 
       1 1 .9];	 % wheat

x = linspace(1,n,length(map));
xi = 1:n;
map = interp1(x,map,xi);

map(map>1)=1;
% keyboard