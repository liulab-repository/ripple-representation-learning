function COLOR = f_deduction_color(scale)
% color map for dedution results visulization
% colorblind safe
 


COLOR = paper_method1;
% COLOR = original_method;

if nargin>0
    COLOR.Variables = COLOR.Variables/256;
end


function COLOR = paper_method1
COLOR = table;
COLOR.white = [1 1 1 ]*256;
 
COLOR.('correct') = [0.4000    0.7608    0.6471]*256;
COLOR.('wrong') = [215,48,39];
% 
COLOR.('ripple') = [174,1,126] ;
% 
COLOR.gray  = [155 155 155];
COLOR.gray_light = COLOR.gray*1.3;
COLOR.gray_dark = COLOR.gray*0.7;
COLOR.sixfold = [237,123,108];
COLOR.rest  = [136 197 192];
 
% 7 network
COLOR.Y7_Visual  =[0.46875    0.070312     0.52344]*256;
COLOR.Y7_Default =[0.80078     0.24219     0.30469]*256;
COLOR.Y7_DorsalAttention  =[      0     0.46094    0.054688]*256;
COLOR.Y7_Frontoparietal  =[0.89844     0.57812     0.13281]*256;
COLOR.Y7_VentralAttention  =[0.76562     0.22656     0.97656]*256;
COLOR.Y7_Limbic  =[0.85938     0.96875     0.64062]*256;
COLOR.Y7_Somatomotor  =[0.27344     0.50781     0.70312]*256;

% set by figure 
COLOR.FcompInfer = [215,25,28];
COLOR.FcompMemory= [253,174,97];
COLOR.FeleInfer = [255,255,191];
COLOR.FeleMemory= [171,217,233];
COLOR.Flearning = [44,123,182];
COLOR.Falign    = [237,123,108];

COLOR.FadjacentInfer = [237,123,108];
COLOR.FnonadjacentInfer = [249,212,200];
 
% 
cols2 = [142 207 201
         255 190 122
         250 127 111
         130 176 210];
 
COLOR.FcompInfer = [255 137 121];
COLOR.FcompMemory= cols2(3,:);
COLOR.FeleInfer = cols2(2,:);
COLOR.FeleMemory= cols2(4,:);
COLOR.Flearning = cols2(1,:);
 