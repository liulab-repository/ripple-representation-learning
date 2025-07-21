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
 

% COLOR.infer_1d =  [255 190 122];
% COLOR.infer_2d = [250 127 111];
% COLOR.memory  = [130 176 210];
% COLOR.other = [142 207 201];
% COLOR = struct2table(COLOR);
% 
% COLOR.('1D inference')=  [255 190 122];
% COLOR.('2D inference')= [250 127 111];
% COLOR.('pairwise memory')= [130 176 210];
% COLOR.('compound')= COLOR.('2D inference');
% 
% COLOR.('before') = [0.6510  0.8902 0.8078 ]*256;
% COLOR.('after') = [0.1216 0.7059 0.4706 ]*256;
% 
COLOR.('correct') = [0.4000    0.7608    0.6471]*256;
% % COLOR.('wrong') = [0.9569    0.4275    0.2627]*256;
COLOR.('wrong') = [215,48,39];
% 
COLOR.('ripple') = [174,1,126] ;
% 
% COLOR.('inference_1D') = COLOR.infer_1d;
% COLOR.('inference_2D') = COLOR.infer_2d;
% COLOR.('pairwise_memroy') = COLOR.memory;
COLOR.gray  = [155 155 155];
COLOR.gray_light = COLOR.gray*1.3;
COLOR.gray_dark = COLOR.gray*0.7;
% 
% COLOR.eleMemory = COLOR.memory;
% COLOR.eleInfer  = COLOR.infer_1d;
% 
% % COLOR.sixfold = [.2 .8 .1 ]*256;
COLOR.sixfold = [237,123,108];
% COLOR.significantTime = [178,24,43];
% COLOR.align = [103,169,207;];
% COLOR.misalign = [239,138,98];
% COLOR.feedback = [123,50,148];
% COLOR.interval = [0,136,55];
% COLOR.rest2 = [0.1 .7 .1]*256;
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
% 
% cols1 = [178,24,43
% 239,138,98
% 253,219,199
% 209,229,240
% 103,169,207
% 33,102,172];
% 
cols2 = [142 207 201
         255 190 122
         250 127 111
         130 176 210];

% COLOR.FcompInfer = cols(1,:);
% COLOR.FcompMemory= cols(2,:);
% COLOR.FeleInfer = cols(4,:);
% COLOR.FeleMemory= cols(5,:);
% COLOR.Flearning = cols(6,:);

COLOR.FcompInfer = [255 137 121];
COLOR.FcompMemory= cols2(3,:);
COLOR.FeleInfer = cols2(2,:);
COLOR.FeleMemory= cols2(4,:);
COLOR.Flearning = cols2(1,:);

% COLOR.FcompInfer = ;
% COLOR.FcompMemory= ;
% COLOR.FeleInfer = ;
% COLOR.FeleMemory= ;
% COLOR.Flearning = ;
% 
% function COLOR = original_method
% COLOR.infer_1d =  [255 190 122];
% COLOR.infer_2d = [250 127 111];
% COLOR.memory  = [130 176 210];
% COLOR.other = [142 207 201];
% COLOR = struct2table(COLOR);
% 
% COLOR.('1D inference')=  [255 190 122];
% COLOR.('2D inference')= [250 127 111];
% COLOR.('pairwise memory')= [130 176 210];
% COLOR.('compound')= COLOR.('2D inference');
% 
% COLOR.('before') = [0.6510  0.8902 0.8078 ]*256;
% COLOR.('after') = [0.1216 0.7059 0.4706 ]*256;
% 
% COLOR.('correct') = [0.4000    0.7608    0.6471]*256;
% % COLOR.('wrong') = [0.9569    0.4275    0.2627]*256;
% COLOR.('wrong') = [215,48,39];
% 
% COLOR.('ripple') = [174,1,126] ;
% 
% COLOR.('inference_1D') = COLOR.infer_1d;
% COLOR.('inference_2D') = COLOR.infer_2d;
% COLOR.('pairwise_memroy') = COLOR.memory;
% COLOR.gray  = [155 155 155];
% COLOR.gray_light = COLOR.gray*1.3;
% COLOR.gray_dark = COLOR.gray*0.7;
% 
% COLOR.eleMemory = COLOR.memory;
% COLOR.eleInfer  = COLOR.infer_1d;
% 
% COLOR.sixfold = [.2 .8 .1 ]*256;
% COLOR.significantTime = [178,24,43];
% COLOR.align = [103,169,207;];
% COLOR.misalign = [239,138,98];
% COLOR.feedback = [123,50,148];
% COLOR.interval = [0,136,55];
% COLOR.rest2 = [0.1 .7 .1]*256;
% 
% % 7 network
% COLOR.Y7_Visual  =[0.46875    0.070312     0.52344]*256;
% COLOR.Y7_Default =[0.80078     0.24219     0.30469]*256;
% COLOR.Y7_DorsalAttention  =[      0     0.46094    0.054688]*256;
% COLOR.Y7_Frontoparietal  =[0.89844     0.57812     0.13281]*256;
% COLOR.Y7_VentralAttention  =[0.76562     0.22656     0.97656]*256;
% COLOR.Y7_Limbic  =[0.85938     0.96875     0.64062]*256;
% COLOR.Y7_Somatomotor  =[0.27344     0.50781     0.70312]*256;

 
 



% [130 176 210]/256  ; [255 190 122]/256 ; [250 127 111]/256

