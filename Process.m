function [output,w,s,ref] = Process(ref_num,i,imscc,TE) %,R2starmap,fatfraction2] = Process(ref_num,i,imscc,TE)
Size = size(imscc);
len = Size(1); width = Size(2);
Slice_num = Size(4);  Rnum = Size(3); TEnum = Size(6);

% weight threshold?
ref = squeeze(imscc(:,:,ref_num,i,1,:));
%%  h
ss = 10;
slice_window = 3;
block_size = 3;
search_window = 3;
block_num = 1;
[output,w,s] = nlmslice(ref,imscc,i,block_size,search_window,block_num,slice_window,ref_num,ss);
% [outParams2,fatfraction2] = processFatWater(reshape(output,len,width,1,1,1,TEnum),TE(1),TE(2) - TE(1),TEnum);
% R2starmap = outParams2.r2starmap;
