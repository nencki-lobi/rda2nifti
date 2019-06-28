function rda2nii(rda_filename)
[dir,fn,ext]=fileparts(rda_filename);
read_rda102;
clc
voxel=[rda.FoVWidth,rda.FoVHeight,rda.SliceThickness];
r0=rda.PositionVector-cross(rda.RowVector,rda.ColumnVector)*voxel(3)/2;
center=[rda.VOIPositionSag rda.VOIPositionCor rda.VOIPositionTra];

center_spm=center.*[1,1,-1];
r0_spm=r0.*[1,1,-1];
M=-[rda.RowVector.*[1,1,-1];rda.ColumnVector.*[1,1,-1];cross(rda.RowVector,rda.ColumnVector).*[1,1,-1]];
nii = make_nii(ones(voxel), [1,1,1], r0_spm+1,[]); % img , voxel size [mm], center [vox], datatype 
save_nii(nii,'tmp.nii')
nii=load_untouch_nii('tmp.nii');
nii.hdr.hist.srow_x(1:3)=M(:,1);
nii.hdr.hist.srow_y(1:3)=M(:,2);
nii.hdr.hist.srow_z(1:3)=M(:,3);
save_untouch_nii(nii,fullfile(dir,[fn,'.nii']));

%check
diff_v=r0_spm-center_spm
dir_v=M(1,:)*voxel(1)/2+M(2,:)*voxel(2)/2+M(3,:)*voxel(3)/2
if(diff_v-dir_v<0.001)
    fprintf('\tvertex to center - OK\n')
end

end