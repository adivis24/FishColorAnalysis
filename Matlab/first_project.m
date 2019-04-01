
img_dir = 'C:\Users\adith\Downloads\FrogAnalysis';%image folde directory
cd (img_dir)

img_list_v=[dir(['*ropped*/*v.*']); dir(['*ropped*/*Ven*']);dir(['*ropped*/* V *'])]
img_list_d=[dir(['*ropped*/*d.*']); dir(['*ropped*/*Dor*']);dir(['*ropped*/* D *'])]
v={img_list_v.name} %list all "ventral" images.
v_num= numel(v)
d={img_list_d.name} %list all "dorsal" images.
d_num= numel(d)
