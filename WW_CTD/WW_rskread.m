function WW_rskread(WWmeta)
% read .rsk file downloaded from rbr
% bz, June 15, 2021


if length(WWmeta.rbrfile)>2;
    fprintf('Watch out \nThere is more than one rsk file\n')
    for j=1:length(filedir); disp(WWmeta.rbrfile(j).name);end
end
fprintf('read rbr file is %s\n',WWmeta.rbrfile(1).name)

disp('RSK_wrapper--- It is may take a while --- coffee time buddy?')

RSKfile= fullfile(WWmeta.rbrpath,WWmeta.rbrfile(1).name);
RSKdb=RSKopen(RSKfile);
RSKread=RSKreaddata(RSKdb);


save(fullfile(WWmeta.rbrpath,[WWmeta.name_rbr '.mat']),'RSKread');
end