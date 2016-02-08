function [] = write_brainMap(imgFile, userOptions, writeOpts)
 %write file to standard.
returnHere=pwd; 
readFile = writeOpts.template;
subjectMetadataStruct = spm_vol(readFile);

rMapMetadataStruct_nS = subjectMetadataStruct;
rMapMetadataStruct_nS.fname = fullfile(userOptions.rootPath, ['Maps/' writeOpts.name '.img']);
rMapMetadataStruct_nS.descrip =  writeOpts.description;
rMapMetadataStruct_nS.dim = size(imgFile);
gotoDir(userOptions.rootPath, 'Maps/');

spm_write_vol(rMapMetadataStruct_nS, imgFile);
cd(returnHere);
end