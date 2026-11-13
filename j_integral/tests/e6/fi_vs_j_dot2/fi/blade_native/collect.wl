$HistoryLength=0;out=DirectoryName[$InputFileName];base=DirectoryName[DirectoryName[out]];
fd=Import[FileNameJoin[{base,"input","fi_family.json"}],"RawJSON"];
AppendTo[$Path,FileNameJoin[{fd["blade_root"],"BLAddOns","BladeIBP"}]];Get["BladeIBP.m"];
files=FileNames["indep",FileNameJoin[{out,"cache"}],Infinity];
records=Table[Module[{dir,eqfiles,eqs,selected,ints,selectedints,info,allfiles,all,indep},
dir=DirectoryName[DirectoryName[file]];eqfiles=FileNames["ids_*.mx",FileNameJoin[{dir,"ibps"}]];
If[eqfiles==={},Return[<|"directory"->dir,"error"->"no numeric equation files"|>,Module]];
eqs=DeleteCases[Join@@(Import/@eqfiles),0];indep=Get[file];selected=eqs[[indep]];
ints=IntegralsIn[eqs];selectedints=IntegralsIn[selected];info=Get[FileNameJoin[{dir,"results","datainfo"}]];
allfiles=FileNames["ints_*.mx",FileNameJoin[{dir,"ibps"}]];all=DeleteDuplicates[Join@@(Import/@allfiles)];
<|"directory"->FileNameTake[dir],"raw_equations"->Length[eqs],"raw_FI_unknowns"->Length[ints],"all_integral_file_entries"->Length[all],"selected_equations"->Length[indep],"selected_FI_unknowns"->Length[selectedints],"returned_dependent_integrals"->Length[info[[1]]],"masters"->Length[info[[2]]]|>
],{file,files}];
Export[FileNameJoin[{out,"system_sizes.json"}],<|"systems"->records,"system_count"->Length[records]|>,"RawJSON"];Print[records];Exit[];
