$HistoryLength=0;
out=DirectoryName[$InputFileName];base=DirectoryName[DirectoryName[out]];fd=Import[FileNameJoin[{base,"input","fi_family.json"}],"RawJSON"];
AppendTo[$Path,FileNameJoin[{fd["blade_root"],"BLAddOns","BladeIBP"}]];Get["BladeIBP.m"];
p=FileNameJoin[{out,"cache","e6NativeDot2_dot2_1","ibps"}];
f=FileNames["ids_*.mx",p];n=Total[Length[DeleteCases[Import[#],0]]&/@f];
a=DeleteDuplicates[Join@@(Import/@FileNames["ints_*.mx",p])];
Export[FileNameJoin[{out,"raw_size.json"}],<|"equations"->n,"FI_unknowns"->Length[a],"file_count"->Length[f]|>,"RawJSON"];Print[n," ",Length[a]];Exit[];
