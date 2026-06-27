(* ::Package:: *)

If[Length[$ScriptCommandLine] != 3,
  Print["Usage: wolframscript -file generate_symmetry.wl <family_input.wl> <output_dir>"];
  Exit[1]
];

inputPath = $ScriptCommandLine[[2]];
outputDir = $ScriptCommandLine[[3]];

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..",
  "symmetry", "wolfram", "SectorSymmetry.wl"}]];

Get[inputPath];

If[! ValueQ[polynomials] || Length[polynomials] != 2,
  Print["Input must define polynomials = {U,F}."];
  Exit[1]
];

If[! ValueQ[parameterHead], parameterHead = x];
If[! ValueQ[branches],
  Print["Input must define branches = {b1,b2,...}."];
  Exit[1]
];

If[! DirectoryQ[outputDir], CreateDirectory[outputDir, CreateIntermediateDirectories -> True]];

Clear[fullBranchQ, mapPairs, crossBranchQ, twoBranchOnlyQ, vecText, ruleText];

fullBranchQ[sector_List] := Module[{active},
  active = DeleteDuplicates[Pick[branches, sector, 1]];
  Length[active] == Length[DeleteDuplicates[branches]]
];

mapPairs[rules_List] := Cases[
  rules,
  HoldPattern[Rule[a_[i_Integer], b_[j_Integer]]] /; a === parameterHead && b === parameterHead :> (i -> j)
];

crossBranchQ[rules_List] := Module[{pairs = mapPairs[rules]},
  AnyTrue[pairs, branches[[#[[1]]]] =!= branches[[#[[2]]]] &]
];

twoBranchOnlyQ[src_List, dst_List] := ! fullBranchQ[src] || ! fullBranchQ[dst];

vecText[v_List] := "{" <> StringRiffle[ToString /@ v, ", "] <> "}";
ruleText[pairs_List] := "{" <> StringRiffle[(ToString[#[[1]]] <> " -> " <> ToString[#[[2]]]) & /@ pairs, ", "] <> "}";

result = SectorSymmetry`FindSectorRelations[polynomials, parameterHead];
relations = Select[result["Relations"], #["Sector"] =!= #["Representative"] &];

accepted = {};
rejected = {};

Do[
  src = rel["Sector"];
  dst = rel["Representative"];
  rules = rel["ParameterRules"];
  pairs = mapPairs[rules];
  reason = Which[
    ! TrueQ[rel["UVerified"]] || ! TrueQ[rel["FVerified"]], "verification_failed",
    twoBranchOnlyQ[src, dst], "two_branch_only",
    crossBranchQ[rules], "cross_branch_propagator_map",
    True, ""
  ];
  If[reason === "",
    AppendTo[accepted, {src, dst, pairs}],
    AppendTo[rejected, <|
      "Source" -> src,
      "Target" -> dst,
      "Map" -> pairs,
      "Reason" -> reason
    |>]
  ],
  {rel, relations}
];

Export[
  FileNameJoin[{outputDir, "sectormap"}],
  StringRiffle[
    ("{" <> vecText[#[[1]]] <> " -> " <> vecText[#[[2]]] <> ", " <> ruleText[#[[3]]] <> "}") & /@ accepted,
    "\n"
  ] <> If[Length[accepted] > 0, "\n", ""],
  "Text"
];

Export[
  FileNameJoin[{outputDir, "symmetry_rejected"}],
  StringRiffle[
    ("source=" <> vecText[#["Source"]] <>
      " target=" <> vecText[#["Target"]] <>
      " map=" <> ruleText[#["Map"]] <>
      " reason=" <> #["Reason"]) & /@ rejected,
    "\n"
  ] <> If[Length[rejected] > 0, "\n", ""],
  "Text"
];

orbits = GatherBy[
  Join[
    ({#[[1]], #[[2]]} &) /@ accepted,
    ({#[[2]], #[[2]]} &) /@ accepted
  ],
  Last
];

Export[
  FileNameJoin[{outputDir, "symmetry_summary"}],
  StringRiffle[
    {
      "accepted=" <> ToString[Length[accepted]],
      "rejected=" <> ToString[Length[rejected]],
      "parameters=" <> ToString[result["Parameters"] // InputForm],
      "representatives=" <> ToString[result["Representatives"] // InputForm]
    },
    "\n"
  ] <> "\n",
  "Text"
];

Print["accepted = ", Length[accepted]];
Print["rejected = ", Length[rejected]];
