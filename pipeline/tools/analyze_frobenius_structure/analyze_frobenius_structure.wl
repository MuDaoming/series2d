BeginPackage["LocalFrobeniusStructure`"];

AnalyzeLocalFrobeniusStructure::usage =
  "AnalyzeLocalFrobeniusStructure[A, x, x0] uses DERun's " <>
  "FindAsymptoticBehavior and CalcTaylor2Num to determine the true local " <>
  "power-log structures and the first appearances of all free parameters.";

CompareSimplePoleWithResidue::usage =
  "CompareSimplePoleWithResidue[A, x, x0, analysis] compares the true " <>
  "Frobenius exponents found by AnalyzeLocalFrobeniusStructure with the " <>
  "eigenvalues of Res[(x-x0) A, x=x0].";

Options[AnalyzeLocalFrobeniusStructure] = {
  "Masters" -> Automatic,
  "XOrder" -> 8,
  "ExtraXOrder" -> 12,
  "WorkingPrecision" -> 120,
  "ChopPrecision" -> 50,
  "EvaluationOffset" -> 1/100
};

Begin["`Private`"];

toolDirectory = DirectoryName[$InputFileName];
repoRoot = ExpandFileName[
  FileNameJoin[{toolDirectory, "..", "..", ".."}]
];

If[!NameQ["DESolver`FindAsymptoticBehavior"],
  Get[FileNameJoin[{repoRoot, "ref", "DERun", "DESolver.m"}]]
];
If[!NameQ["DESolver`FindPowerLog"],
  Get[FileNameJoin[{repoRoot, "ref", "DERun", "frobenius.m"}]]
];

ClearAll[numericallyZeroQ, numericalRank, observationOrder];

numericallyZeroQ[value_, tolerance_] :=
  Max[Abs[N[Flatten[{value}]]]] <= tolerance;

numericalRank[rows_, tolerance_] := If[
  rows === {},
  0,
  Module[{singularValues, scale},
    singularValues = SingularValueList[N[rows]];
    scale = Max[singularValues];
    Count[singularValues, value_ /; value > tolerance Max[1, scale]]
  ]
];

observationOrder[item_] := {
  N[item["Power"]],
  -item["LogPower"],
  item["IntegralIndex"],
  item["CandidateIndex"],
  item["SeriesOrder"]
};

AnalyzeLocalFrobeniusStructure[
  matrix_,
  variable_Symbol,
  point_,
  OptionsPattern[]
] := Module[
  {
    masters, xOrder, extraXOrder, workingPrecision, chopPrecision,
    evaluationOffset, tolerance, dimension, localMatrix, nheq, nheqn,
    candidateBehavior, expansions, outputIntegralOrder, metadata,
    observations, vector,
    recurrenceFreeParameters, totalOrder, sp, reduction, unsolved,
    basisIndices, block, beh, logk, fid,
    orderedObservations, pivotObservations = {}, pivotRows = {},
    oldRank = 0, newRank, completeQ, pivotMatrix, inversePivotMatrix,
    transformedObservations, freeParameters, parameterProfiles,
    perIntegralStructures, leadingForParameterAndIntegral, candidates,
    minimumPower, leadingCandidates, chosen, parameterIndex,
    integralIndex, parameterCount, actualExponentCounts
  },

  dimension = Length[matrix];
  masters = Replace[
    OptionValue["Masters"],
    Automatic :> ("M" <> ToString[#] & /@ Range[dimension])
  ];
  If[Length[masters] =!= dimension,
    Return[Failure[
      "InvalidMasters",
      <|"Message" -> "The number of master labels must equal the matrix dimension."|>
    ]]
  ];

  xOrder = OptionValue["XOrder"];
  extraXOrder = OptionValue["ExtraXOrder"];
  workingPrecision = OptionValue["WorkingPrecision"];
  chopPrecision = OptionValue["ChopPrecision"];
  evaluationOffset = OptionValue["EvaluationOffset"];
  (* CalcTaylor2Num contains a few machine-precision structural decisions.
     Do not interpret their roundoff as additional Frobenius modes. *)
  tolerance = 10^-Min[chopPrecision, 12];
  localMatrix = matrix /. variable -> DESolver`eta;

  Block[
    {
      DESolver`Private`WorkingPre = workingPrecision,
      DESolver`Private`ChopPre = chopPrecision,
      DESolver`Private`SilentMode = True,
      DESolver`Private`XOrder = xOrder,
      DESolver`Private`ExtraXOrder = extraXOrder,
      $MinPrecision = workingPrecision,
      $MaxPrecision = Infinity
    },
    DESolver`Private`SetWorkingPrecision[workingPrecision];
    nheq = DESolver`Private`NHEquations2[localMatrix, point, "Taylor"];
    nheqn = DESolver`NHEquationsNum[nheq];
    candidateBehavior =
      DESolver`FindAsymptoticBehavior[localMatrix, point];
    expansions = Table[
      DESolver`Private`CalcTaylor2Num[
        nheq,
        nheqn,
        UnitVector[dimension, boundaryIndex],
        N[evaluationOffset, workingPrecision],
        candidateBehavior
      ],
      {boundaryIndex, dimension}
    ];

    (* Repeat CalcTaylor2Num's homogeneous-equation basis/fid analysis
       exactly, but retain the information that CalcTaylor2Num normally
       keeps private.  These are the canonical recurrence-level first
       appearances; unlike a pivot basis chosen from the completed
       fundamental matrix, they are not changed by propagation between
       blocks. *)
    totalOrder = xOrder + extraXOrder;
    recurrenceFreeParameters = Flatten[
      Table[
        block = nheq[[blockIndex, 4]];
        beh = candidateBehavior[[blockIndex, behaviorIndex]];
        logk = beh[[2]];
        sp = DESolver`Private`GenHomoEqsDemo[
          nheq[[blockIndex, 1]],
          nheq[[blockIndex, 2]],
          beh,
          totalOrder
        ];
        reduction = First[
          DESolver`Private`SparseGaussianNum[
            sp,
            ConstantArray[0, Total[Length /@ sp[[2]]]]
          ]
        ];
        fid[index_] := {
          block[[Mod[index - 1, Length[block]] + 1]],
          Mod[
            Quotient[index - 1, Length[block]],
            logk + 1
          ],
          totalOrder + 1 -
            Quotient[
              index - 1,
              (logk + 1) Length[block]
            ]
        };
        unsolved = Complement[
          Range[sp[[1, 2]]],
          Keys[reduction][[All, 1]]
        ];
        basisIndices = Select[
          unsolved,
          0 <= fid[#][[3]] <= xOrder &
        ];
        Map[
          Function[index,
            With[{position = fid[index]},
              <|
                "IntegralIndex" -> position[[1]],
                "Master" -> masters[[position[[1]]]],
                "LogPower" -> position[[2]],
                "CandidateMu" -> beh[[1]],
                "CandidateSeriesOrder" -> position[[3]],
                "TrueMu" -> beh[[1]] + position[[3]]
              |>
            ]
          ],
          basisIndices
        ],
        {blockIndex, Length[nheq]},
        {
          behaviorIndex,
          Length[candidateBehavior[[blockIndex]]]
        }
      ],
      2
    ];
  ];

  (* CalcTaylor2Num returns integrals in diagonal-block order, not in the
     original matrix order. *)
  outputIntegralOrder = Flatten[nheqn[[All, 5]]];

  metadata = Flatten[
    Table[
      With[
        {
          integralIndex = outputIntegralOrder[[outputPosition]],
          rule = expansions[[1, outputPosition, candidateIndex]],
          candidateMu =
            First[expansions[[1, outputPosition, candidateIndex]]],
          logSeries =
            Last[expansions[[1, outputPosition, candidateIndex]]]
        },
        Flatten[
          Table[
            <|
              "OutputPosition" -> outputPosition,
              "IntegralIndex" -> integralIndex,
              "Master" -> masters[[integralIndex]],
              "CandidateIndex" -> candidateIndex,
              "CandidateMu" -> candidateMu,
              "LogPower" -> logPower,
              "SeriesOrder" -> seriesOrder,
              "Power" -> candidateMu + seriesOrder
            |>,
            {logPower, 0, Length[logSeries] - 1},
            {seriesOrder, 0, Length[logSeries[[logPower + 1]]] - 1}
          ],
          1
        ]
      ],
      {outputPosition, dimension},
      {candidateIndex, Length[expansions[[1, outputPosition]]]}
    ],
    2
  ];

  observations = Map[
    Function[item,
      vector = Table[
        Last[
          expansions[[
            boundaryIndex,
            item["OutputPosition"],
            item["CandidateIndex"]
          ]]
        ][[
          item["LogPower"] + 1,
          item["SeriesOrder"] + 1
        ]],
        {boundaryIndex, dimension}
      ];
      Append[item, "BoundaryCoefficientVector" -> vector]
    ],
    metadata
  ];
  observations = Select[
    observations,
    !numericallyZeroQ[#["BoundaryCoefficientVector"], tolerance] &
  ];

  orderedObservations = SortBy[observations, observationOrder];
  Do[
    newRank = numericalRank[
      Append[pivotRows, item["BoundaryCoefficientVector"]],
      tolerance
    ];
    If[newRank > oldRank,
      AppendTo[pivotRows, item["BoundaryCoefficientVector"]];
      AppendTo[pivotObservations, item];
      oldRank = newRank;
    ],
    {item, orderedObservations}
  ];

  parameterCount = Length[pivotRows];
  completeQ = parameterCount === dimension;
  If[!completeQ,
    Return[<|
      "Point" -> point,
      "Dimension" -> dimension,
      "CandidateBehavior" -> candidateBehavior,
      "FreeParameterCount" -> parameterCount,
      "CompleteQ" -> False,
      "Message" ->
        "The retained coefficient rows do not yet span the full solution space. " <>
        "Increase XOrder, ExtraXOrder, or WorkingPrecision.",
      "PivotObservations" -> pivotObservations
    |>]
  ];

  pivotMatrix = pivotRows;
  inversePivotMatrix = Inverse[pivotMatrix];
  transformedObservations = Map[
    Append[
      #,
      "FreeParameterCoefficientVector" ->
        Chop[#["BoundaryCoefficientVector"] . inversePivotMatrix, tolerance]
    ] &,
    observations
  ];

  freeParameters = MapIndexed[
    Function[{item, index},
      <|
        "Parameter" -> ("C" <> ToString[index[[1]]]),
        "TrueMu" -> item["Power"],
        "LogPower" -> item["LogPower"],
        "FirstAppearanceIntegralIndex" -> item["IntegralIndex"],
        "FirstAppearanceMaster" -> item["Master"],
        "CandidateMu" -> item["CandidateMu"],
        "CandidateSeriesOrder" -> item["SeriesOrder"]
      |>
    ],
    pivotObservations
  ];

  leadingForParameterAndIntegral[parameterIndex_, integralIndex_] :=
    Module[{available, minPower, atMinPower, maxLog, atMaxLog},
      available = Select[
        transformedObservations,
        #["IntegralIndex"] === integralIndex &&
        !numericallyZeroQ[
          #["FreeParameterCoefficientVector"][[parameterIndex]],
          tolerance
        ] &
      ];
      If[available === {}, Return[Missing["Absent"]]];
      minPower = Min[available[[All, "Power"]]];
      atMinPower = Select[available, #["Power"] === minPower &];
      maxLog = Max[atMinPower[[All, "LogPower"]]];
      atMaxLog = Select[atMinPower, #["LogPower"] === maxLog &];
      chosen = First@SortBy[atMaxLog, observationOrder];
      <|
        "Parameter" -> ("C" <> ToString[parameterIndex]),
        "Mu" -> minPower,
        "LogPower" -> maxLog,
        "CandidateMu" -> chosen["CandidateMu"],
        "FirstSeriesOrder" -> chosen["SeriesOrder"]
      |>
    ];

  parameterProfiles = Table[
    <|
      "Parameter" -> ("C" <> ToString[parameterIndex]),
      "GlobalStructure" -> <|
        "Mu" -> freeParameters[[parameterIndex, "TrueMu"]],
        "LogPower" -> freeParameters[[parameterIndex, "LogPower"]]
      |>,
      "IntegralLeadingStructures" -> Association@Table[
        masters[[integralIndex]] ->
          leadingForParameterAndIntegral[parameterIndex, integralIndex],
        {integralIndex, dimension}
      ]
    |>,
    {parameterIndex, parameterCount}
  ];

  perIntegralStructures = Table[
    candidates = DeleteMissing[
      leadingForParameterAndIntegral[#, integralIndex] & /@
        Range[parameterCount]
    ];
    <|
      "IntegralIndex" -> integralIndex,
      "Master" -> masters[[integralIndex]],
      "Structures" -> Map[
        Function[group,
          <|
            "Mu" -> group[[1, "Mu"]],
            "LogPower" -> group[[1, "LogPower"]],
            "Parameters" -> group[[All, "Parameter"]]
          |>
        ],
        GatherBy[candidates, {#["Mu"], #["LogPower"]} &]
      ]
    |>,
    {integralIndex, dimension}
  ];

  actualExponentCounts = KeySort@Counts[
    freeParameters[[All, "TrueMu"]]
  ];

  <|
    "Point" -> point,
    "Dimension" -> dimension,
    "Masters" -> masters,
    "XOrder" -> xOrder,
    "ExtraXOrder" -> extraXOrder,
    "CandidateBehavior" -> candidateBehavior,
    "FreeParameterCount" -> parameterCount,
    "CompleteQ" -> completeQ,
    "ActualExponentCounts" -> actualExponentCounts,
    "RecurrenceFreeParameters" ->
      MapIndexed[
        Prepend[
          #1,
          "Parameter" -> ("R" <> ToString[#2[[1]]])
        ] &,
        recurrenceFreeParameters
      ],
    "RecurrenceExponentCounts" ->
      KeySort@Counts[recurrenceFreeParameters[[All, "TrueMu"]]],
    "FreeParameters" -> freeParameters,
    "PerIntegralStructures" -> perIntegralStructures,
    "ParameterProfiles" -> parameterProfiles
  |>
];

CompareSimplePoleWithResidue[
  matrix_,
  variable_Symbol,
  point_,
  analysis_Association
] := Module[
  {residue, residueEigenvalues, residueCounts, frobeniusCounts},
  residue = FullSimplify[
    Limit[(variable - point) matrix, variable -> point]
  ];
  residueEigenvalues = Eigenvalues[residue];
  residueCounts = KeySort@Counts[residueEigenvalues];
  frobeniusCounts = analysis["ActualExponentCounts"];
  <|
    "Point" -> point,
    "Residue" -> residue,
    "ResidueEigenvalueCounts" -> residueCounts,
    "FrobeniusExponentCounts" -> frobeniusCounts,
    "MatchQ" -> (residueCounts === frobeniusCounts)
  |>
];

End[];
EndPackage[];
