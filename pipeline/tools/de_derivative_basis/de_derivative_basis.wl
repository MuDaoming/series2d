BeginPackage["DEDerivativeBasis`"];

ReadApartReductionRows::usage =
  "ReadApartReductionRows[path, masters, x, p] reconstructs rational rows.";
AnalyzeDerivativeBasisGrid::usage =
  "AnalyzeDerivativeBasisGrid[A,x,targets,rValues,sValues] tests " <>
  "polynomial derivative-basis ansatzes over a finite field.";

Options[AnalyzeDerivativeBasisGrid] = {
  "Modulus" -> 0,
  "TrainingExtraPoints" -> 4,
  "HoldoutPoints" -> 12,
  "SamplePointStart" -> 101,
  "Verbose" -> True
};

Begin["`Private`"];

ClearAll[
  fieldScalar, polynomialFromCoefficients, apartRowExpression,
  trimCoefficients, rationalData, polynomialTaylor, rationalTaylor,
  matrixData, matrixTaylorAtPoint, targetValuesAtPoint,
  derivativeMatricesAtPoint, goodPointRecord, collectPointRecords,
  designMatrix, solveTarget, analyzePair
];

fieldScalar[value_, modulus_Integer] := Module[
  {rational = Together[value], numerator, denominator},
  numerator = Mod[Numerator[rational], modulus];
  denominator = Mod[Denominator[rational], modulus];
  If[denominator === 0, Return[$Failed]];
  Mod[numerator PowerMod[denominator, -1, modulus], modulus]
];

polynomialFromCoefficients[coefficients_List, variable_] :=
  Sum[coefficients[[index + 1]] variable^index,
    {index, 0, Length[coefficients] - 1}];

apartRowExpression[row_Association, variable_] := Module[
  {numerator, factor},
  numerator = polynomialFromCoefficients[
    ToExpression[row["numerator_coeffs"]], variable
  ];
  Switch[
    row["kind"],
    "POLY", numerator,
    "POLE",
      factor = polynomialFromCoefficients[
        ToExpression[row["factor_coeffs"]], variable
      ];
      numerator/factor^ToExpression[row["exponent"]],
    _, 0
  ]
];

ReadApartReductionRows[
  path_String, masters_List, variable_Symbol, modulus_Integer
] := Module[{table, rows, selected, names},
  table = Import[path, "TSV"];
  If[Length[table] < 2, Return[<||>]];
  rows = AssociationThread[First[table], #] & /@ Rest[table];
  selected = Select[rows, MemberQ[masters, #["master"]] &];
  names = DeleteDuplicates[Lookup[selected, "name"]];
  Association@Table[
    name -> Table[
      Together@Total[
        apartRowExpression[#, variable] & /@
          Select[selected,
            #["name"] === name && #["master"] === master &]
      ],
      {master, masters}
    ],
    {name, names}
  ]
];

trimCoefficients[list_List] := Module[{result = list},
  While[Length[result] > 1 && Last[result] === 0, result = Most[result]];
  result
];

rationalData[expression_, variable_Symbol, modulus_Integer] := Module[
  {rational, numerator, denominator},
  rational = Together[expression];
  numerator = trimCoefficients[
    fieldScalar[#, modulus] & /@
      CoefficientList[Numerator[rational], variable]
  ];
  denominator = trimCoefficients[
    fieldScalar[#, modulus] & /@
      CoefficientList[Denominator[rational], variable]
  ];
  <|"Numerator" -> numerator, "Denominator" -> denominator|>
];

matrixData[matrix_, variable_Symbol, modulus_Integer] :=
  Map[rationalData[#, variable, modulus] &, matrix, {2}];

polynomialTaylor[
  coefficients_List, point_Integer, order_Integer, modulus_Integer
] := Table[
  Mod[
    Sum[
      coefficients[[degree + 1]] Binomial[degree, power] *
        PowerMod[point, degree - power, modulus],
      {degree, power, Length[coefficients] - 1}
    ],
    modulus
  ],
  {power, 0, order}
];

rationalTaylor[
  data_Association, point_Integer, order_Integer, modulus_Integer
] := Module[{numerator, denominator, quotient, degree},
  numerator = polynomialTaylor[
    data["Numerator"], point, order, modulus
  ];
  denominator = polynomialTaylor[
    data["Denominator"], point, order, modulus
  ];
  If[denominator[[1]] === 0, Return[$Failed]];
  quotient = ConstantArray[0, order + 1];
  quotient[[1]] = Mod[
    numerator[[1]] PowerMod[denominator[[1]], -1, modulus], modulus
  ];
  Do[
    quotient[[degree + 1]] = Mod[
      (
        numerator[[degree + 1]] -
        Sum[
          denominator[[index + 1]] *
            quotient[[degree - index + 1]],
          {index, 1, degree}
        ]
      ) PowerMod[denominator[[1]], -1, modulus],
      modulus
    ],
    {degree, 1, order}
  ];
  quotient
];

matrixTaylorAtPoint[
  data_, point_Integer, order_Integer, modulus_Integer
] := Module[{entryJets, dimensions},
  dimensions = Dimensions[data];
  entryJets = Map[
    rationalTaylor[#, point, order, modulus] &, data, {2}
  ];
  If[!FreeQ[entryJets, $Failed], Return[$Failed]];
  Table[
    Table[
      entryJets[[row, column, power + 1]],
      {row, dimensions[[1]]}, {column, dimensions[[2]]}
    ],
    {power, 0, order}
  ]
];

targetValuesAtPoint[
  data_Association, point_Integer, modulus_Integer
] := Module[{values},
  values = Map[
    Function[dataRow,
      Map[
        Function[entry,
          With[{jet = rationalTaylor[entry, point, 0, modulus]},
            If[jet === $Failed, $Failed, First[jet]]
          ]
        ],
        dataRow
      ]
    ],
    data
  ];
  If[FreeQ[values, $Failed], values, $Failed]
];

derivativeMatricesAtPoint[
  aJets_List, maximumDerivative_Integer, modulus_Integer
] := Module[{dimension, taylorMatrices, next, n, k},
  dimension = Length[First[aJets]];
  taylorMatrices = {IdentityMatrix[dimension]};
  Do[
    next = ConstantArray[0, {dimension, dimension}];
    Do[
      next = Mod[
        next + aJets[[k + 1]] . taylorMatrices[[n - k + 1]],
        modulus
      ],
      {k, 0, n}
    ];
    next = Mod[next PowerMod[n + 1, -1, modulus], modulus];
    AppendTo[taylorMatrices, next],
    {n, 0, maximumDerivative - 1}
  ];
  Table[
    Mod[Factorial[n] taylorMatrices[[n + 1]], modulus],
    {n, 0, maximumDerivative}
  ]
];

goodPointRecord[
  aData_, targetData_Association, point_Integer,
  maximumDerivative_Integer, modulus_Integer
] := Module[{aJets, targetValues},
  aJets = matrixTaylorAtPoint[
    aData, point, Max[0, maximumDerivative - 1], modulus
  ];
  If[aJets === $Failed, Return[$Failed]];
  targetValues = targetValuesAtPoint[targetData, point, modulus];
  If[targetValues === $Failed, Return[$Failed]];
  <|
    "Point" -> point,
    "DerivativeMatrices" ->
      derivativeMatricesAtPoint[aJets, maximumDerivative, modulus],
    "TargetValues" -> targetValues
  |>
];

collectPointRecords[
  aData_, targetData_Association, count_Integer, start_Integer,
  maximumDerivative_Integer, modulus_Integer, verbose_
] := Module[{records = {}, candidate = start, record, lastReport = 0},
  While[Length[records] < count,
    record = goodPointRecord[
      aData, targetData, Mod[candidate, modulus],
      maximumDerivative, modulus
    ];
    If[record =!= $Failed, AppendTo[records, record]];
    candidate++;
    If[
      verbose && Length[records] >= lastReport + 10,
      lastReport = 10 Floor[Length[records]/10];
      Print["  prepared points: ", Length[records], "/", count]
    ];
  ];
  records
];

designMatrix[
  records_List, derivativeOrder_Integer, polynomialDegree_Integer,
  modulus_Integer
] := Join @@ Map[
  Function[record,
    Module[{basisRows = {}, point, matrices},
      point = record["Point"];
      matrices = record["DerivativeMatrices"];
      Do[
        AppendTo[
          basisRows,
          Mod[
            PowerMod[point, power, modulus] *
              matrices[[order + 1, row]],
            modulus
          ]
        ],
        {order, 0, derivativeOrder},
        {row, Length[First[matrices]]},
        {power, 0, polynomialDegree}
      ];
      Transpose[basisRows]
    ]
  ],
  records
];

solveTarget[
  design_, holdoutDesign_, trainingRecords_List, holdoutRecords_List,
  targetName_, derivativeOrder_Integer, polynomialDegree_Integer,
  dimension_Integer, variable_Symbol, modulus_Integer
] := Module[
  {rhs, holdoutRhs, solution, trainingResidual, holdoutResidual,
   tensor, polynomialRows, passed},
  rhs = Flatten[
    Lookup[#["TargetValues"], targetName] & /@ trainingRecords
  ];
  holdoutRhs = Flatten[
    Lookup[#["TargetValues"], targetName] & /@ holdoutRecords
  ];
  solution = Quiet[
    Check[LinearSolve[design, rhs, Modulus -> modulus], $Failed]
  ];
  If[
    solution === $Failed || !VectorQ[solution, IntegerQ],
    Return[<|"PassedQ" -> False, "Status" -> "no_solution"|>]
  ];
  solution = Mod[solution, modulus];
  trainingResidual = Mod[design . solution - rhs, modulus];
  holdoutResidual = Mod[holdoutDesign . solution - holdoutRhs, modulus];
  passed = And @@ (# === 0 & /@
    Join[trainingResidual, holdoutResidual]);
  tensor = ArrayReshape[
    solution,
    {derivativeOrder + 1, dimension, polynomialDegree + 1}
  ];
  polynomialRows = Table[
    Table[
      Sum[
        tensor[[order + 1, row, power + 1]] variable^power,
        {power, 0, polynomialDegree}
      ],
      {row, dimension}
    ],
    {order, 0, derivativeOrder}
  ];
  <|
    "PassedQ" -> passed,
    "Status" -> If[passed, "passed", "holdout_failed"],
    "NonzeroCoordinateCount" ->
      Count[solution, value_ /; value =!= 0],
    "CoordinateTensor" -> tensor,
    "PolynomialRows" -> polynomialRows
  |>
];

analyzePair[
  records_List, targetNames_List, derivativeOrder_Integer,
  polynomialDegree_Integer, trainingExtra_Integer,
  holdoutCount_Integer, dimension_Integer, variable_Symbol,
  modulus_Integer, verbose_
] := Module[
  {unknownCount, trainingCount, trainingRecords, holdoutRecords,
   design, holdoutDesign, rank, fits, passedNames},
  unknownCount =
    dimension (derivativeOrder + 1) (polynomialDegree + 1);
  trainingCount = Ceiling[unknownCount/dimension] + trainingExtra;
  trainingRecords = Take[records, trainingCount];
  holdoutRecords = Take[
    records, {trainingCount + 1, trainingCount + holdoutCount}
  ];
  design = designMatrix[
    trainingRecords, derivativeOrder, polynomialDegree, modulus
  ];
  holdoutDesign = designMatrix[
    holdoutRecords, derivativeOrder, polynomialDegree, modulus
  ];
  rank = MatrixRank[design, Modulus -> modulus];
  If[verbose,
    Print[
      "  (r,s)=(", derivativeOrder, ",", polynomialDegree,
      ") unknowns=", unknownCount, " rank=", rank
    ]
  ];
  fits = AssociationMap[
    solveTarget[
      design, holdoutDesign, trainingRecords, holdoutRecords, #,
      derivativeOrder, polynomialDegree, dimension, variable, modulus
    ] &,
    targetNames
  ];
  passedNames = Select[targetNames, TrueQ[fits[#]["PassedQ"]] &];
  <|
    "DerivativeOrder" -> derivativeOrder,
    "PolynomialDegree" -> polynomialDegree,
    "UnknownCount" -> unknownCount,
    "DesignRank" -> rank,
    "KernelDimension" -> unknownCount - rank,
    "TrainingPointCount" -> trainingCount,
    "HoldoutPointCount" -> holdoutCount,
    "PassedCount" -> Length[passedNames],
    "FailedCount" -> Length[targetNames] - Length[passedNames],
    "PassedNames" -> passedNames,
    "FailedNames" -> Complement[targetNames, passedNames],
    "Fits" -> fits
  |>
];

AnalyzeDerivativeBasisGrid[
  matrix_, variable_Symbol, targets_Association,
  derivativeOrders_List, polynomialDegrees_List, OptionsPattern[]
] := Module[
  {modulus, trainingExtra, holdoutCount, start, verbose, dimension,
   targetNames, maximumDerivative, maximumDegree, maximumUnknowns,
   maximumTraining, totalPointCount, aData, targetData, records,
   pairs, results},
  modulus = OptionValue["Modulus"];
  If[!IntegerQ[modulus] || modulus <= 1,
    Print["AnalyzeDerivativeBasisGrid currently requires Modulus -> p."];
    Return[$Failed]
  ];
  trainingExtra = OptionValue["TrainingExtraPoints"];
  holdoutCount = OptionValue["HoldoutPoints"];
  start = OptionValue["SamplePointStart"];
  verbose = TrueQ[OptionValue["Verbose"]];
  dimension = Length[matrix];
  targetNames = Keys[targets];
  maximumDerivative = Max[derivativeOrders];
  maximumDegree = Max[polynomialDegrees];
  maximumUnknowns =
    dimension (maximumDerivative + 1) (maximumDegree + 1);
  maximumTraining =
    Ceiling[maximumUnknowns/dimension] + trainingExtra;
  totalPointCount = maximumTraining + holdoutCount;
  If[verbose,
    Print[
      "Preparing data: N=", dimension,
      ", targets=", Length[targetNames],
      ", max r=", maximumDerivative,
      ", max s=", maximumDegree
    ]
  ];
  aData = matrixData[matrix, variable, modulus];
  targetData = Map[
    Map[rationalData[#, variable, modulus] &, #] &,
    targets
  ];
  If[verbose, Print["Collecting ", totalPointCount, " good points"]];
  records = collectPointRecords[
    aData, targetData, totalPointCount, start,
    maximumDerivative, modulus, verbose
  ];
  pairs = Tuples[{derivativeOrders, polynomialDegrees}];
  results = Association@Table[
    ToString[pair[[1]]] <> "," <> ToString[pair[[2]]] ->
      analyzePair[
        records, targetNames, pair[[1]], pair[[2]], trainingExtra,
        holdoutCount, dimension, variable, modulus, verbose
      ],
    {pair, pairs}
  ];
  <|
    "Dimension" -> dimension,
    "TargetCount" -> Length[targetNames],
    "TargetNames" -> targetNames,
    "DerivativeOrders" -> derivativeOrders,
    "PolynomialDegrees" -> polynomialDegrees,
    "Modulus" -> modulus,
    "SamplePoints" -> Lookup[records, "Point"],
    "Results" -> results
  |>
];

End[];
EndPackage[];
