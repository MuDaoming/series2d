BeginPackage["DEReductionPoleConstraints`"];

ComputePoleConstraintStructures::usage =
  "ComputePoleConstraintStructures[A, x, x0, maxDI] constructs the exact " <>
  "Taylor-jet solution space at x0 and the joint pole-numerator structures " <>
  "for dI=1,...,maxDI. The working assumption is that all negative integer " <>
  "powers and logarithmic integer-region modes have been removed.";

PoleTupleContainedQ::usage =
  "PoleTupleContainedQ[tuple, structure] tests whether a joint numerator " <>
  "tuple (c_dI,...,c_1) belongs to a computed pole-constraint structure.";

Options[ComputePoleConstraintStructures] = {
  "Modulus" -> 0,
  "EquationOrderStart" -> Automatic,
  "EquationOrderStep" -> 2,
  "EquationOrderLimit" -> 24,
  "StableRepeats" -> 2
};

Begin["`Private`"];

ClearAll[
  fieldScalar, fieldMatrix, fieldRowReduce, fieldNullSpace, nonzeroRowQ,
  rowBasis, polynomialValuation, localPoleOrder, buildTaylorSystem,
  projectedJetBasis, buildToeplitz
];

fieldScalar[value_, 0] := Together[value];
fieldScalar[value_, modulus_Integer] := Module[
  {rational = Together[value], numerator, denominator},
  numerator = Mod[Numerator[rational], modulus];
  denominator = Mod[Denominator[rational], modulus];
  Mod[
    numerator PowerMod[denominator, -1, modulus],
    modulus
  ]
];

fieldMatrix[matrix_, modulus_] :=
  Map[fieldScalar[#, modulus] &, matrix, {2}];

fieldRowReduce[matrix_, 0] := RowReduce[matrix];
fieldRowReduce[matrix_, modulus_Integer] :=
  RowReduce[matrix, Modulus -> modulus];

fieldNullSpace[matrix_, 0] := NullSpace[matrix];
fieldNullSpace[matrix_, modulus_Integer] :=
  NullSpace[matrix, Modulus -> modulus];

nonzeroRowQ[row_, 0] := !And @@ (PossibleZeroQ /@ row);
nonzeroRowQ[row_, modulus_Integer] :=
  !And @@ (Mod[#, modulus] === 0 & /@ row);

rowBasis[matrix_, modulus_] := If[
  matrix === {} || Length[matrix] === 0,
  {},
  Select[
    fieldRowReduce[matrix, modulus],
    nonzeroRowQ[#, modulus] &
  ]
];

polynomialValuation[0, variable_] := Infinity;
polynomialValuation[polynomial_, variable_] :=
  Exponent[Expand[polynomial], variable, Min];

localPoleOrder[matrix_, variable_, point_, localVariable_] := Module[
  {local, valuations},
  local = Together[matrix /. variable -> point + localVariable];
  valuations = Map[
    Function[entry,
      If[
        entry === 0,
        Infinity,
        polynomialValuation[Numerator[Together[entry]], localVariable] -
          polynomialValuation[
            Denominator[Together[entry]],
            localVariable
          ]
      ]
    ],
    local,
    {2}
  ];
  Max[0, -Min[Flatten[valuations]]]
];

buildTaylorSystem[
  coefficients_Association,
  dimension_Integer,
  poleOrder_Integer,
  equationOrder_Integer,
  modulus_
] := Module[
  {identity, zero, maximumJet, coefficient, blockRows},
  identity = IdentityMatrix[dimension];
  zero = ConstantArray[0, {dimension, dimension}];
  maximumJet = equationOrder + Max[poleOrder, 1];
  coefficient[index_] := Lookup[coefficients, index, zero];
  blockRows = Table[
    ArrayFlatten[{
      Table[
        fieldMatrix[
          coefficient[power - jet] -
            If[
              power >= 0 && jet === power + 1,
              (power + 1) identity,
              zero
            ],
          modulus
        ],
        {jet, 0, maximumJet}
      ]
    }],
    {power, -poleOrder, equationOrder}
  ];
  <|
    "MaximumJetOrder" -> maximumJet,
    "Matrix" -> Join @@ blockRows
  |>
];

projectedJetBasis[
  systemMatrix_,
  retainedColumns_Integer,
  modulus_
] := Module[{null, projected},
  null = fieldNullSpace[systemMatrix, modulus];
  If[null === {}, Return[{}]];
  projected = null[[All, 1 ;; retainedColumns]];
  rowBasis[projected, modulus]
];

buildToeplitz[jetBasisRows_, dimension_, dI_, modulus_] := Module[
  {restrictedJetBasis, solutionColumns, parameterCount, blocks},
  restrictedJetBasis = rowBasis[
    If[
      jetBasisRows === {},
      {},
      jetBasisRows[[All, 1 ;; dI dimension]]
    ],
    modulus
  ];
  parameterCount = Length[restrictedJetBasis];
  If[parameterCount === 0,
    Return[ConstantArray[0, {dI dimension, 0}]]
  ];
  solutionColumns = Transpose[restrictedJetBasis];
  blocks = Table[
    solutionColumns[[
      order dimension + 1 ;; (order + 1) dimension,
      All
    ]],
    {order, 0, dI - 1}
  ];
  ArrayFlatten@Table[
    If[column >= row, blocks[[column - row + 1]],
      ConstantArray[0, {dimension, parameterCount}]
    ],
    {row, 1, dI},
    {column, 1, dI}
  ]
];

ComputePoleConstraintStructures[
  matrix_,
  variable_Symbol,
  point_,
  maximumDI_Integer?Positive,
  OptionsPattern[]
] := Module[
  {
    modulus, equationStart, equationStep, equationLimit, stableRepeats,
    dimension, localVariable, poleOrder, localMatrix, seriesMatrix,
    coefficients, equationOrder, system, jetBasis, previousBasis = Missing[],
    stableCount = 0, history = {}, stabilizedQ = False, retainedColumns,
    structures, toeplitz, allowedBasis, rank, dI,
    singleLayerBasis, productBasis, productRank, unionRank
  },

  modulus = OptionValue["Modulus"];
  equationStep = OptionValue["EquationOrderStep"];
  equationLimit = OptionValue["EquationOrderLimit"];
  stableRepeats = OptionValue["StableRepeats"];
  dimension = Length[matrix];
  localVariable = Unique["localX"];
  poleOrder = localPoleOrder[
    matrix, variable, point, localVariable
  ];
  equationStart = Replace[
    OptionValue["EquationOrderStart"],
    Automatic :>
      Max[0, maximumDI - 1 - Max[poleOrder, 1]]
  ];
  retainedColumns = maximumDI dimension;

  localMatrix = Together[
    matrix /. variable -> point + localVariable
  ];
  seriesMatrix = Map[
    Normal@Series[#, {localVariable, 0, equationLimit}] &,
    localMatrix,
    {2}
  ];
  coefficients = Association@Table[
    order -> fieldMatrix[
      Map[
        Coefficient[#, localVariable, order] &,
        seriesMatrix,
        {2}
      ],
      modulus
    ],
    {order, -poleOrder, equationLimit}
  ];

  jetBasis = {};
  Do[
    system = buildTaylorSystem[
      coefficients,
      dimension,
      poleOrder,
      equationOrder,
      modulus
    ];
    jetBasis = projectedJetBasis[
      system["Matrix"],
      retainedColumns,
      modulus
    ];
    AppendTo[
      history,
      <|
        "EquationOrder" -> equationOrder,
        "MaximumJetOrder" -> system["MaximumJetOrder"],
        "RetainedJetDimension" -> Length[jetBasis]
      |>
    ];
    If[
      !MissingQ[previousBasis] && jetBasis === previousBasis,
      stableCount++,
      stableCount = 0
    ];
    If[stableCount >= stableRepeats,
      stabilizedQ = True;
      Break[]
    ];
    previousBasis = jetBasis,
    {
      equationOrder,
      equationStart,
      equationLimit,
      equationStep
    }
  ];

  structures = Association@Table[
    toeplitz = buildToeplitz[
      jetBasis,
      dimension,
      dI,
      modulus
    ];
    rank = If[
      Dimensions[toeplitz][[2]] === 0,
      0,
      Length[rowBasis[toeplitz, modulus]]
    ];
    allowedBasis = If[
      Dimensions[toeplitz][[2]] === 0,
      IdentityMatrix[dI dimension],
      fieldNullSpace[Transpose[toeplitz], modulus]
    ];
    dI -> <|
      "dI" -> dI,
      "Modulus" -> modulus,
      "TaylorJetDimension" ->
        Length[rowBasis[
          If[
            jetBasis === {},
            {},
            jetBasis[[All, 1 ;; dI dimension]]
          ],
          modulus
        ]],
      "ToeplitzMatrix" -> toeplitz,
      "ToeplitzRank" -> rank,
      "AmbientDimension" -> dI dimension,
      "AllowedDimension" -> Length[allowedBasis],
      "AllowedJointNumeratorBasis" -> allowedBasis
    |>,
    {dI, 1, maximumDI}
  ];

  singleLayerBasis =
    structures[1]["AllowedJointNumeratorBasis"];
  Do[
    productBasis = Join @@ Table[
      Map[
        Function[basisRow,
          Flatten@Table[
            If[
              block === activeBlock,
              basisRow,
              ConstantArray[0, dimension]
            ],
            {block, 1, dI}
          ]
        ],
        singleLayerBasis
      ],
      {activeBlock, 1, dI}
    ];
    productRank = Length[rowBasis[productBasis, modulus]];
    unionRank = Length[rowBasis[
      Join[
        structures[dI]["AllowedJointNumeratorBasis"],
        productBasis
      ],
      modulus
    ]];
    structures[dI] = Append[
      structures[dI],
      "EqualsLayerwiseD1ProductQ" ->
        (
          productRank === structures[dI]["AllowedDimension"] &&
          unionRank === structures[dI]["AllowedDimension"]
        )
    ],
    {dI, 1, maximumDI}
  ];

  <|
    "Point" -> point,
    "Dimension" -> dimension,
    "MaximumDI" -> maximumDI,
    "Modulus" -> modulus,
    "DEPoleOrder" -> poleOrder,
    "StabilizedQ" -> stabilizedQ,
    "StabilityHistory" -> history,
    "TaylorJetBasisRows" -> jetBasis,
    "Structures" -> structures
  |>
];

PoleTupleContainedQ[tuple_, structure_Association] := Module[
  {modulus, residual},
  modulus = Lookup[structure, "Modulus", 0];
  residual = tuple . structure["ToeplitzMatrix"];
  If[
    modulus === 0,
    And @@ (PossibleZeroQ /@ residual),
    And @@ (Mod[#, modulus] === 0 & /@ residual)
  ]
];

End[];
EndPackage[];
