BeginPackage["DELocalStandardization`"];

StandardizeDEAtPoint::usage =
  "StandardizeDEAtPoint[A,x,a] maps x=a to t=0, reduces the system to " <>
  "Fuchs form, and shifts residue eigenvalues to the standard strip. " <>
  "Use a=Infinity for t=1/x.";

Options[StandardizeDEAtPoint] = {
  "LocalVariable" -> Automatic,
  "DESolverPath" -> Automatic
};

Begin["`Private`"];

ClearAll[zeroMatrixQ];
zeroMatrixQ[matrix_] := And @@ (PossibleZeroQ /@ Flatten[Together[matrix]]);

StandardizeDEAtPoint[
  matrix_, variable_Symbol, point_, OptionsPattern[]
] := Module[
  {localVariable, solverPath, localMatrix, coordinateDefinition,
   transformation, inverseTransformation, standardizedMatrix,
   residue, originalRank, standardizedRank},
  localVariable = Replace[
    OptionValue["LocalVariable"],
    Automatic :> DESolver`eta
  ];
  solverPath = OptionValue["DESolverPath"];
  If[solverPath =!= Automatic && !NameQ["DESolver`NormalizeMat"],
    Get[solverPath]
  ];
  If[!NameQ["DESolver`NormalizeMat"],
    Return[Failure[
      "DESolverNotLoaded",
      <|"Message" -> "Load ref/DERun/DESolver.m or pass DESolverPath."|>
    ]]
  ];
  If[localVariable =!= DESolver`eta,
    Return[Failure[
      "UnsupportedLocalVariable",
      <|"Message" -> "The current DERun implementation requires DESolver`eta."|>
    ]]
  ];
  If[point === Infinity,
    localMatrix = Together[
      -(matrix /. variable -> 1/localVariable)/localVariable^2
    ];
    coordinateDefinition = HoldForm[localVariable == 1/variable],
    localMatrix = Together[
      matrix /. variable -> point + localVariable
    ];
    coordinateDefinition = HoldForm[localVariable == variable - point]
  ];
  originalRank = DESolver`Private`PoincareRank[localMatrix];
  {
    transformation,
    inverseTransformation,
    standardizedMatrix
  } = DESolver`NormalizeMat[localMatrix];
  standardizedRank = DESolver`Private`PoincareRank[standardizedMatrix];
  residue = Together[localVariable standardizedMatrix] /.
    localVariable -> 0;
  <|
    "Point" -> point,
    "LocalVariable" -> localVariable,
    "CoordinateDefinition" -> coordinateDefinition,
    "LocalMatrix" -> localMatrix,
    "OriginalPoincareRank" -> originalRank,
    "Transformation" -> transformation,
    "InverseTransformation" -> inverseTransformation,
    "StandardizedMatrix" -> standardizedMatrix,
    "StandardizedPoincareRank" -> standardizedRank,
    "Residue" -> residue,
    "ResidueEigenvalues" -> Eigenvalues[residue],
    "ResidueEigenvalueCounts" -> Counts[Eigenvalues[residue]],
    "InverseCheckQ" -> zeroMatrixQ[
      inverseTransformation . transformation -
        IdentityMatrix[Length[matrix]]
    ],
    "DifferentialEquationCheckQ" -> zeroMatrixQ[
      inverseTransformation . (
        localMatrix . transformation -
          D[transformation, localVariable]
      ) - standardizedMatrix
    ]
  |>
];

End[];
EndPackage[];
