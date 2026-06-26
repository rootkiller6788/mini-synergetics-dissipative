-- Self-Organization Far-from-Equilibrium -- Lean 4 Formalization
-- Formalizes: thermodynamic regimes, entropy production, bifurcations,
--   order parameters, reaction networks, pattern formation conditions

inductive Regime : Type where
  | equilibrium | near_equilibrium | nonlinear_noneq
  | far_equilibrium | critical | chaotic_noneq
  deriving Repr, DecidableEq

def Regime.toString : Regime -> String
  | .equilibrium => "Equilibrium"
  | .near_equilibrium => "Near-Equilibrium"
  | .nonlinear_noneq => "Nonlinear Non-Equilibrium"
  | .far_equilibrium => "Far-Equilibrium"
  | .critical => "Critical"
  | .chaotic_noneq => "Chaotic Non-Equilibrium"

inductive BifurcationType : Type where
  | none | saddle_node | pitchfork_super | pitchfork_sub
  | hopf_super | hopf_sub | transcritical | turing
  deriving Repr, DecidableEq

inductive StabilityClass : Type where
  | stable_node | stable_focus | unstable_node
  | unstable_focus | saddle | center
  deriving Repr, DecidableEq

inductive SymmetryType : Type where
  | none | z2 | o2 | so3 | translational | gauge | chiral | scale
  deriving Repr, DecidableEq

structure EntropyProduction where
  totalProduction : Float
  entropyFlow : Float
  internalProduction : Float
  nonNegativity : Bool := totalProduction >= 0.0
  deriving Repr

structure OnsagerMatrix (n : Nat) where
  coefficients : List (List Float)
  symmetric : Bool := True
  deriving Repr

structure ChemicalSpecies where
  name : String
  concentration : Float
  diffusionCoefficient : Float
  isAutocatalytic : Bool := false
  isInhibitor : Bool := false
  deriving Repr

structure Reaction where
  name : String
  rateConstant : Float
  activationEnergy : Float
  reactants : List (Nat * Float)
  products : List (Nat * Float)
  kineticOrder : Nat := 1
  deriving Repr

structure ReactionNetwork where
  name : String
  species : List ChemicalSpecies
  reactions : List Reaction
  temperature : Float := 298.15
  volume : Float := 1.0
  deriving Repr

structure SpatialGrid where
  nx : Nat
  ny : Nat := 1
  nz : Nat := 1
  dx : Float
  dy : Float := 1.0
  dz : Float := 1.0
  periodicX : Bool := true
  periodicY : Bool := true
  periodicZ : Bool := true
  deriving Repr

structure OrderParameter where
  name : String
  nComponents : Nat
  amplitude : List Float
  brokenSymmetry : SymmetryType
  isSlowVariable : Bool := true
  relaxationTime : Float := 1.0
  deriving Repr

def landauFreeEnergy (psi a b Tc T : Float) : Float :=
  let r := a * (T - Tc)
  0.5 * r * psi * psi + 0.25 * b * psi * psi * psi * psi

def landauEquilibriumAmplitude (a b Tc T : Float) : Float :=
  let r := a * (T - Tc)
  if r < 0.0 && b > 0.0 then Float.sqrt (-r / b) else 0.0

theorem landau_amplitude_satisfies_equilibrium (psi a b Tc T : Float)
    (h_eq : psi = landauEquilibriumAmplitude a b Tc T)
    (h_b_positive : b > 0.0) : psi * (a * (T - Tc) + b * psi * psi) = 0.0 := by
  subst h_eq
  unfold landauEquilibriumAmplitude
  simp
  split
  · exact rfl
  · exact rfl

structure GPStabilityCriterion where
  excessEntropyProduction : Float
  excessEntropy : Float
  isStable : Bool := excessEntropyProduction <= 0.0
  deriving Repr

def entropyProductionDecomposition (dS_dt deS_dt diS_dt : Float) : Bool :=
  dS_dt == deS_dt + diS_dt && diS_dt >= 0.0

def massActionRate (k : Float) (concentrations stoichiometries : List Float) : Float :=
  let pairs := List.zip concentrations stoichiometries
  List.foldl (fun acc (c, nu) => acc * (c ^ nu.toUInt64.toNat)) k pairs

def stoichiometricEntry (productStoich reactantStoich : Float) : Float :=
  productStoich - reactantStoich

def brusselatorRates (A B x y : Float) : Float * Float :=
  (A - (B + 1.0) * x + x * x * y, B * x - x * x * y)

def schloglRate (k1 k2 k3 k4 A B x : Float) : Float :=
  k1 * A * x * x - k2 * x * x * x - k3 * x + k4 * B

structure TuringCondition where
  homogeneousStable : Bool
  spatialUnstable : Bool
  activatorSlowerThanInhibitor : Bool
  isTuringUnstable : Bool := homogeneousStable && spatialUnstable && activatorSlowerThanInhibitor
  deriving Repr

def oregonatorRates (k1 k2 k3 k4 k5 A B f x y z : Float) : Float * Float * Float :=
  (k1 * A * y - k2 * x * y + k3 * A * x - 2.0 * k4 * x * x,
   -k1 * A * y - k2 * x * y + f * k5 * B * z,
   2.0 * k3 * A * x - k5 * B * z)

def ginzburgLandauRate (Gamma a b Tc T psi : Float) : Float :=
  -Gamma * (a * (T - Tc) * psi + b * psi * psi * psi)

def isChaotic (lambdaMax : Float) : Bool := lambdaMax > 0.0

def kaplanYorkeDimension (exponents : List Float) : Float :=
  let rec go (remaining : List Float) (acc_sum : Float) (k : Nat) : Float :=
    match remaining with
    | [] => k.toFloat
    | lam :: rest =>
      let new_sum := acc_sum + lam
      if new_sum >= 0.0 then go rest new_sum (k + 1)
      else k.toFloat + acc_sum / (-lam)
  go exponents 0.0 0

def fluctuationDissipationRelation (kT omega imag_chi : Float) : Float :=
  (2.0 * kT / omega) * imag_chi

def jarzynskiEqualityEstimate (workValues : List Float) (kT : Float) : Float :=
  let n := workValues.length.toFloat
  let avg_exp := List.foldl (fun acc w => acc + Float.exp (-w / kT)) 0.0 workValues / n
  -kT * Float.log avg_exp

def thermodynamicUncertaintyRelation
    (currentMean currentVariance entropyProdRate kB : Float) : Bool :=
  let lhs := currentVariance / (currentMean * currentMean)
  let rhs := 2.0 * kB / entropyProdRate
  lhs >= rhs
