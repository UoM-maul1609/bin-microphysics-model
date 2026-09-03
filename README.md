# Bin Microphysics Model (BMM)

BMM is a detailed sectional/bin cloud microphysics model for studying aerosol activation, warm-cloud evolution, ice formation, mixed-phase microphysics, collisions, secondary ice production, sedimentation/fallout, particle optics, entrainment and cloud-chamber processes.

The model is written primarily in Fortran and writes NetCDF output. Python tools under `python/` provide plotting and experiment-specific workflows.

BMM is intended as a **process model and numerical laboratory** rather than a bulk microphysics parameterisation. A central design goal is to retain aerosol composition and particle history while hydrometeors grow, evaporate, freeze, sublimate, collide and are remapped between numerical bins.

> **Research-code status**
>
> BMM is actively developed. The model contains several alternative numerical and physical treatments specifically so they can be compared. Results should therefore be reported with the relevant namelist settings and BMM revision.

---

## Contents

- [Overview](#overview)
- [Main scientific capabilities](#main-scientific-capabilities)
- [Particle and aerosol representation](#particle-and-aerosol-representation)
- [Liquid water, hygroscopic growth and activation](#liquid-water-hygroscopic-growth-and-activation)
- [Ice nucleation](#ice-nucleation)
- [Ice growth and habit](#ice-growth-and-habit)
- [Collisions and secondary ice](#collisions-and-secondary-ice)
- [Bin representations](#bin-representations)
- [Aerosol release during evaporation and sublimation](#aerosol-release-during-evaporation-and-sublimation)
- [Entrainment and cloud dynamics](#entrainment-and-cloud-dynamics)
- [Cloud-chamber model](#cloud-chamber-model)
- [AIDAd/iSKYLAB interpretation of the chamber treatment](#aidadiskylab-interpretation-of-the-chamber-treatment)
- [Particle loss to surfaces and fallout](#particle-loss-to-surfaces-and-fallout)
- [Optical properties](#optical-properties)
- [Building BMM](#building-bmm)
- [Running BMM](#running-bmm)
- [Namelist structure](#namelist-structure)
- [NetCDF output](#netcdf-output)
- [Python analysis](#python-analysis)
- [Repository layout](#repository-layout)
- [Development and numerical checking](#development-and-numerical-checking)
- [Selected references](#selected-references)

---

# Overview

BMM follows populations of aerosol-containing particles through changes in water content and phase.

Aerosol material is not discarded when a cloud droplet forms or an ice particle nucleates. Instead, nonvolatile component masses remain attached to the particle and are transferred through:

- hygroscopic growth;
- cloud activation;
- condensation and evaporation;
- freezing and melting;
- vapour deposition and sublimation;
- collision/coalescence;
- aggregation;
- secondary ice production;
- numerical remapping;
- entrainment;
- complete hydrometeor evaporation or sublimation.

This makes it possible to ask questions that are difficult to address with a bulk scheme, for example:

- Which dry aerosol sizes activated?
- Which aerosol populations subsequently froze?
- How does collision/coalescence alter aerosol provenance?
- What happens to residual aerosol after complete droplet evaporation?
- How sensitive is ice production to the numerical bin treatment?
- How do homogeneous and inhomogeneous mixing alter the particle-size distribution for the same bulk liquid-water loss?
- How rapidly must chamber air mix to reproduce laboratory cloud evolution?
- How important are loss processes to chamber walls, fan-driven circulation and the chamber floor?

The main executable can be used for atmospheric parcel/bubble/plume calculations or for prescribed chamber experiments.

---

# Main scientific capabilities

| Process | BMM capability |
|---|---|
| Aerosol | Multiple externally mixed modes, lognormal internal submodes and multiple composition components |
| Hygroscopic growth | Molecular/Raoult Köhler and κ-Köhler |
| Adsorption | FHH adsorption treatment |
| Cloud activation | Critical point of the current Köhler/FHH equilibrium curve |
| Warm microphysics | Condensation and evaporation with coupled thermodynamics |
| Ice nucleation | Koop homogeneous freezing, INAS/IASD, DeMott and Daily/DCMEX options |
| INP memory | Prognostic cumulative freezing-threshold moments |
| DeMott aerosol basis | Prognostic number of aerosol cores with dry diameter > 0.5 µm |
| Ice growth | Vapour deposition/sublimation with evolving habit and density |
| Collisions | Stochastic collection equation (SCE) |
| Secondary ice | Hallett-Mossop, mode-1/mode-2 treatments and collisional breakup |
| Bin numerics | Fully moving, moving-centre and Chen-Lamb conservative remapping |
| Entrainment | Adiabatic, lateral/vertical entrainment, homogeneous and inhomogeneous mixing |
| Chamber mode | Prescribed pressure/temperature/water forcing and optional wall-BL processing |
| Chamber evaporation | Homogeneous shrinkage or inhomogeneous complete-particle evaporation |
| Chamber wall water | Legacy closure plus finite-reservoir and finite-rate mass-transfer alternatives |
| Surface losses | Fan-associated particle loss, non-gravitational wall deposition and fallout |
| Optics | Geometric approximation or anomalous diffraction theory |
| Output | NetCDF bulk and sectional diagnostics |

---

# Particle and aerosol representation

## External modes and internal submodes

The aerosol initial condition is specified as `n_mode` externally mixed modes. Each external mode can contain `n_intern` lognormal submodes.

Typical aerosol inputs include:

```fortran
n_aer1
d_aer1
sig_aer1
mass_frac_aer1
```

where:

- `n_aer1` is aerosol number mixing ratio;
- `d_aer1` is number-median dry diameter;
- `sig_aer1` is `ln(sigma_g)`;
- `mass_frac_aer1` defines the nonvolatile composition.

Component properties include:

```fortran
molw_core1
density_core1
nu_core1
kappa_core1
afhh_core1
bfhh_core1
inp_category
```

A particle can therefore contain several nonvolatile aerosol components while remaining part of one dynamically evolving hydrometeor population.

## Extensive and inherited moments

In addition to particle number and water/ice mass, BMM carries several particle moments.

Some are ordinary extensive quantities that add during collisions, while others use specialised inheritance rules. Examples include:

- nonvolatile aerosol component masses;
- ice habit/aggregation properties;
- cumulative INAS threshold populations;
- number of DeMott-origin primary ice monomers;
- number of DeMott-eligible aerosol cores with dry diameter greater than 0.5 µm.

This moment-based architecture is what allows aerosol and ice-nucleation history to survive growth, remapping and collisions.

---

# Liquid water, hygroscopic growth and activation

## Köhler formulations

`kappa_flag` selects the hygroscopic activity treatment:

```text
0 = molecular/Raoult Köhler
1 = κ-Köhler
```

The molecular treatment uses component molecular weight and effective dissociation number.

The κ-Köhler treatment uses component `kappa_core1`.

Both formulations include the Kelvin term and can be combined with FHH adsorption where configured.

## FHH adsorption

Components with non-zero:

```fortran
afhh_core1
bfhh_core1
```

can take up water through the Frenkel-Halsey-Hill adsorption treatment.

This is useful for aerosol such as mineral dust where adsorption may be important before conventional Köhler activation.

BMM distinguishes adsorption-grown haze from activated droplets.

## Activation

Activation is diagnosed from the maximum of the current equilibrium growth curve.

A wet aerosol particle is therefore not considered a cloud droplet merely because it contains water. It becomes activated after passing the critical point of the appropriate Köhler/FHH equilibrium curve.

This activation state is also used by heterogeneous immersion-freezing treatments.

---

# Ice nucleation

Ice microphysics is enabled with:

```fortran
ice_flag = 1
```

The four primary ice mechanisms are selected with:

```fortran
ice_nucleation_mech(1:4)
```

in the order:

```text
1 = Koop homogeneous freezing
2 = INAS / discrete active-site spectrum
3 = DeMott
4 = Daily/DCMEX
```

For one mixed aerosol population, heterogeneous schemes use the precedence:

```text
INAS > DeMott > Daily
```

so a population is not independently counted by several empirical heterogeneous schemes.

Koop homogeneous freezing is treated separately.

---

## Discrete INAS / IASD treatment

BMM can represent a continuous temperature-dependent ice-active-site spectrum with a set of prognostic cumulative threshold classes.

The number of classes is controlled by:

```fortran
n_inp_classes
```

and their temperatures by:

```fortran
inp_temp
```

For example:

```fortran
n_inp_classes = 8
inp_temp(1:8) = -5.0, -8.0, -11.0, -14.0, -17.0, -20.0, -23.0, -26.0
```

There is no requirement that the model use only the historical 8- or 16-class grids. The number of classes is configurable.

The threshold temperatures are common, but the fraction of aerosol assigned to each cumulative class depends on:

- aerosol component;
- component surface area;
- the configured `inp_category`;
- the corresponding `n_s(T)` or IASD formulation.

Supported categories in the current source include:

```text
none
niemand12
kaolinite_murray11
kfeldspar_atkinson13
demott
daily25
```

The active-site spectrum is created for fresh aerosol and then transported prognostically.

It is **not recalculated from the current mixed composition after every collision**.

### INAS inheritance through collisions

For a cumulative threshold with active fractions `F1` and `F2` in two colliding parent populations, the collision product uses:

\[
F_{\rm new}=1-(1-F_1)(1-F_2).
\]

Thus a collision product inherits an active site if either parent carried one.

This is handled through a specialised inherited-moment type rather than simple mass averaging.

---

## DeMott treatment and the >0.5 µm aerosol moment

The DeMott parameterisation depends on the concentration of eligible aerosol with dry diameter greater than 0.5 µm.

BMM carries this quantity explicitly as a **prognostic extensive particle moment**.

For every hydrometeor population, the model tracks the number of DeMott-eligible dry aerosol cores with:

\[
D_{\rm dry}>0.5~\mu{\rm m}.
\]

This is preferable to reconstructing the number from the current hydrometeor-bin centre because the aerosol cores can subsequently experience:

- activation;
- condensation;
- coalescence;
- freezing;
- remapping;
- evaporation;
- sublimation;
- entrainment.

The >0.5 µm core count follows those processes.

BMM also carries a separate extensive `n_demott` provenance moment for currently frozen DeMott-origin primary monomers.

This allows the cumulative DeMott target to account for ice already nucleated without losing provenance after aggregation.

---

# Ice growth and habit

Ice particles carry additional state describing habit and aggregation.

The model contains prognostic moments used to reconstruct quantities such as:

- aspect ratio;
- monomer number;
- effective ice volume;
- effective density;
- rime mass;
- depositional/unrimed mass.

Vapour growth includes:

- capacitance;
- thermal and vapour-transfer corrections;
- ventilation;
- evolving particle shape;
- updated fall speed.

The ice-growth machinery includes Chen-Lamb-based habit development.

Complete ice sublimation can return the associated nonvolatile aerosol and INP residual to the warm/aerosol population.

---

# Collisions and secondary ice

`sce_flag` controls stochastic collection:

```text
0 = no SCE collection
1 = SCE collection
2 = SCE collection with configured secondary-ice processing
```

The collection calculation transports particle number, composition and prognostic moments into the collision product.

Secondary-ice options include:

```fortran
hm_flag
break_flag
mode1_flag
mode2_flag
```

with:

```text
hm_flag       Hallett-Mossop rime splintering
break_flag=1 Vardiman breakup
break_flag=2 Phillips-style collisional breakup
mode1_flag    mode-1 secondary-ice/freezing treatment
mode2_flag    mode-2 secondary-ice treatment
```

Because BMM retains aerosol and ice-property moments through the collision calculation, secondary products can preserve considerably more particle history than is possible in a bulk scheme.

---

# Bin representations

BMM currently supports three principal warm-particle bin treatments through:

```fortran
bin_scheme_flag
```

| Value | Scheme | Description |
|---:|---|---|
| `0` | Fully moving | Representative particle water masses move continuously |
| `1` | Moving centre | Fixed sectional water grid with moving-centre redistribution |
| `2` | Chen-Lamb | Conservative remapping onto a fixed sectional grid |

The alternatives allow numerical sensitivity to sectional representation to be studied directly.

---

## Fully moving scheme

With:

```fortran
bin_scheme_flag = 0
```

particle populations retain moving water-mass coordinates during diffusional growth.

This has the advantage that condensational growth does not itself require repeated projection to fixed water-mass bins.

When additional receiving categories are needed for SCE or population splitting, the model allocates extra full-moving categories.

The full-moving treatment also has explicit handling for aerosol returned by complete hydrometeor evaporation/sublimation.

---

## Hybrid fixed grid for moving-centre and Chen-Lamb

For `bin_scheme_flag=1` and `2`, the default:

```fortran
fixed_grid_mode = 1
```

constructs a hybrid grid.

The first `n_bins` bins are aerosol-adapted:

1. the dry aerosol PSD between `dmina` and `dmaxa` is divided into equal-number intervals;
2. representative aerosol mass uses the correct lognormal third moment;
3. dry boundaries are mapped to equilibrium water-mass boundaries at the initial thermodynamic state.

The model then appends `n_binsc` geometric cloud bins extending to `dmaxc`.

For example:

```text
n_bins  = 60
n_binsc = 80
```

gives:

```text
60 aerosol-adapted bins + 80 cloud/collision bins
```

per external mode.

The cloud-bin ratio is determined from `n_binsc` and `dmaxc`, so the final boundary reaches the requested maximum cloud-drop size.

The historical fixed-grid construction remains available with:

```fortran
fixed_grid_mode = 0
```

for reproducibility.

---

# Aerosol release during evaporation and sublimation

BMM can return nonvolatile aerosol when a liquid or ice hydrometeor disappears completely.

This is controlled by:

```fortran
release_aerosol
```

and is used in inhomogeneous entrainment and chamber-processing pathways.

The returned aerosol carries:

- particle number;
- nonvolatile component masses;
- relevant INP moments;
- DeMott >0.5 µm aerosol-core provenance.

This is important in repeated activation/evaporation cycles because evaporation does not imply destruction of the aerosol particle.

---

## Full-moving residual insertion

For `bin_scheme_flag=0`:

```fortran
full_moving_release_mode
```

selects how released aerosol is returned.

### Mode 0: same numerical population

```text
full_moving_release_mode = 0
```

returns the residual to the same numerical population from which the hydrometeor evaporated or sublimated.

This is the default and is useful when a simple number-conserving response is desired.

### Mode 1: water-coordinate splitting

```text
full_moving_release_mode = 1
```

places the released residual according to its water-mass coordinate.

If the source lies between occupied moving pivots, number and extensive moments are partitioned between the neighbouring populations.

If a new outside-range population is required, for example a dry residual below all currently wet pivots, BMM keeps that population distinct where possible.

Multiple residual sources at the same water coordinate can share one pending receiving population.

If all categories are occupied, a nearby existing pair can be conservatively compacted to free a category rather than immediately averaging a newly dry population into a wet endpoint.

The aim is to retain the physically important distinction between:

```text
wet survivor
```

and:

```text
dry aerosol residual
```

while conserving the carried extensive quantities.

---

# Entrainment and cloud dynamics

BMM can be run as an adiabatic parcel or with entrainment.

`updraft_type` supports several velocity histories, including:

```text
1 = prescribed constant w
2 = prescribed w until t_thresh, then zero
3 = prescribed/oscillatory history
4 = prognostic/model-specific buoyant motion
```

`bubble_flag` selects the geometry used by the entraining dynamics:

```text
.true.  = bubble
.false. = plume / jet
```

Entrainment controls include:

```fortran
adiabatic_prof
entrain_period
vert_ent
ent_rate
thresh_to_start_hom_mix
entrain_aerosol
```

The model can represent both homogeneous and inhomogeneous mixing.

In inhomogeneous mixing, complete evaporation/sublimation can explicitly return aerosol residuals rather than deleting the particles.

Environmental aerosol can also be entrained.

For the fully moving scheme, incoming aerosol is inserted according to its hydrated water state rather than assuming the same source-array index remains the correct receiving population.

---

# Cloud-chamber model

BMM includes a zero-dimensional chamber mode for laboratory expansion-cloud experiments.

The chamber can be driven by measured/prescribed time series of:

- pressure;
- gas temperature;
- total water;
- wall temperature.

The chamber machinery is intentionally modular. It separates:

1. imposed chamber thermodynamic forcing;
2. unresolved wall/boundary-layer processing;
3. homogeneous versus inhomogeneous particle evaporation;
4. wall vapour exchange;
5. fan-associated particle loss;
6. non-gravitational wall deposition;
7. gravitational fallout.

This separation is useful because laboratory chamber measurements constrain some of these processes better than others.

---

## Chamber forcing

The chamber can independently force:

```fortran
chamber_force_pressure
chamber_force_temperature
chamber_force_qtot
```

from arrays in the chamber specification.

Pressure and temperature forcing do not require total-water forcing.

This is important when wall and particle loss processes are treated explicitly: imposing measured total water while also calculating water loss internally can otherwise double count the chamber water sink.

---

## Boundary-layer / recirculation operator

The modern chamber interface uses:

```fortran
chamber_bl_mix
```

as a simple on/off switch.

When enabled, the fraction of chamber air processed over an interval follows:

\[
f_{\rm mix}=1-\exp(-\Delta t/\tau_{\rm BL}),
\]

where:

```fortran
chamber_bl_tau
```

is the effective wall-BL/recirculation timescale.

The nonlinear chamber operator is internally subcycled so that very short `tau` values do not produce one unrealistically large BL event per outer model timestep.

This timescale is best interpreted as an **effective exchange/recirculation time**, not necessarily the residence time of an individual droplet next to a wall.

---

## Effective processed-air temperature

Before latent phase change, the chamber BL operator constructs:

\[
T_{\rm sens}
=
T_{\rm gas}
+
\alpha_T(T_{\rm wall}-T_{\rm gas})
+
\Delta T_{\rm offset}.
\]

The corresponding controls are:

```fortran
chamber_bl_alpha_t
chamber_bl_temp_offset
```

Useful limits are:

```text
alpha_T = 0, offset = 0  -> processed air follows bulk gas temperature
alpha_T = 1, offset = 0  -> processed air follows measured wall temperature
```

A small offset can be used to represent unresolved thermodynamic heterogeneity.

---

## Homogeneous and inhomogeneous chamber evaporation

The wall/BL thermodynamic calculation first determines the amount of liquid water that must evaporate.

Only **after that** does:

```fortran
chamber_bl_evap_mode
```

determine how the same bulk target appears in the PSD.

### Mode 1: homogeneous

```text
chamber_bl_evap_mode = 1
```

retains particle number and decreases liquid mass/size.

### Mode 2: extreme inhomogeneous

```text
chamber_bl_evap_mode = 2
```

completely evaporates selected wet particles to aerosol residuals while surviving particles retain their size.

Every warm population containing liquid can participate; activation status does not gate the chamber evaporation operator.

The homogeneous and inhomogeneous options are required to produce the same bulk liquid-water loss for a given thermodynamic event. They differ only in the resulting particle-size distribution.

---

## Common evaporation size exponent

Both chamber evaporation modes use:

```fortran
chamber_bl_evap_size_exp = p
```

For homogeneous evaporation:

```text
p = 0  equal fractional liquid-mass shrinkage
p = 2  common finite D^2-like decrement
```

For inhomogeneous evaporation, selection is weighted approximately as:

\[
w_i\propto m_{w,i}^{-p/3}.
\]

Thus:

```text
p = 0  uniform number-fraction removal
p = 2  inverse-D^2 lifetime weighting
```

and larger values increasingly favour complete evaporation of smaller wet particles.

This provides a compact way to represent unresolved mixing structure while leaving the bulk liquid-water target unchanged.

---

## Wall-water closures

The wall-vapour treatment is selected with:

```fortran
chamber_bl_wall_water_mode
```

### Mode 0: legacy saturation-cap closure

```text
chamber_bl_wall_water_mode = 0
```

is the historical chamber closure.

It uses the effective processed-air temperature and an instantaneous saturation adjustment. If the processed air becomes supersaturated, vapour is removed according to the legacy saturation-cap treatment; the resulting drying can drive cloud-particle evaporation on remixing.

This mode does **not** use a prognostic finite wall-water reservoir.

It is retained both for reproducibility and because it remains scientifically useful in comparisons with chamber observations.

### Mode 1: finite reservoir with fractional relaxation

```text
chamber_bl_wall_water_mode = 1
```

introduces prognostic liquid/frost reservoirs on the wall.

The exchange toward wall equilibrium is controlled by a fractional relaxation and is coupled to the chamber BL mixing fraction.

Water can be returned from the wall only if it is present in the stored reservoir.

### Mode 2: finite reservoir with physical mass-transfer velocity

```text
chamber_bl_wall_water_mode = 2
```

uses the same finite reservoir but computes an explicit vapour flux from the wall-air vapour-pressure disequilibrium:

\[
J_v
=
k_m
\frac{e_{\rm eq}(T_{\rm wall})-e_{\rm air}}
{R_vT_{\rm gas}}.
\]

The transfer velocity is:

```fortran
chamber_wall_vapour_transfer_velocity
```

and the flux is integrated over the chamber surface area.

Unlike mode 1, the bulk wall-water transfer strength in mode 2 is not set by `chamber_bl_tau`.

The particle evaporation and wall-vapour calculation are coupled so vapour released by cloud evaporation can contribute to the wall deposition flux during the same BL substep.

---

# AIDAd/iSKYLAB interpretation of the chamber treatment

Comparisons with AIDAd/iSKYLAB chamber experiments currently indicate an important result:

> **The AIDAd experiments appear to be reproduced better by the legacy wall-water closure (`chamber_bl_wall_water_mode=0`) than by the newer finite wall-water-reservoir treatments.**

This should not be interpreted simply as evidence that the chamber has no wall-water reservoir. The more useful interpretation is that the laboratory cloud is revealing unresolved **spatial structure and mixing** that a zero-dimensional chamber model cannot represent explicitly.

In the AIDAd cases explored so far, successful simulations commonly require:

- relatively rapid effective chamber mixing/recirculation;
- a small effective thermal perturbation associated with wall-influenced air;
- inhomogeneous evaporation behaviour that changes particle number/PSD structure;
- important loss of water and/or particles to chamber surfaces.

The sign of the fitted effective temperature perturbation is also physically suggestive:

- in warm experiments, where the walls can be colder than the chamber gas, the successful perturbation tends to be negative;
- in colder experiments, where the walls can be warmer than the gas, the successful perturbation can be positive.

The fitted perturbation does not necessarily equal the full measured wall-gas temperature difference. It is better viewed as the temperature signature of **partially mixed wall-conditioned air** in a strongly stirred chamber.

The resulting physical picture is therefore:

```text
bulk chamber air
      |
      v
wall-influenced thermodynamic heterogeneity
      |
      v
rapid fan-driven recirculation / mixing
      |
      +----> homogeneous or inhomogeneous cloud evaporation
      |
      +----> water/particle loss to chamber surfaces
```

The chamber comparisons therefore tell us more than simply how much liquid water is lost.

They provide information about:

1. **cloud structure**  
   The difference between homogeneous and inhomogeneous evaporation constrains whether the observed cloud is behaving like one uniformly mixed population or a mixture of differently processed air.

2. **mixing timescale**  
   The effective `chamber_bl_tau` constrains how rapidly wall-conditioned air is redistributed through the chamber.

3. **surface losses**  
   The inability of a closed-water calculation to reproduce observed LWC shows that losses to the chamber surfaces are important.

4. **particle-number evolution**  
   Fan, wall and fallout losses are separate from reversible droplet evaporation and aerosol release, allowing observed number loss to be distinguished from simple phase cycling.

The finite-reservoir modes remain valuable because they provide a more explicit water budget and are useful for mechanistic sensitivity tests. However, they should not automatically be assumed to be superior simply because they contain more explicit wall physics.

For AIDAd/iSKYLAB, the present evidence suggests that the legacy closure is acting as an effective parameterisation of unresolved wall-conditioned heterogeneity and surface exchange.

This is an empirical conclusion from current chamber comparisons and should be reassessed as additional experiments and diagnostics are added.

---

# Particle loss to surfaces and fallout

The chamber implementation includes three distinct particle-loss pathways.

Keeping these separate is important because they represent different mechanisms and have different size dependence.

---

## Fan-associated particle loss

The optional fan/circulation loss uses a saturating size-dependent first-order rate:

\[
k_{\rm fan}(D)
=
\frac{k_{\max}}
{1+(D_{50}/D)^n}.
\]

Controls include:

```fortran
chamber_fan_loss
chamber_fan_loss_kmax
chamber_fan_loss_d50_ref
chamber_fan_loss_exp
chamber_fan_rpm
chamber_fan_rpm_ref
```

The characteristic diameter is adjusted with fan RPM.

This term should be interpreted as an **effective fan/circulation-associated particle sink**. It can represent losses associated with fan blades and the circulation pattern without requiring every loss event to be literal blade impaction.

The same survival fraction is applied to number and associated extensive moments.

---

## Non-gravitational wall deposition

An independent chamber wall-loss calculation represents non-gravitational deposition to chamber surfaces.

Controls include:

```fortran
chamber_wall_loss
chamber_wall_ustar
chamber_diameter
chamber_height
```

The current implementation uses a smooth-wall deposition formulation based on the Lai-Nazaroff treatment.

This is independent of wall vapour transfer.

---

## Fallout / gravitational sedimentation

With:

```fortran
fallout_flag = .true.
```

BMM applies a residence-time sink based on particle terminal velocity.

For a chamber, the chamber geometry provides the relevant well-mixed vertical scale.

For atmospheric parcel calculations, `residence_depth` supplies the characteristic depth.

The model writes cumulative mass/number fallout and precipitation-flux diagnostics.

---

# Optical properties

BMM can calculate extinction/absorption using:

```fortran
use_adt_optics = .true.
```

for anomalous diffraction theory (ADT).

The optical wavelength is set with:

```fortran
optics_wavelength
```

Refractive-index support is in:

```text
opt/
```

With ADT disabled, the simpler geometric approximation uses an extinction efficiency close to:

```text
Qext = 2
```

with zero absorption in the basic treatment.

---

# Building BMM

## Requirements

The main model requires:

- a Fortran compiler;
- NetCDF C;
- NetCDF-Fortran;
- `make`;
- standard Unix build tools.

The supplied Makefiles use `gfortran` by default.

## NetCDF paths

The top-level Makefile uses:

```make
NETCDF_FOR
NETCDF_C
NETCDF_LIB
```

Set these for your system, for example:

```bash
export NETCDF_FOR=/path/to/netcdf-fortran
export NETCDF_C=/path/to/netcdf-c
```

The normal Fortran library link is:

```make
NETCDF_LIB=-lnetcdff
```

Adjust the local linker flags if your installation also requires explicit NetCDF C linkage.

## Compile

```bash
make -j4
```

The main executable is:

```text
main.exe
```

For development/runtime checking:

```bash
make debug
```

Clean the build with:

```bash
make cleanall
```

or the relevant clean target in the current Makefile.

---

# Running BMM

Run with a namelist path:

```bash
./main.exe namelist.in
```

Alternative examples are under:

```text
configs/
```

for example:

```bash
./main.exe configs/namelist_inp_example.in
```

The NetCDF output filename is specified by:

```fortran
outputfile
```

in the namelist.

The standalone SCE model under `sce/` has its own driver and configuration.

---

# Namelist structure

The core atmospheric/parcel namelist contains groups such as:

```text
&run_vars
&aerosol_setup
&sounding_spec
&aerosol_spec
```

Chamber configurations additionally use chamber-specific groups containing:

- chamber process switches;
- chamber geometry;
- pressure/temperature/water forcing;
- wall temperature;
- wall-BL parameters;
- surface-loss parameters.

The main groups have the following roles.

| Group | Purpose |
|---|---|
| `run_vars` | runtime, timestep, dynamics, ice, bin scheme, SCE, entrainment, fallout and optics |
| `aerosol_setup` | number of modes, bins, components and INP classes |
| `sounding_spec` | environmental sounding/profile |
| `aerosol_spec` | aerosol number, size, composition, hygroscopicity, FHH and INP category |
| chamber options | chamber forcing, BL/wall model, geometry and surface-loss settings |
| chamber specification | time-dependent chamber forcing arrays |
| `cloud_setup` / `cloud_spec` | standalone/fixed SCE cloud-grid configuration |

The root `namelist.in` is a general model configuration.

Additional examples are in:

```text
configs/
```

Experiment-specific configurations are kept with the Python workflows that generate/use them.

---

# NetCDF output

The exact fields depend on enabled physics.

## Parcel/environment

Typical fields include:

```text
time
z
p
t
rh
w
rad_par
```

## Warm particles

Typical warm-phase diagnostics include:

```text
ql
ndrop
deff
nwat
mwat
dwet
maer
```

The sectional fields allow the full wet aerosol/cloud distribution to be reconstructed.

## Ice

With ice enabled, output can include:

```text
qi
nice
nicem
mice
dmaxice
phi
nmon
rhoi
maeri
```

together with additional ice-state moments.

## Grid information

Fixed sectional edges are written as:

```text
mbinedges
```

where relevant.

The output also identifies the bin scheme so analysis tools can apply the correct interpretation.

## Fallout and surface-loss diagnostics

Depending on enabled processes, BMM can write quantities such as:

```text
qfall_liq
qfall_ice
nfall_liq
nfall_ice

qfan_liq
qfan_ice
nfan_liq
nfan_ice

qwall_liq
qwall_ice
nwall_liq
nwall_ice
```

## Chamber wall/BL diagnostics

Current chamber diagnostics include quantities such as:

```text
qchamber_bl
qchamber_bl_evap

qchamber_wall_liq_evap
qchamber_wall_liq_cond
qchamber_wall_ice_subl
qchamber_wall_ice_dep

chamber_wall_liquid_water
chamber_wall_ice_water

chamber_wall_rh
chamber_wall_vapour_flux
```

`qchamber_bl_evap` is a liquid-to-vapour phase-transfer diagnostic and should not automatically be interpreted as net chamber water loss.

The wall-reservoir variables are physical chamber water masses.

---

# Python analysis

The `python/` directory contains general plotting examples and experiment-specific workflows.

Useful general examples include:

```text
python/example_plot_bmm.py
python/example_plot_bmm2.py
python/example_plot_bmm3.py
python/example_PSD_bmm.py
python/example_PSD_bmm2.py
```

Research workflows include directories/scripts for:

```text
AIDA / AIDAd / iSKYLAB
MICC / REFLECT
DCMEX
marine cloud brightening
cloud seeding
cirrus
```

Some workflows require campaign data that are not distributed with BMM.

The AIDAd/iSKYLAB analysis tools can generate experiment-specific chamber namelists, run BMM and compare:

- LWC;
- droplet number;
- ice number;
- effective diameter;
- relative dispersion;
- complete size distributions;
- wall/fan/fallout diagnostics.

---

# Repository layout

```text
.
├── bin_microphysics_module.f90   main microphysics, chamber physics and NetCDF I/O
├── main.f90                      main executable driver
├── namelist.in                   general model configuration
├── configs/                      alternative/example namelists
├── sce/                          stochastic collection / SIP code
├── osnf/                         numerical support routines
├── opt/                          refractive indices and ADT optics
├── python/                       plotting and research workflows
├── batch_runs.sh                 historical/example batch launcher
└── Makefile
```

---

# Development and numerical checking

BMM contains several alternative numerical representations of the same physical processes. Conservation and cross-scheme comparison are therefore important development tests.

Changes affecting any of the following should be checked particularly carefully:

- particle number;
- liquid/ice water mass;
- aerosol component mass;
- INP threshold moments;
- DeMott provenance;
- >0.5 µm aerosol-core number;
- collision-product mapping;
- complete aerosol release;
- homogeneous versus inhomogeneous bulk-LWC consistency.

Use:

```bash
make debug
```

when developing changes that affect indexing, remapping or conservation.

For production/sensitivity calculations, retain:

- the exact BMM revision;
- the namelist used;
- any external chamber/sounding data;
- the NetCDF output;
- analysis configuration.

---

# Notes on interpretation

BMM is deliberately detailed enough that several model variables have a more specific meaning than similarly named bulk quantities.

A few important examples are:

- `ndrop` is the BMM activated-liquid diagnostic, not necessarily identical to number above an instrument diameter threshold;
- `nwat` contains the complete warm population, including aerosol/haze as well as activated droplets;
- `qchamber_bl_evap` is phase transfer from liquid to vapour, not necessarily permanent chamber water loss;
- inhomogeneous evaporation can remove cloud droplets while returning aerosol residuals, so droplet loss is not the same thing as aerosol-number loss;
- fan/wall/fallout terms are permanent particle sinks, whereas complete evaporation is generally a phase/state transition.

These distinctions are especially important when comparing BMM with chamber instruments.

---

# Selected references

BMM incorporates or tests treatments from a broad cloud-microphysics literature. Particularly relevant references include:

- Chen, J.-P. and Lamb, D. (1994), *The Theoretical Basis for the Parameterization of Ice Crystal Habits: Growth by Vapor Deposition*, Journal of the Atmospheric Sciences, 51, 1206–1222.
- Koop, T., Luo, B., Tsias, A. and Peter, T. (2000), *Water activity as the determinant for homogeneous ice nucleation in aqueous solutions*, Nature, 406, 611–614.
- DeMott, P. J. et al. (2010), *Predicting global atmospheric ice nuclei distributions and their impacts on climate*, PNAS, 107, 11217–11222.
- Murray, B. J. et al. (2011), *Heterogeneous freezing of water droplets containing kaolinite particles*, Atmospheric Chemistry and Physics, 11, 4191–4207.
- Niemand, M. et al. (2012), *A Particle-Surface-Area-Based Parameterization of Immersion Freezing on Desert Dust Particles*, Journal of the Atmospheric Sciences, 69, 3077–3092.
- Atkinson, J. D. et al. (2013), *The importance of feldspar for ice nucleation by mineral dust in mixed-phase clouds*, Nature, 498, 355–358.

Further references for secondary-ice and breakup treatments are documented next to the corresponding source routines.

---

# Citation and use

BMM is a research code. If it is used for scientific work, cite the relevant model/application papers and document the configuration used.

For reproducibility, the generated or archived namelist should accompany any reported simulation because the model supports multiple choices for:

- aerosol thermodynamics;
- ice nucleation;
- bin numerics;
- collection/SIP;
- entrainment;
- chamber mixing;
- wall-water exchange;
- particle losses.
