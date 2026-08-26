# base_test_full_run — setup notes and differences

Comparison of this setup directory against the reference GFED setup:

| | Path |
|---|---|
| **A** (reference) | `/projects/0/prjs1474/aarifi/DASEH_R/HAM/T63_echam6-ham_new_v0002_PD_GFED_PE/` |
| **B** (this dir) | `/projects/0/prjs1474/aarifi/PROJECTS/PPE_Project/MM_PPE_ECHAM/setup_templates/mmppe/base_test_full_run/` |

## File inventory

| File | A (GFED_PE) | B (base_test_full_run) |
|---|---|---|
| `emi_spec_*.txt` | yes (119 lines) | yes (126 lines) |
| `symlinks_*.sh` | yes (39302 B) | yes (40199 B) |
| `settings_*` | yes | **missing** |
| `echam_jobscript_*.sh` | no | yes (generated) |
| `readme` | `readme.txt` | `readme.md` (this file) |
| `Optimal_Constraint_Ensembles_1290_non_norm.csv` | yes | no |

The two dirs are at different stages of the toolkit workflow: A holds the *input* settings file
that `jobsubm_echam.sh` sources; B holds the *generated* SLURM runscript that `jobsubm_echam.sh`
emits. They are not directly comparable — see the settings section below.

---

## 1. `emi_spec` — emission specification

**A:** `emi_spec_T63_echam6-ham_new_v0002_PD_GFED_PE.txt`
**B:** `emi_spec_base_test_full_run.txt`

### Summary

Same HAMMOZ template (identical header block; `DUST`, `OCEANI`, `SEASALT`, `VOLCE`, `VOLCC`,
`BIOGENIC` and `TERR` lines unchanged). Three substantive changes: the fire inventory, the CEDS
anthropogenic time handling, and the emission matrix / ALIAS table.

### 1.1 Fire inventory: GFED → GFAS, two sectors → one

| | A | B |
|---|---|---|
| Sectors | `FFIRE` + `GFIRE` | `TFIRE` |
| Source | `GFED5.1ext_daily/%T0/daily/%Y4/emiss_GFED_%C0_wildfire_%Y4_%T0.nc` | `../../aarifi_GFAS_v0006_emissions_inventories_GFAS_CAMS_%T0/daily/%Y4/emiss_GFAS_%C0_wildfire_%Y4_%T0.nc` |
| Path style | relative to `emi_basepath` | relative **escape** out of the versioned tree |

In A, `FFIRE` and `GFIRE` pointed at the *same* GFED wildfire file; the matrix set `FFIRE=1` and
`GFIRE=0`, so `GFIRE` contributed nothing. B collapses this into a single `TFIRE` sector.

Two commented-out `TFIRE` alternatives are kept in B for reference (CMIP6 biomass burning, and the
versioned `gfas/%G0/monthly/...` path).

**Verified:** `emi_basepath` = `${input_basepath}/emissions_inventories`, so `../../` resolves to
`/projects/0/prjs1474/aarifi/INPUT/`, where
`aarifi_GFAS_v0006_emissions_inventories_GFAS_CAMS_T63/daily/` exists and contains
`emiss_GFAS_{bc,c2h6s,oc,so2}_wildfire_<year>_T63.nc`.

**Constraint:** the GFAS daily archive holds **2024 and 2025 only** (plus a `2025_V20250323`
variant). A's GFED source covered a much longer record. Any run outside 2024–2025 will fail to
find `TFIRE` input.

### 1.2 Agricultural waste burning

`AWB` is active in A (`CMIP6biomassburning`, year 2015, `EF_IGNOREYEAR`) and **commented out** in B.

### 1.3 CEDS anthropogenic: transient year → frozen year

Every CEDS sector (`AGR`, `AIRC`, `DOM`, `ENE`, `IND`, `SHIPS`, `SLV`, `TRA`, `WST`) changed from

```
ceds/%T0/%Y4/em-anthropogenic_..._%C0_%Y4_%T0.nc
```

to

```
ceds/%T0/IGNORE_YEAR/em-anthropogenic_..._%C0_2022_2023_IGNORE_YEAR.nc
```

B reads a single pre-merged 2022-06→2023-05 file from an `IGNORE_YEAR` directory instead of
resolving the model year. `EF_IGNOREYEAR` was already set in both, so this hard-freezes the input
*file* rather than just the cycling. See A's `readme.txt` for the rationale (2023 seasonality
looked wrong; data merged from 2022+2023 and the year floored to 2022 by
`/projects/0/prjs1474/aarifi/INPUT/merge_IGNOREYEAR.py`).

**Verified:** `ceds/T63/IGNORE_YEAR/` exists under both `INPUT/v0002` and `INPUT/v0006`.

### 1.4 Emission matrix

Columns `AWB`, `FFIRE`, `GFIRE` removed; `TFIRE` added. Species rows changed as follows:

| Sector | A (BC / DMS / OC / SO2 / SS / DU) | B |
|---|---|---|
| `AGR` | off for all species | **on** for BC, OC, SO2 |
| `SLV` | off for all species | **on** for BC, OC, SO2 |
| `AIRC` | BC only | **BC, OC, SO2** |
| `FFIRE`/`GFIRE` → `TFIRE` | BC, DMS, OC, SO2 (via `FFIRE`) | BC, DMS, OC, SO2 (via `TFIRE`) |

Everything else (`BIOGENIC`, `DOM`, `ENE`, `IND`, `SHIPS`, `TERR`, `TRA`, `WST`, `DUST`, `OCEANI`,
`SEASALT`, `VOLCE`, `VOLCC`) is unchanged.

**Verified:** `em-aircraft_..._OC_2022_2023_IGNORE_YEAR.nc` and the `SO2` equivalent both exist in
`v0002` and `v0006`, so turning `AIRC` on for OC and SO2 will not fail at read time.

### 1.5 ALIAS table

A had one entry, using the uppercase species token found in GFED filenames:

```
DMS   C2H6S   FFIRE GFIRE
```

B has four entries, using the **lowercase** tokens found in GFAS filenames:

```
BC    bc      TFIRE
DMS   c2h6s   TFIRE
OC    oc      TFIRE
SO2   so2     TFIRE
```

This is required — the GFAS files are named `emiss_GFAS_bc_...`, `emiss_GFAS_c2h6s_...`, etc.

---

## 2. `symlinks_*.sh` — input file linking

**A:** `symlinks_T63_echam6-ham_new_v0002_PD_GFED_PE.sh`
**B:** `symlinks_base_test_full_run.sh`

### Summary

Same toolkit script. Five differences: the two base paths, the nudging vertical resolution and
filename convention, the nudging edge-month handling, the ERA5 SST/SIC fallback, and a
faked-restart block present only in B.

### 2.1 Base paths

| | A | B |
|---|---|---|
| `input_basepath` | hardcoded `INPUT/v0002` | `INPUT/${input_files_version}` (from settings) |
| `nudg_basepath` | `process_NDG/monthly_T63L31/era5_T63L31_netcdf` | `process_NDG/monthly_T63L47/` |

B's parameterised `input_basepath` is the more maintainable form, but it **requires a settings file
to define `input_files_version`** — which B currently lacks (section 3).

### 2.2 Nudging: L31 → L47

B points at the **L47** nudging archive, A at **L31**. Any settings file used with B must therefore
set `vres="L47"`, not `"L31"`.

The filename pattern also differs, and each matches its own dataset:

| | A | B |
|---|---|---|
| `nudg_name` | `${prefix}${hres}${vres}_${idate}` | `${prefix}_${hres}${vres}_${idate}` |
| Actual files | `era5T63L31_202404.nc` | `era5_T63L47_202404.nc` |

Note that `monthly_T63L47/era5_T63L47_netcdf/` contains *both* conventions — 32 files with the
underscore and 4 legacy files without. B's underscore form matches the 32-file majority.

**Month coverage:**

| Archive | Coverage |
|---|---|
| L47 (B) | 200605–200607, then **202308–202512** |
| L31 (A) | **202308–202601** |

B's nudging therefore stops one month earlier than A's.

### 2.3 Nudging edge months

ECHAM with nudging needs one month before `exp_start` and one after `stop` for interpolation.

- **A** looks up the *actual* adjacent-month file (with year rollover handling) and falls back to a
  fake link at the boundary month only if that file is absent.
- **B** always creates the fake link, pointing the edge month at the boundary-month file itself.

A's version is the more careful one. B's is the older toolkit default and is still valid — nudging
files contain a few timesteps either side of their calendar month.

### 2.4 ERA5 SST/SIC fallback

A contains a guard that substitutes the 2024 SST/SIC file when a 2026 file is requested but
missing. B does not, and links `era5_sst_${year}_T63.nc` unconditionally.

**Verified:** ERA5 SST/SIC files exist through **2025** in both `v0002` and `v0006` — there is no
2026 file. A run in B whose `stop_year_p1` reaches 2026 will produce a dangling `sst2026` link.

### 2.5 Faked restart (B only)

B adds a `if $faked_restart` block that links restart files from a parent experiment
(hardcoded `parent_exp="PPE_ENS_Control"`), renaming them to the current `exp`. A has no such
block.

B also narrows the pre-run symlink cleanup to preserve `flxatm*`, `sstoce*` and `hdrestart*`:

```bash
find . -type l -and -not -name "${prefix_rerun_file}_${exp}_*${suffix_rerun_file}" \
               -and -not -name "flxatm*" -and -not -name "sstoce*" -and -not -name "hdrestart*" \
     -exec \rm -f {} \;
```

---

## 3. Settings and jobscript

### 3.1 `settings_*` — present in A, missing in B

A has `settings_T63_echam6-ham_new_v0002_PD_GFED_PE`; **B has no settings file at all.**

This is the main gap in B. The settings file is what `jobsubm_echam.sh` sources to define
`exp`, `exp_dir`, `hres`/`vres`, `input_files_version`, `scenario`, the date range, and the full
ECHAM/JSBACH namelists. Without it, `symlinks_base_test_full_run.sh` cannot resolve
`${input_files_version}` (section 2.1) and the run cannot be submitted.

Key values from A, for reference when writing B's settings file:

| Setting | A |
|---|---|
| `model_dir` | `MODEL/aarifi/ECHAM-HAMMOZ-PERTURBATION` |
| `exp` | `ens001` |
| `hres` / `vres` / `oceres` | `T63` / `L31` / `GR15` — **B needs `L47`** |
| `ntiles` | 11 |
| `date_start` / `date_stop` | `2024,04,01` → `2025,08,01` |
| `input_files_version` | `v0002` |
| `scenario` | `ssp245` |
| `queue` / `ncpus` / `walltime` | `genoa` / 192 / `12:00:00` |
| `account` | `srsei10308` |
| `nproma` | 16 |
| `lamip` / `lnudge` / `lmidatm` | `.true.` / `.true.` / `.false.` |
| AeroCom flags | `lHEmon`, `lHEaci`, `lMMPPE`, `lAP3M`, `lAP3D` all `.TRUE.` |
| `&HAMMOZ_PERTURBATIONS` | full PPE parameter block, all at source defaults |

### 3.2 `echam_jobscript_*.sh` — present in B, not in A

B has `echam_jobscript_base_test_full_run.sh`, a **generated** SLURM runscript (the output of
`jobsubm_echam.sh`, not an input). It runs `srun ./echam6`, tars restart files, flips `lresume` to
`.true.`, chains the next job, and calls a CDO post-processing script.

| Setting | B jobscript |
|---|---|
| `--job-name` | `BAse_test_full_run` |
| `--nodes` / `--ntasks` | 1 / 192 |
| `--time` | `14:00:00` |
| `--partition` | `genoa` |
| `--account` | `srs25001` |

**Two things to fix before reusing it:**

1. **Wrong account.** `srs25001` here vs `srsei10308` in A's settings.
2. **Wrong user tree.** Every path points into another user's directories —
   `/home/ybhatti2/prjs1474/Pace_PPE_Output/...` for output, `.err`/`.out` and post-processing, and
   `/gpfs/work3/0/prjs1474/ybhatti/PPE_Project/...` for the job-chaining `cd` and `sbatch`. This
   file was copied from a `ybhatti` MM_PPE experiment and has not been re-pointed at `aarifi`.

Since this file is regenerated by `jobsubm_echam.sh` from the settings file, the right fix is to
write B's settings file (section 3.1) and let the toolkit regenerate the jobscript, rather than
editing it by hand.

---

## 4. Net effect

**A (`T63_echam6-ham_new_v0002_PD_GFED_PE`)** — transient GFED fire emissions, narrower
anthropogenic sector set, L31 nudging, complete and runnable.

**B (`base_test_full_run`)** — fixed-year GFAS fire emissions, broader anthropogenic sector set
(`AGR`, `SLV`, and aircraft OC/SO2 additionally on), L47 nudging, faked-restart support. Not yet
runnable: no settings file, and the generated jobscript points at another user's tree.

## 5. Open items

- [ ] Write `settings_base_test_full_run` (must set `vres="L47"` and `input_files_version`).
- [ ] Regenerate the jobscript so account and output paths point at `aarifi`, not `ybhatti2`.
- [ ] Confirm the run window stays inside **2024–2025** (GFAS daily) and **202308–202512** (L47 nudging).
- [ ] Decide whether dropping `AWB` is intended, given `AGR` was turned on.
- [ ] Consider porting A's ERA5 2026 fallback and adjacent-month nudging lookup into B's symlinks script.

## 6. Context from A's `readme.txt`

Reproduced here since it documents the input preparation both setups depend on:

> The input files here used are only applicable for 2024 to 2025.
> We applied following changes to increase the accuracy and credibility.
> 1. New GFAS data for 2025 has been produced, check `/projects/0/prjs1474/aarifi/INPUT/process_GFAS`
> 2. For the ceds data we use the time period 2022-06 to 2023-05 and it's sourced from v0006. The
>    reason for not using the entire year 2023 is that the seasonality seems to be not correct. We
>    merged the data together from 2022 and 2023 and floor the year from 2023 to 2022 after that.
>    Final data in `/projects/0/prjs1474/aarifi/INPUT/v0006/emissions_inventories/ceds/T63/IGNORE_YEAR`,
>    processed by `/projects/0/prjs1474/aarifi/INPUT/merge_IGNOREYEAR.py`
> 3. Use ssp scenarios instead of rcp scenarios for GHG and ozone, copy the data over from v0006.
>    Note that before, `lamip` was false, using it as a climatology.
>
> Ask Hailing if this method is suitable or not
