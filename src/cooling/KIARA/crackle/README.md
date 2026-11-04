# Crackle in KIARA

This directory holds the crackle cooling and chemistry solver, together with
the parts of Grackle it needs (rate tables, Cloudy metal-cooling tables, UV
background, units and parameters). It is compiled into SWIFT for the KIARA
cooling and KIARA-RT modules, so KIARA no longer links an external Grackle.

Imported on 2026-10-29 from:

- crackle: https://github.com/romeeld/crackle, commit 6937733
  ("Do not attenuate dust ISRF heating by H2 self-shielding; double sol").
  That repository is frozen; crackle is now developed here.
- Grackle: https://github.com/EwanBJones98/grackle, branch `kiara`, commit
  908cc83 (the C part of src/clib only; none of the Fortran is used by
  crackle).

Licences: LICENSE.grackle (Grackle, BSD), LICENSE.crackle.
