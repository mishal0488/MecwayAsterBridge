# Mecway-route: vergelijking met de gewone run

Gemaakt door `run_verificatie.py --mecway-opslaan` op 2026-10-07. Elk model is eerst door Mecway geopend en opgeslagen (`TestFiles/Tools/mw_save_as.py`), daarna door de bridge en Code-Aster gerekend en vergeleken met `../resultaten.json`.

**1611 controles identiek** (zelfde status, FEM-waarde gelijk tot 1e-6 relatief, of tot 1e-6 x de absolute tolerantie bij grootheden rond 0), **0 anders**, **0 niet mogelijk** (270 testcases).

.comm (zonder commentaarregels en paden) anders dan in de gewone run: geen; alleen in de laatste cijfers van getallen (Mecway schrijft 15 significante cijfers): DY-01a, DY-01b, DY-01c, DY-01d, DY-02, HR-01d, HR-02a, HR-02b, HR-02c, MC-01b, MC-02b, MD-03, MP-01, MP-02, MP-03.

