# Code-Aster-validatietests (installatie)

Gemaakt door `run_code_aster_tests.py` op 2026-10-06 met `%LOCALAPPDATA%\code_aster\v2025\17.4\bin\run_aster.bat`. Elke test draait als `run_aster.bat --ctest <test>.export` (zoals ctest) in een tijdelijke map. Toelichting: `docs/Verificatie.md` (hoofdstuk Code-Aster-tests).

**19 van 19 tests OK**; 1598 TEST_RESU-controles, waarvan 0 NOOK.

| Test | Bridge-onderdeel | Onderwerp | Wat geverifieerd wordt | Referentie | Controles (type: aantal) | Max. fout t.o.v. referentie | Uitkomst | Tijd [s] |
|---|---|---|---|---|---|---|---|---|
| sslv04b | MECA_STATIQUE 3D (hex20/penta15) | Dikwandige holle cilinder onder inwendige druk 60 MPa | radiale verplaatsing en spanningen (Lame), 3D kwadratisch, PRES_REP | analytisch (Lame) | ANALYTIQUE 77, NON_REGRESSION 86 | ANALYTIQUE 4.8% (tol. 6%) | OK | 3 |
| ssll102a | MECA_STATIQUE POU_D_E | Ingeklemde balk met eenheidskrachten en -momenten op het vrije eind | verplaatsingen, rotaties en snedekrachten van Euler-Bernoulli-balken | analytisch (balktheorie) | ANALYTIQUE 318, NON_REGRESSION 286 | ANALYTIQUE 8.6e-05% (tol. 0.0001%) | OK | 7 |
| ssls101b | MECA_STATIQUE DKT | Ronde plaat (kwart model) onder druk, F3 en eigen gewicht | doorbuiging en momenten in DKT-schaalelementen; FORCE_COQUE, PESANTEUR | analytisch (plaattheorie) | ANALYTIQUE 35, NON_DEFINI 4, NON_REGRESSION 36 | ANALYTIQUE 5.6% (tol. 6%), NON_DEFINI 0.12% (tol. 0.15%) | OK | 3 |
| ssls105a | MECA_STATIQUE DKT (dubbel gekromd) | NAFEMS LE3: halve bol, dubbel geknepen (Morley) | radiale verplaatsing in het belastingspunt; grote rotaties in het vlak | NAFEMS LE3 (SOURCE_EXTERNE) | SOURCE_EXTERNE 2, NON_REGRESSION 2 | SOURCE_EXTERNE 0.66% (tol. 1%) | OK | 3 |
| ssls114a | MECA_STATIQUE COQUE_3D | Cilinderschaal onder inwendige druk (PRES_REP en FORCE_COQUE) | verplaatsingen, membraankrachten en rekken in COQUE_3D | analytisch (membraantheorie) | ANALYTIQUE 48, NON_REGRESSION 48 | ANALYTIQUE 27% (tol. 28%) | OK | 3 |
| ssnp15a | STAT_NON_LINE VMIS_ISOT_LINE | Vierkante plaat, niet-proportionele trek + afschuiving, 3D | spanningen en plastische rek bij isotrope versteviging (Von Mises) | analytisch (integratie van de constitutieve wet) | ANALYTIQUE 44, AUTRE_ASTER 34, NON_REGRESSION 94 | ANALYTIQUE 2.4% (tol. 3%), AUTRE_ASTER 0.064% (tol. 0.1%) | OK | 8 |
| ssnp14b | STAT_NON_LINE VMIS_CINE_LINE | NAFEMS: vierkante plaat, trek + afschuiving, kinematische versteviging | spanningen en plastische rek bij kinematische versteviging, penta6 | NAFEMS (in de .comm: NON_DEFINI) | NON_REGRESSION 13 |  | OK | 3 |
| hsnv100g | STAT_NON_LINE VMIS_ISOT_TRAC + thermische rek | Thermoplasticiteit: element met opgelegde temperatuur en trekkromme | spanning en plastische rek; VMIS_ISOT_TRAC, ALPHA, AFFE_VARC TEMP | analytisch | ANALYTIQUE 12, NON_REGRESSION 12 | ANALYTIQUE 0.001% (tol. 0.1%) | OK | 3 |
| ssnp159a | STAT_NON_LINE GDEF_LOG | Blok in trek tot grote rek, VMIS_ISOT_LINE, PETIT / GDEF_LOG / SIMO_MIEHE | spanning en rek bij grote vervorming (logaritmische rek) | analytisch | ANALYTIQUE 6, NON_DEFINI 9, NON_REGRESSION 15 | ANALYTIQUE 0.15% (tol. 1%), NON_DEFINI 1% (tol. 1.5%) | OK | 4 |
| ssnv122a | STAT_NON_LINE GREEN_LAGRANGE, meebewegende druk | Starre rotatie van een blok + druk in de vervormde configuratie | objectiviteit van GREEN_LAGRANGE; PRES_REP met TYPE_CHARGE='SUIV' | analytisch | ANALYTIQUE 18, NON_REGRESSION 12 | ANALYTIQUE 0.004% (tol. 0.1%) | OK | 3 |
| ssnv167d | contact DEFI_CONTACT CONTINUE | Twee blokken in contact, opgelegde verplaatsing, 3D | contactdruk (LAGS_C) en verplaatsingen; continue formulering, NEWTON | analytisch | ANALYTIQUE 20, NON_REGRESSION 20 | ANALYTIQUE 4.2e-14% (tol. 1e-06%) | OK | 2 |
| sdls115a | DYNA_NON_LINE DIFF_CENT en HHT | Plaatstrook (DKT) met sinusvormige oplossing u = q0 sin(kx) sin(wt) | verplaatsing in de tijd met expliciet (DIFF_CENT) en HHT compleet | analytisch | ANALYTIQUE 8, AUTRE_ASTER 3, NON_REGRESSION 12 | ANALYTIQUE 5.5% (tol. 10%), AUTRE_ASTER 0.00069% (tol. 0.001%) | OK | 4 |
| sdlv121d | DYNA_NON_LINE HHT / NEWMARK (continuum) | Vlakke golf in een elastische kolom: voortplanting, reflectie, absorptie | verplaatsing/snelheid van de golf, D_PLAN, HHT en NEWMARK | analytisch (golfvoortplanting) | ANALYTIQUE 8, NON_REGRESSION 8 | ANALYTIQUE 2.1% (tol. 2.2%) | OK | 9 |
| tpll01e | THER_LINEAIRE stationair (3D) | Wand met opgelegde temperatuur, flux en convectie (1D-oplossing in 3D) | temperatuur; TEMP_IMPO, FLUX_REP, ECHANGE | analytisch | ANALYTIQUE 7, NON_REGRESSION 7 | ANALYTIQUE 4.3e-07% (tol. 0.1%) | OK | 3 |
| ttlv01a | THER_LINEAIRE transient (3D) | Afkoelend/opwarmend lichaam met convectie (penta6/tetra4) | temperatuur in de tijd; ECHANGE, THETA-schema, ETAT_INIT | analytisch (reeksoplossing) | ANALYTIQUE 50, NON_REGRESSION 50 | ANALYTIQUE 3.3% (tol. 5%) | OK | 5 |
| tpna01b | THER_NON_LINE straling (3D) | Stationaire geleiding met straling naar de omgeving | temperatuur en flux; RAYONNEMENT (T+273,15 intern), THER_NL | analytisch | ANALYTIQUE 5, NON_REGRESSION 7 | ANALYTIQUE 0.087% (tol. 0.2%) | OK | 4 |
| tplv101c | THER_NON_LINE PLAN (thermische schaal in de bridge) | Vlak model met bronnen, flux, opgelegde temperatuur en wandcontact | temperatuur; MODELISATION='PLAN', SOURCE, FLUX_REP, TEMP_IMPO | analytisch | ANALYTIQUE 70, NON_REGRESSION 70 | ANALYTIQUE 1.3e-12% (tol. 0.1%) | OK | 4 |
| sslv153a | CALC_ERREUR ERME (foutschatter) | Blok van hexa/penta/pyram/tetra (LINE_QUAD), druk, vrije randen met druk ~0 | ERREST van ERME_NOEU en QIRE_NOEU; ook via de bridge nagebouwd (FS-03) | NON_REGRESSION | NON_REGRESSION 8 |  | OK | 2 |
| sslv111f | CALC_ERREUR ERME en ERZ1 (3D) | Plaat met gat (Kirsch), hexa8/pyram5/tetra4, FORCE_FACE | ERME_ELEM/ERME_ELNO en ZZ1-spanningen in 3D | analytisch (Kirsch) / NON_REGRESSION | ANALYTIQUE 17, NON_REGRESSION 17 | ANALYTIQUE 12% (tol. 15%) | OK | 9 |

## Gebruikte Code-Aster-opties (uit de `.code`)

| Test | Opties |
|---|---|
| sslv04b | mod 3D; solver MUMPS |
| ssll102a | mod POU_D_E; solver MUMPS |
| ssls101b | mod DKT; solver MUMPS |
| ssls105a | mod DKT |
| ssls114a | mod COQUE_3D; solver MUMPS |
| ssnp15a | mod 3D; rel VMIS_ECMI_LINE/VMIS_ISOT_LINE; def PETIT; solver MUMPS |
| ssnp14b | mod 3D; rel VMIS_CINE_LINE; def PETIT; solver MUMPS |
| hsnv100g | mod 3D; rel VMIS_ISOT_TRAC; def PETIT; schema THETA; solver MUMPS |
| ssnp159a | mod 3D; rel VMIS_ISOT_LINE; def GDEF_LOG/PETIT/SIMO_MIEHE; solver MUMPS |
| ssnv122a | mod 3D; rel ELAS; def GREEN_LAGRANGE; solver MUMPS |
| ssnv167d | mod 3D; rel ELAS; def PETIT; form CONTINUE; solver MUMPS |
| sdls115a | mod DKT; rel ELAS; def PETIT; schema DIFF_CENT/HHT; form ACCELERATION/DEPLACEMENT; solver MUMPS |
| sdlv121d | mod D_PLAN/D_PLAN_ABSO; rel ELAS; def PETIT; schema HHT/NEWMARK; form DEPLACEMENT; solver LDLT/MUMPS |
| tpll01e | mod 3D; schema THETA; solver LDLT |
| ttlv01a | mod 3D; schema THETA; solver MUMPS |
| tpna01b | mod 3D; rel THER_NL; solver MUMPS |
| tplv101c | mod PLAN; rel THER_NL; schema THETA; solver MUMPS |
| sslv153a | mod 3D; solver MUMPS |
| sslv111f | mod 3D; solver MUMPS |

## Alarmen (diagnose <A>_ALARM, test geslaagd)

| Test | Alarmen |
|---|---|
| ssls114a | MODELISA8_14, UTILITAI_4 |
| ssnv122a | CHARGES9_53 |
| sdls115a | CATAMESS_41, MODELISA8_14 |
| sdlv121d | MECANONLINE_23 |
| sslv153a | CATAMESS_41, CHARGES2_88 |
