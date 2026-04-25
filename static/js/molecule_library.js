/**
 * molecule_library.js — shared curated SMILES catalog.
 *
 * Source of truth for "Pick from library" pickers across every page.
 * Each entry is { name, smiles, [synonyms] }. Synonyms broaden the
 * search index without cluttering the visible label.
 *
 * Coverage targets:
 *   - Common FDA drugs (top-prescribed + landmark molecules)
 *   - Solvents / reagents bench chemists actually reach for
 *   - Heterocyclic & carbocyclic scaffolds
 *   - All 20 proteinogenic amino acids
 *   - Carbohydrates / sugars
 *   - Fatty acids / lipids
 *   - Vitamins
 *   - Hormones & neurotransmitters
 *   - Natural products (terpenes / alkaloids / flavonoids)
 *   - Polymer monomers
 *   - Pesticides / herbicides
 *   - Astex / Rule-of-3 fragment exemplars
 *   - Stereochemistry / functional-group test molecules
 *
 * Public API (window.NuGenLibrary):
 *   list()                         → flat [{category, name, smiles, ...}]
 *   categories()                   → ["Common drugs", ...]
 *   byCategory()                   → { category: [entry, ...] }
 *   search(q, opts)                → ranked entries (fuzzy)
 *   count()                        → integer
 */
(function () {
  "use strict";

  const CATALOG = {
    "Common drugs": [
      { name: "Aspirin",          smiles: "CC(=O)OC1=CC=CC=C1C(=O)O", synonyms: "acetylsalicylic NSAID analgesic" },
      { name: "Caffeine",         smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", synonyms: "stimulant xanthine" },
      { name: "Ibuprofen",        smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", synonyms: "NSAID analgesic advil" },
      { name: "Paracetamol",      smiles: "CC(=O)NC1=CC=C(O)C=C1", synonyms: "acetaminophen tylenol" },
      { name: "Naproxen",         smiles: "COc1ccc2cc(C(C)C(=O)O)ccc2c1", synonyms: "NSAID aleve" },
      { name: "Diclofenac",       smiles: "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl", synonyms: "NSAID voltaren" },
      { name: "Tamoxifen",        smiles: "CCC(=C(C1=CC=CC=C1)C2=CC=C(C=C2)OCCN(C)C)C3=CC=CC=C3", synonyms: "SERM estrogen breast cancer" },
      { name: "Penicillin G",     smiles: "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O", synonyms: "antibiotic beta-lactam" },
      { name: "Sildenafil",       smiles: "CCCc1nn(C)c2c1NC(c1cc(S(=O)(=O)N3CCN(C)CC3)ccc1OCC)=NC2=O", synonyms: "viagra ED PDE5" },
      { name: "Atorvastatin",     smiles: "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O", synonyms: "lipitor statin cholesterol" },
      { name: "Metformin",        smiles: "CN(C)C(=N)NC(=N)N", synonyms: "diabetes biguanide" },
      { name: "Omeprazole",       smiles: "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1", synonyms: "PPI prilosec heartburn" },
      { name: "Amoxicillin",      smiles: "CC1(C)S[C@@H]2[C@H](NC(=O)[C@@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O", synonyms: "antibiotic beta-lactam" },
      { name: "Ciprofloxacin",    smiles: "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O", synonyms: "quinolone antibiotic" },
      { name: "Lisinopril",       smiles: "NCCCC[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O", synonyms: "ACE inhibitor blood pressure" },
      { name: "Atenolol",         smiles: "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1", synonyms: "beta-blocker hypertension" },
      { name: "Levothyroxine",    smiles: "Oc1cc(I)c(Oc2cc(I)c(C[C@H](N)C(=O)O)cc2I)c(I)c1", synonyms: "T4 thyroid synthroid" },
      { name: "Warfarin",         smiles: "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O", synonyms: "anticoagulant blood thinner" },
      { name: "Diazepam",         smiles: "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21", synonyms: "valium benzodiazepine anxiety" },
      { name: "Fluoxetine",       smiles: "CNCCC(c1ccccc1)Oc1ccc(C(F)(F)F)cc1", synonyms: "prozac SSRI antidepressant" },
      { name: "Sertraline",       smiles: "CN[C@H]1CC[C@@H](c2ccc(Cl)c(Cl)c2)c2ccccc21", synonyms: "zoloft SSRI antidepressant" },
      { name: "Loratadine",       smiles: "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc32)CC1", synonyms: "claritin antihistamine" },
      { name: "Albuterol",        smiles: "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1", synonyms: "salbutamol asthma bronchodilator" },
      { name: "Morphine",         smiles: "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", synonyms: "opioid analgesic" },
      { name: "Codeine",          smiles: "COc1ccc2C[C@H]3N(C)CC[C@@]45[C@H](Oc1c24)[C@H](O)C=C[C@@H]35", synonyms: "opioid antitussive" },
      { name: "Imatinib",         smiles: "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1", synonyms: "gleevec kinase inhibitor leukemia" },
      { name: "Sofosbuvir",       smiles: "CC(C)OC(=O)[C@@H](C)NP(=O)(OC[C@H]1O[C@@H](n2ccc(=O)[nH]c2=O)[C@](C)(F)[C@@H]1O)Oc1ccccc1", synonyms: "HCV antiviral" },
      { name: "Remdesivir",       smiles: "CCC(CC)COC(=O)[C@H](C)NP(=O)(OC[C@H]1O[C@](C#N)([C@H](O)[C@@H]1O)c1ccc2c(N)ncnn12)Oc1ccccc1", synonyms: "antiviral covid" },
    ],

    "Solvents & reagents": [
      { name: "Water",            smiles: "O" },
      { name: "Methanol",         smiles: "CO" },
      { name: "Ethanol",          smiles: "CCO", synonyms: "EtOH" },
      { name: "Isopropanol",      smiles: "CC(C)O", synonyms: "IPA" },
      { name: "n-Butanol",        smiles: "CCCCO" },
      { name: "tert-Butanol",     smiles: "CC(C)(C)O", synonyms: "tBuOH" },
      { name: "Acetone",          smiles: "CC(=O)C" },
      { name: "Acetonitrile",     smiles: "CC#N", synonyms: "MeCN ACN" },
      { name: "DMSO",             smiles: "CS(=O)C", synonyms: "dimethyl sulfoxide" },
      { name: "DMF",              smiles: "CN(C)C=O", synonyms: "dimethylformamide" },
      { name: "DMAc",             smiles: "CC(=O)N(C)C", synonyms: "dimethylacetamide" },
      { name: "NMP",              smiles: "CN1CCCC1=O", synonyms: "N-methyl pyrrolidone" },
      { name: "THF",              smiles: "C1CCOC1", synonyms: "tetrahydrofuran" },
      { name: "1,4-Dioxane",      smiles: "C1COCCO1" },
      { name: "Diethyl ether",    smiles: "CCOCC", synonyms: "Et2O" },
      { name: "Acetic acid",      smiles: "CC(=O)O", synonyms: "AcOH" },
      { name: "Formic acid",      smiles: "OC=O" },
      { name: "Trifluoroacetic acid", smiles: "OC(=O)C(F)(F)F", synonyms: "TFA" },
      { name: "DCM",              smiles: "ClCCl", synonyms: "dichloromethane CH2Cl2" },
      { name: "Chloroform",       smiles: "ClC(Cl)Cl", synonyms: "CHCl3" },
      { name: "Carbon tet",       smiles: "ClC(Cl)(Cl)Cl", synonyms: "CCl4" },
      { name: "Toluene",          smiles: "Cc1ccccc1", synonyms: "PhMe" },
      { name: "Benzene",          smiles: "c1ccccc1" },
      { name: "Xylene (p-)",      smiles: "Cc1ccc(C)cc1" },
      { name: "Hexane",           smiles: "CCCCCC" },
      { name: "Cyclohexane",      smiles: "C1CCCCC1" },
      { name: "Pentane",          smiles: "CCCCC" },
      { name: "Ethyl acetate",    smiles: "CCOC(C)=O", synonyms: "EtOAc" },
      { name: "Pyridine",         smiles: "c1ccncc1" },
      { name: "Triethylamine",    smiles: "CCN(CC)CC", synonyms: "TEA Et3N" },
      { name: "Diisopropylethylamine", smiles: "CC(C)N(CC)C(C)C", synonyms: "DIPEA Hünig" },
    ],

    "Scaffolds & rings": [
      { name: "Benzene",          smiles: "c1ccccc1" },
      { name: "Cyclopropane",     smiles: "C1CC1" },
      { name: "Cyclobutane",      smiles: "C1CCC1" },
      { name: "Cyclopentane",     smiles: "C1CCCC1" },
      { name: "Cyclohexane",      smiles: "C1CCCCC1" },
      { name: "Cycloheptane",     smiles: "C1CCCCCC1" },
      { name: "Adamantane",       smiles: "C1C2CC3CC1CC(C2)C3" },
      { name: "Norbornane",       smiles: "C1CC2CCC1C2", synonyms: "bicyclic" },
      { name: "Naphthalene",      smiles: "c1ccc2ccccc2c1" },
      { name: "Anthracene",       smiles: "c1ccc2cc3ccccc3cc2c1" },
      { name: "Phenanthrene",     smiles: "c1ccc2c(c1)ccc1ccccc12" },
      { name: "Biphenyl",         smiles: "c1ccc(-c2ccccc2)cc1" },
      { name: "Pyridine",         smiles: "c1ccncc1" },
      { name: "Pyrazine",         smiles: "c1cnccn1" },
      { name: "Pyrimidine",       smiles: "c1cncnc1" },
      { name: "Pyridazine",       smiles: "c1ccnnc1" },
      { name: "Triazine",         smiles: "c1cnncn1" },
      { name: "Pyrrole",          smiles: "c1cc[nH]c1" },
      { name: "Furan",            smiles: "c1ccoc1" },
      { name: "Thiophene",        smiles: "c1ccsc1" },
      { name: "Imidazole",        smiles: "c1nc[nH]c1" },
      { name: "Oxazole",          smiles: "c1ocnc1" },
      { name: "Thiazole",         smiles: "c1csnc1" },
      { name: "Pyrazole",         smiles: "c1cn[nH]c1" },
      { name: "Triazole (1,2,3)", smiles: "c1cn[nH]n1" },
      { name: "Tetrazole",        smiles: "c1nnn[nH]1" },
      { name: "Indole",           smiles: "c1ccc2[nH]ccc2c1" },
      { name: "Benzofuran",       smiles: "c1ccc2occc2c1" },
      { name: "Benzothiophene",   smiles: "c1ccc2sccc2c1" },
      { name: "Benzimidazole",    smiles: "c1ccc2[nH]cnc2c1" },
      { name: "Quinoline",        smiles: "c1ccc2ncccc2c1" },
      { name: "Isoquinoline",     smiles: "c1ccc2cnccc2c1" },
      { name: "Quinazoline",      smiles: "c1ccc2ncncc2c1" },
      { name: "Purine",           smiles: "c1[nH]cnc2ncncc12" },
      { name: "Indazole",         smiles: "c1ccc2[nH]ncc2c1" },
      { name: "Carbazole",        smiles: "c1ccc2c(c1)[nH]c1ccccc12" },
      { name: "Piperidine",       smiles: "C1CCNCC1" },
      { name: "Piperazine",       smiles: "C1CNCCN1" },
      { name: "Morpholine",       smiles: "O1CCNCC1" },
      { name: "Pyrrolidine",      smiles: "C1CCNC1" },
      { name: "Tetrahydrofuran",  smiles: "C1CCOC1", synonyms: "THF saturated" },
      { name: "Tetrahydropyran",  smiles: "C1CCOCC1" },
    ],

    "Amino acids": [
      { name: "Glycine (Gly)",       smiles: "NCC(=O)O" },
      { name: "Alanine (Ala)",       smiles: "C[C@@H](N)C(=O)O" },
      { name: "Valine (Val)",        smiles: "CC(C)[C@@H](N)C(=O)O" },
      { name: "Leucine (Leu)",       smiles: "CC(C)C[C@@H](N)C(=O)O" },
      { name: "Isoleucine (Ile)",    smiles: "CC[C@H](C)[C@@H](N)C(=O)O" },
      { name: "Proline (Pro)",       smiles: "OC(=O)[C@@H]1CCCN1" },
      { name: "Phenylalanine (Phe)", smiles: "N[C@@H](Cc1ccccc1)C(=O)O" },
      { name: "Tryptophan (Trp)",    smiles: "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O" },
      { name: "Methionine (Met)",    smiles: "CSCC[C@@H](N)C(=O)O" },
      { name: "Serine (Ser)",        smiles: "N[C@@H](CO)C(=O)O" },
      { name: "Threonine (Thr)",     smiles: "C[C@@H](O)[C@@H](N)C(=O)O" },
      { name: "Cysteine (Cys)",      smiles: "N[C@@H](CS)C(=O)O" },
      { name: "Tyrosine (Tyr)",      smiles: "N[C@@H](Cc1ccc(O)cc1)C(=O)O" },
      { name: "Asparagine (Asn)",    smiles: "NC(=O)C[C@@H](N)C(=O)O" },
      { name: "Glutamine (Gln)",     smiles: "NC(=O)CC[C@@H](N)C(=O)O" },
      { name: "Lysine (Lys)",        smiles: "NCCCC[C@@H](N)C(=O)O" },
      { name: "Arginine (Arg)",      smiles: "N=C(N)NCCC[C@@H](N)C(=O)O" },
      { name: "Histidine (His)",     smiles: "N[C@@H](Cc1cnc[nH]1)C(=O)O" },
      { name: "Aspartate (Asp)",     smiles: "N[C@@H](CC(=O)O)C(=O)O" },
      { name: "Glutamate (Glu)",     smiles: "N[C@@H](CCC(=O)O)C(=O)O" },
    ],

    "Sugars & carbohydrates": [
      { name: "Glucose (D)",       smiles: "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O" },
      { name: "Fructose (D)",      smiles: "OCC(=O)[C@@H](O)[C@H](O)[C@H](O)CO" },
      { name: "Galactose (D)",     smiles: "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O" },
      { name: "Mannose (D)",       smiles: "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O" },
      { name: "Ribose (D)",        smiles: "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O" },
      { name: "Sucrose",           smiles: "OC[C@H]1O[C@H](O[C@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", synonyms: "table sugar" },
      { name: "Lactose",           smiles: "OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@H]1O", synonyms: "milk sugar" },
      { name: "Maltose",           smiles: "OC[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O" },
      { name: "Cellobiose",        smiles: "OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O" },
      { name: "Glycerol",          smiles: "OCC(O)CO", synonyms: "glycerin" },
      { name: "Xylitol",           smiles: "OC[C@@H](O)C(O)[C@H](O)CO" },
      { name: "Sorbitol",          smiles: "OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO" },
    ],

    "Lipids & fatty acids": [
      { name: "Acetic acid",       smiles: "CC(=O)O" },
      { name: "Butyric acid",      smiles: "CCCC(=O)O" },
      { name: "Hexanoic acid",     smiles: "CCCCCC(=O)O", synonyms: "caproic" },
      { name: "Octanoic acid",     smiles: "CCCCCCCC(=O)O", synonyms: "caprylic" },
      { name: "Lauric acid",       smiles: "CCCCCCCCCCCC(=O)O", synonyms: "C12:0" },
      { name: "Myristic acid",     smiles: "CCCCCCCCCCCCCC(=O)O", synonyms: "C14:0" },
      { name: "Palmitic acid",     smiles: "CCCCCCCCCCCCCCCC(=O)O", synonyms: "C16:0 saturated" },
      { name: "Stearic acid",      smiles: "CCCCCCCCCCCCCCCCCC(=O)O", synonyms: "C18:0 saturated" },
      { name: "Oleic acid",        smiles: "CCCCCCCC/C=C\\CCCCCCCC(=O)O", synonyms: "C18:1 omega-9" },
      { name: "Linoleic acid",     smiles: "CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)O", synonyms: "C18:2 omega-6" },
      { name: "Alpha-linolenic",   smiles: "CC/C=C\\C/C=C\\C/C=C\\CCCCCCCC(=O)O", synonyms: "ALA C18:3 omega-3" },
      { name: "Arachidonic acid",  smiles: "CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(=O)O", synonyms: "AA C20:4" },
      { name: "EPA",               smiles: "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(=O)O", synonyms: "eicosapentaenoic omega-3" },
      { name: "DHA",               smiles: "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)O", synonyms: "docosahexaenoic omega-3" },
      { name: "Cholesterol",       smiles: "CC(C)CCC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(C)CC[C@H]1[C@H]2CC=C2C[C@@H](O)CC[C@]12C" },
    ],

    "Vitamins": [
      { name: "Vitamin C",         smiles: "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O", synonyms: "ascorbic acid" },
      { name: "Vitamin B1",        smiles: "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1", synonyms: "thiamine" },
      { name: "Vitamin B2",        smiles: "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)[C@H](O)CO)c2cc1C", synonyms: "riboflavin" },
      { name: "Vitamin B3",        smiles: "OC(=O)c1cccnc1", synonyms: "niacin nicotinic acid" },
      { name: "Vitamin B6",        smiles: "Cc1ncc(CO)c(CO)c1O", synonyms: "pyridoxine" },
      { name: "Vitamin B12",       smiles: "[Co]", synonyms: "cobalamin (placeholder)" },
      { name: "Folic acid",        smiles: "Nc1nc2ncc(CNc3ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc3)nc2c(=O)[nH]1", synonyms: "vitamin B9" },
      { name: "Biotin",            smiles: "O=C1N[C@@H]2CS[C@@H](CCCCC(=O)O)[C@@H]2N1", synonyms: "vitamin B7" },
      { name: "Vitamin A",         smiles: "OC/C=C(\\C)/C=C/C=C(\\C)/C=C/C1=C(C)CCCC1(C)C", synonyms: "retinol" },
      { name: "Vitamin D3",        smiles: "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2/C(=C\\C=C3\\C[C@@H](O)CCC3=C)CCC[C@]12C", synonyms: "cholecalciferol" },
      { name: "Vitamin E",         smiles: "Cc1c(C)c2c(c(C)c1O)CC[C@](C)(CCC[C@H](C)CCC[C@H](C)CCCC(C)C)O2", synonyms: "alpha-tocopherol" },
      { name: "Vitamin K1",        smiles: "CC1=C(C/C=C(\\C)CCCC(C)CCCC(C)CCCC(C)C)C(=O)c2ccccc2C1=O", synonyms: "phylloquinone" },
    ],

    "Hormones & neurotransmitters": [
      { name: "Dopamine",          smiles: "NCCc1ccc(O)c(O)c1" },
      { name: "Serotonin",         smiles: "NCCc1c[nH]c2ccc(O)cc12", synonyms: "5-HT" },
      { name: "Adrenaline",        smiles: "CNC[C@H](O)c1ccc(O)c(O)c1", synonyms: "epinephrine" },
      { name: "Noradrenaline",     smiles: "NC[C@H](O)c1ccc(O)c(O)c1", synonyms: "norepinephrine" },
      { name: "Acetylcholine",     smiles: "CC(=O)OCC[N+](C)(C)C" },
      { name: "GABA",              smiles: "NCCCC(=O)O", synonyms: "gamma-aminobutyric" },
      { name: "Glutamate",         smiles: "N[C@@H](CCC(=O)O)C(=O)O" },
      { name: "Histamine",         smiles: "NCCc1cnc[nH]1" },
      { name: "Melatonin",         smiles: "COc1ccc2[nH]cc(CCNC(C)=O)c2c1" },
      { name: "Testosterone",      smiles: "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O" },
      { name: "Estradiol",         smiles: "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc43)[C@@H]1CC[C@@H]2O" },
      { name: "Progesterone",      smiles: "CC(=O)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C" },
      { name: "Cortisol",          smiles: "O=C(CO)[C@@]1(O)CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3[C@@H](O)C[C@@]21C", synonyms: "hydrocortisone" },
      { name: "Thyroxine",         smiles: "Oc1cc(I)c(Oc2cc(I)c(C[C@H](N)C(=O)O)cc2I)c(I)c1", synonyms: "T4" },
    ],

    "Natural products": [
      { name: "Quinine",           smiles: "C=C[C@H]1CN2CC[C@H]1C[C@H]2[C@H](O)c1ccnc2ccc(OC)cc12", synonyms: "antimalarial" },
      { name: "Nicotine",          smiles: "CN1CCC[C@H]1c1cccnc1" },
      { name: "Capsaicin",         smiles: "COc1cc(CNC(=O)CCCCC/C=C/C(C)C)ccc1O" },
      { name: "Curcumin",          smiles: "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O" },
      { name: "Resveratrol",       smiles: "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1" },
      { name: "Quercetin",         smiles: "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", synonyms: "flavonoid" },
      { name: "Catechin",          smiles: "Oc1cc(O)c2C[C@@H](O)[C@H](c3ccc(O)c(O)c3)Oc2c1" },
      { name: "Cinnamaldehyde",    smiles: "O=C/C=C/c1ccccc1" },
      { name: "Vanillin",          smiles: "COc1cc(C=O)ccc1O" },
      { name: "Eugenol",           smiles: "C=CCc1ccc(O)c(OC)c1", synonyms: "clove" },
      { name: "Menthol",           smiles: "C[C@@H]1CC[C@@H](C(C)C)[C@H](O)C1" },
      { name: "Camphor",           smiles: "CC1(C)C2CCC1(C)C(=O)C2" },
      { name: "Limonene (R)",      smiles: "CC(=C)[C@@H]1CCC(C)=CC1", synonyms: "citrus terpene" },
      { name: "α-Pinene",          smiles: "CC1=CCC2CC1C2(C)C", synonyms: "pine terpene" },
      { name: "Linalool",          smiles: "C=CC(C)(O)CC/C=C(C)/C" },
      { name: "Geraniol",          smiles: "CC(C)=CCC/C(C)=C/CO" },
      { name: "Salicin",           smiles: "OC[C@H]1O[C@@H](Oc2ccccc2CO)[C@H](O)[C@@H](O)[C@@H]1O" },
      { name: "Theobromine",       smiles: "Cn1cnc2c1c(=O)[nH]c(=O)n2C", synonyms: "cocoa" },
      { name: "Theophylline",      smiles: "Cn1c(=O)c2[nH]cnc2n(C)c1=O" },
      { name: "Strychnine",        smiles: "O=C1CC2OCC=C3CN4CCC56C4CC3C2N1c1ccccc15" },
    ],

    "Polymer monomers": [
      { name: "Ethylene",          smiles: "C=C" },
      { name: "Propylene",         smiles: "CC=C" },
      { name: "Styrene",           smiles: "C=Cc1ccccc1" },
      { name: "Vinyl chloride",    smiles: "C=CCl" },
      { name: "Vinyl acetate",     smiles: "C=COC(C)=O" },
      { name: "Methyl methacrylate", smiles: "C=C(C)C(=O)OC", synonyms: "MMA acrylic" },
      { name: "Acrylic acid",      smiles: "C=CC(=O)O" },
      { name: "Acrylonitrile",     smiles: "C=CC#N" },
      { name: "Tetrafluoroethylene", smiles: "C(=C(F)F)(F)F", synonyms: "TFE PTFE teflon" },
      { name: "Caprolactam",       smiles: "O=C1CCCCCN1", synonyms: "nylon-6 monomer" },
      { name: "Bisphenol A",       smiles: "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1", synonyms: "BPA polycarbonate" },
      { name: "Terephthalic acid", smiles: "O=C(O)c1ccc(C(=O)O)cc1", synonyms: "PET" },
      { name: "Ethylene glycol",   smiles: "OCCO", synonyms: "PET" },
      { name: "Lactic acid",       smiles: "C[C@H](O)C(=O)O", synonyms: "PLA" },
      { name: "Glycolic acid",     smiles: "OCC(=O)O", synonyms: "PGA" },
    ],

    "Pesticides & herbicides": [
      { name: "Glyphosate",        smiles: "OC(=O)CNCP(=O)(O)O", synonyms: "roundup herbicide" },
      { name: "Atrazine",          smiles: "CCNc1nc(Cl)nc(NC(C)C)n1", synonyms: "herbicide" },
      { name: "DDT",               smiles: "Clc1ccc(C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl)cc1" },
      { name: "Malathion",         smiles: "CCOC(=O)CC(SP(=S)(OC)OC)C(=O)OCC" },
      { name: "Carbaryl",          smiles: "CNC(=O)Oc1cccc2ccccc12" },
      { name: "Paraquat",          smiles: "C[n+]1ccc(-c2cc[n+](C)cc2)cc1", synonyms: "herbicide" },
      { name: "2,4-D",             smiles: "OC(=O)COc1ccc(Cl)cc1Cl", synonyms: "herbicide" },
      { name: "Imidacloprid",      smiles: "O=N/[N-]C1=N/CN(Cc2ccc(Cl)nc2)CC1", synonyms: "neonicotinoid" },
      { name: "Permethrin",        smiles: "CC1(C)C(C=C(Cl)Cl)C1C(=O)OCc1cccc(Oc2ccccc2)c1", synonyms: "insecticide" },
    ],

    "Fragments (Rule of 3)": [
      { name: "Phenol",            smiles: "Oc1ccccc1" },
      { name: "Aniline",           smiles: "Nc1ccccc1" },
      { name: "Benzoic acid",      smiles: "OC(=O)c1ccccc1" },
      { name: "Benzaldehyde",      smiles: "O=Cc1ccccc1" },
      { name: "Acetanilide",       smiles: "CC(=O)Nc1ccccc1" },
      { name: "4-Methoxyphenol",   smiles: "COc1ccc(O)cc1" },
      { name: "Resorcinol",        smiles: "Oc1cccc(O)c1" },
      { name: "Catechol",          smiles: "Oc1ccccc1O" },
      { name: "Salicylic acid",    smiles: "OC(=O)c1ccccc1O" },
      { name: "1,3-Benzodioxole",  smiles: "C1Oc2ccccc2O1" },
      { name: "Indoline",          smiles: "C1Cc2ccccc2N1" },
      { name: "Tetrahydroquinoline", smiles: "C1CCc2ccccc2N1" },
      { name: "1,2,3,4-Tetrahydroisoquinoline", smiles: "C1CNCc2ccccc21" },
      { name: "N-Methylpiperazine", smiles: "CN1CCNCC1" },
      { name: "Methylcyclohexylamine", smiles: "CC1CCCCC1N" },
      { name: "Tropinone",         smiles: "O=C1CC2CCC(C1)N2C" },
    ],

    "Functional-group exemplars": [
      { name: "Methane",           smiles: "C" },
      { name: "Ethanol",           smiles: "CCO", synonyms: "alcohol" },
      { name: "Acetone",           smiles: "CC(=O)C", synonyms: "ketone" },
      { name: "Acetic acid",       smiles: "CC(=O)O", synonyms: "carboxylic" },
      { name: "Ethyl acetate",     smiles: "CCOC(C)=O", synonyms: "ester" },
      { name: "Acetamide",         smiles: "CC(=O)N", synonyms: "amide" },
      { name: "Methylamine",       smiles: "CN", synonyms: "primary amine" },
      { name: "Dimethylamine",     smiles: "CNC", synonyms: "secondary amine" },
      { name: "Trimethylamine",    smiles: "CN(C)C", synonyms: "tertiary amine" },
      { name: "Acetonitrile",      smiles: "CC#N", synonyms: "nitrile" },
      { name: "Nitromethane",      smiles: "C[N+](=O)[O-]", synonyms: "nitro" },
      { name: "Methanethiol",      smiles: "CS", synonyms: "thiol" },
      { name: "Dimethyl sulfide",  smiles: "CSC", synonyms: "sulfide" },
      { name: "Dimethyl sulfoxide", smiles: "CS(=O)C", synonyms: "sulfoxide" },
      { name: "Methanesulfonic acid", smiles: "CS(=O)(=O)O", synonyms: "sulfonic" },
      { name: "Phosgene",          smiles: "ClC(Cl)=O", synonyms: "acid chloride" },
      { name: "Acetyl chloride",   smiles: "CC(Cl)=O", synonyms: "acyl chloride" },
      { name: "Ethyl isocyanate",  smiles: "CCN=C=O", synonyms: "isocyanate" },
      { name: "Hydrazine",         smiles: "NN" },
      { name: "Hydroxylamine",     smiles: "NO" },
    ],

    "Stereo / chirality test": [
      { name: "(R)-CHFClBr",       smiles: "[C@H](F)(Cl)Br" },
      { name: "(S)-CHFClBr",       smiles: "[C@@H](F)(Cl)Br" },
      { name: "(R)-Lactic acid",   smiles: "C[C@H](O)C(=O)O" },
      { name: "(S)-Lactic acid",   smiles: "C[C@@H](O)C(=O)O" },
      { name: "Tartaric (RR)",     smiles: "O[C@H]([C@H](O)C(=O)O)C(=O)O" },
      { name: "meso-Tartaric",     smiles: "O[C@H]([C@@H](O)C(=O)O)C(=O)O" },
      { name: "cis-2-Butene",      smiles: "C/C=C\\C" },
      { name: "trans-2-Butene",    smiles: "C/C=C/C" },
      { name: "cis-1,2-DMCH",      smiles: "C[C@H]1CCCC[C@H]1C", synonyms: "dimethyl cyclohexane" },
      { name: "trans-1,2-DMCH",    smiles: "C[C@H]1CCCC[C@@H]1C" },
    ],
  };

  // -------------------------------------------------------------------- //
  // Public API
  // -------------------------------------------------------------------- //

  function list() {
    const out = [];
    Object.entries(CATALOG).forEach(([category, items]) => {
      items.forEach((entry) => out.push({ category, ...entry }));
    });
    return out;
  }

  function categories() {
    return Object.keys(CATALOG);
  }

  function byCategory() {
    // shallow copy so callers can't mutate the source-of-truth
    const out = {};
    Object.entries(CATALOG).forEach(([cat, items]) => {
      out[cat] = items.map((e) => ({ ...e }));
    });
    return out;
  }

  function score(query, entry) {
    if (!query) return 1;
    const q = query.toLowerCase();
    const name = (entry.name || "").toLowerCase();
    const synonyms = (entry.synonyms || "").toLowerCase();
    const smiles = entry.smiles || "";
    // Match priority:
    //   1. name substring (highest, weighted by position)
    //   2. synonyms substring
    //   3. exact SMILES match
    //   4. SMILES substring (only for queries ≥ 4 chars to avoid C / O hits)
    //   5. name sub-sequence (fallback)
    // Category is intentionally NOT matched: substrings like "caff"
    // accidentally hit "Scaffolds" and dragged in unrelated entries.
    if (name.includes(q)) return 1 - name.indexOf(q) / Math.max(1, name.length);
    if (synonyms.includes(q)) return 0.7;
    if (smiles === query) return 0.9;
    if (q.length >= 4 && smiles.includes(query)) return 0.5;
    let i = 0;
    for (const ch of name) {
      if (ch === q[i]) i += 1;
      if (i === q.length) return 0.25;
    }
    return 0;
  }

  function search(q, { limit = 200 } = {}) {
    const all = list();
    if (!q) return all;
    return all
      .map((e) => ({ e, s: score(q, e) }))
      .filter((r) => r.s > 0)
      .sort((a, b) => b.s - a.s)
      .slice(0, limit)
      .map((r) => r.e);
  }

  function count() {
    return list().length;
  }

  window.NuGenLibrary = { list, categories, byCategory, search, count };
})();
