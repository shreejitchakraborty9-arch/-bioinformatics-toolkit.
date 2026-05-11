import json
from Bio import Restriction as rst

common_names = ["EcoRI", "BamHI", "HindIII", "SalI", "XhoI", "NotI", "XbaI", "SpeI", "NheI", "BglII", "ClaI", "KpnI", "SacI", "MluI", "NcoI", "NdeI", "PstI", "PvuII", "ScaI", "StuI", "SmaI", "AvaI", "SphI", "ApaI", "HaeIII", "TaqI", "HinfI", "MspI", "HpaII", "AluI", "RsaI", "SfcI"]

IUPAC_PROB = {
    'A': 0.25, 'T': 0.25, 'G': 0.25, 'C': 0.25, 'U': 0.25,
    'R': 0.5,  'Y': 0.5,  'S': 0.5,  'W': 0.5,  'K': 0.5,  'M': 0.5,
    'B': 0.75, 'D': 0.75, 'H': 0.75, 'V': 0.75,
    'N': 1.0,
}

def calc_cut_frequency(site):
    prob = 1.0
    for base in site.upper():
        prob *= IUPAC_PROB.get(base, 0.25)
    return prob

enzymes = []
seen = set()

for enz_class in rst.CommOnly:
    try:
        name = enz_class.__name__
        site = enz_class.site
        fst5 = enz_class.fst5
        fst3 = enz_class.fst3
        blunt = enz_class.is_blunt()

        if blunt:
            overhang_type = "blunt"
        elif fst5 < fst3:
            overhang_type = "5' overhang"
        else:
            overhang_type = "3' overhang"

        cType = str(len(site)) + "-cutter"

        tags = []
        if name in common_names:
            tags.append("common")
        if len(site) == 4:
            tags.append("4-cutter")
        elif len(site) == 6:
            tags.append("6-cutter")
        if blunt:
            tags.append("blunt")
        else:
            tags.append("sticky")

        enzymes.append({
            "name": name,
            "site": site,
            "cut": fst5,
            "cut_sense": fst5,
            "cut_antisense": fst3,
            "type": cType,
            "overhang_type": overhang_type,
            "is_common": name in common_names,
            "tags": tags,
            "supplier": "NEB",
            "cut_frequency": calc_cut_frequency(site)
        })
        seen.add(name)
    except Exception as e:
        print("Error on", enz_class, e)
        continue

# Ensure every enzyme in common_names is present even if absent from CommOnly
for name in common_names:
    if name in seen:
        continue
    try:
        enz_class = getattr(rst, name)
        site = enz_class.site
        fst5 = enz_class.fst5
        fst3 = enz_class.fst3
        blunt = enz_class.is_blunt()

        if blunt:
            overhang_type = "blunt"
        elif fst5 < fst3:
            overhang_type = "5' overhang"
        else:
            overhang_type = "3' overhang"

        cType = str(len(site)) + "-cutter"

        tags = ["common"]
        if len(site) == 4:
            tags.append("4-cutter")
        elif len(site) == 6:
            tags.append("6-cutter")
        if blunt:
            tags.append("blunt")
        else:
            tags.append("sticky")

        enzymes.append({
            "name": name,
            "site": site,
            "cut": fst5,
            "cut_sense": fst5,
            "cut_antisense": fst3,
            "type": cType,
            "overhang_type": overhang_type,
            "is_common": True,
            "tags": tags,
            "supplier": "NEB",
            "cut_frequency": calc_cut_frequency(site)
        })
        seen.add(name)
    except Exception as e:
        print("Error adding missing common enzyme", name, e)
        continue

output_path = "js/restriction-enzyme-db.js"
with open(output_path, "w") as f:
    f.write("/* Auto-generated Restriction Enzyme Database */\r\n")
    f.write("window.BioKit = window.BioKit || { utils: {}, core: {}, tools: {}, data: {} };\r\n")
    f.write("window.RESTRICTION_ENZYMES = ")
    json.dump(enzymes, f, indent=2)
    f.write(";\n")
    f.write("\nwindow.BioKit.data = window.BioKit.data || {};\n")
    f.write("window.BioKit.data.RESTRICTION_ENZYMES = window.RESTRICTION_ENZYMES;\n")

print("Generated " + str(len(enzymes)) + " enzymes")
