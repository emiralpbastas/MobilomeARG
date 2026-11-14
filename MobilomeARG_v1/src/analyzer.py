import pandas as pd
import re

MGE_PATTERNS = {
    "Transposon / IS Elementi": r"\b(transposase|insertion sequence|is\d+|iskpn\d+|tn\d+|ista|istb|tnpa|tnpr)\b",
    "İntegron / Rekombinaz": r"\b(integrase|inti\d*|recombinase|resolvase|site-specific)\b",
    "Plazmid / Konjugasyon": r"\b(plasmid|relaxase|mobilization|conjugal|conjugative|tra[a-z0-9]+|mob[a-z0-9]+|virb\d*|rep[a-z0-9]+)\b",
    "Bakteriyofaj (Prophage)": r"\b(phage|prophage|capsid|tail fiber|terminase|portal protein)\b"
}

FALSE_POSITIVE_GENES = {"dnaa", "gyrb", "rpob", "reca", "teta", "tetr", "sul1", "qacedelta1", "aada2", "amra", "amrb", "opra"}

def classify_mge(gene_name: str, product: str) -> str:
    g_clean = str(gene_name).strip().lower()
    if g_clean in FALSE_POSITIVE_GENES:
        return "Hücresel / Diğer"
    text = f"{gene_name} {product}".lower()
    for category, pattern in MGE_PATTERNS.items():
        if re.search(pattern, text, re.IGNORECASE):
            return category
    return "Hücresel / Diğer"

def analyze_hgt_risk(gff_df: pd.DataFrame, rgi_df: pd.DataFrame, window_bp: int = 10000, plasmid_contigs: set = None):
    if plasmid_contigs is None:
        plasmid_contigs = set()

    gff_df = gff_df.copy()
    gff_df["element_type"] = gff_df.apply(lambda r: classify_mge(r["gene_name"], r["product"]), axis=1)
    
    results = []
    neighborhood_records = {}

    for idx, arg in rgi_df.iterrows():
        contig = arg["contig_clean"]
        arg_start = int(arg["Start"])
        arg_stop = int(arg["Stop"])
        arg_name = arg["Best_Hit_ARO"]
        
        win_start = max(1, arg_start - window_bp)
        win_stop = arg_stop + window_bp
        
        neighbors = gff_df[
            (gff_df["contig"] == contig) &
            (gff_df["end"] >= win_start) &
            (gff_df["start"] <= win_stop)
        ].copy().sort_values("start")
        
        neighbors["is_target_arg"] = neighbors.apply(
            lambda r: abs(r["start"] - arg_start) <= 15 and abs(r["end"] - arg_stop) <= 15, axis=1
        )
        
        other_args = rgi_df[
            (rgi_df["contig_clean"] == contig) &
            (rgi_df["Best_Hit_ARO"] != arg_name) &
            (rgi_df["Stop"] >= win_start) &
            (rgi_df["Start"] <= win_stop)
        ]
        other_arg_coords = list(zip(other_args["Start"], other_args["Stop"]))
        
        def tag_visual_role(row):
            if row["is_target_arg"]:
                return "Hedef ARG (Direnç Geni)"
            for s, e in other_arg_coords:
                if abs(row["start"] - s) <= 15 and abs(row["end"] - e) <= 15:
                    return "Komşu ARG (Direnç Geni)"
            return row["element_type"]

        neighbors["visual_role"] = neighbors.apply(tag_visual_role, axis=1)
        mge_neighbors = neighbors[(~neighbors["is_target_arg"]) & (neighbors["element_type"] != "Hücresel / Diğer")].copy()
        
        score = 0
        evidence = []
        
        has_plasmid_gene = (mge_neighbors["element_type"] == "Plazmid / Konjugasyon").any()
        in_plasmidfinder = contig in plasmid_contigs
        if has_plasmid_gene or in_plasmidfinder:
            score += 4
            details = []
            if has_plasmid_gene:
                genes = ", ".join(mge_neighbors[mge_neighbors["element_type"] == "Plazmid / Konjugasyon"]["gene_name"].tolist())
                details.append(f"Konjugasyon/Plazmid genleri ({genes})")
            if in_plasmidfinder:
                details.append("PlasmidFinder pozitif contig")
            evidence.append(f"{' & '.join(details)} [+4]")
            
        transposons = mge_neighbors[mge_neighbors["element_type"] == "Transposon / IS Elementi"]
        if not transposons.empty:
            min_dist = min(
                min(abs(row["end"] - arg_start), abs(row["start"] - arg_stop))
                for _, row in transposons.iterrows()
            )
            t_genes = ", ".join(transposons["gene_name"].tolist())
            if min_dist <= 5000:
                score += 3
                evidence.append(f"Yakın Transpozon/IS ({t_genes}, mesafe: {min_dist} bp) [+3]")
            else:
                score += 2
                evidence.append(f"Orta mesafe Transpozon/IS ({t_genes}, mesafe: {min_dist} bp) [+2]")
            
            upstream_is = transposons[transposons["end"] <= arg_start]
            downstream_is = transposons[transposons["start"] >= arg_stop]
            if not upstream_is.empty and not downstream_is.empty:
                score += 3
                evidence.append("Bileşik Transpozon Yapısı (Sağ ve soldan IS/Transpozon ile çevrili sandviç yapı) [+3]")
                
        integrons = mge_neighbors[mge_neighbors["element_type"] == "İntegron / Rekombinaz"]
        if not integrons.empty:
            score += 2
            i_genes = ", ".join(integrons["gene_name"].tolist())
            evidence.append(f"İntegron/İntegraz komşuluğu ({i_genes}) [+2]")
            
        if not other_args.empty:
            score += 1
            co_genes = ", ".join(other_args["Best_Hit_ARO"].tolist())
            evidence.append(f"Çoklu direnç kümesi ({co_genes}) [+1]")
            
        score = min(10, score)
        if score >= 7:
            risk_level = "YÜKSEK RİSK (Kritik Mobilite)"
        elif score >= 3:
            risk_level = "ORTA RİSK (Potansiyel Mobil)"
        else:
            risk_level = "DÜŞÜK RİSK (Kromozomal / Sabit)"
            if not evidence:
                evidence.append("Komşulukta MGE yok; intrinsik/kromozomal bölge [0]")
                
        results.append({
            "ARG_Geni": arg_name,
            "Contig": contig,
            "Koordinat": f"{arg_start}-{arg_stop}",
            "Ilac_Sinifi": arg.get("Drug Class", "Bilinmiyor"),
            "Mekanizma": arg.get("Resistance Mechanism", "Bilinmiyor"),
            "HGT_Risk_Skoru": score,
            "Risk_Seviyesi": risk_level,
            "Komsu_MGE_Sayisi": len(mge_neighbors),
            "Kanit_ve_Puan_Detayi": " | ".join(evidence)
        })
        neighborhood_records[f"{arg_name}_{contig}"] = neighbors
        
    return pd.DataFrame(results), neighborhood_records
