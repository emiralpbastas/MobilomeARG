import pandas as pd
import re

def parse_gff3(file_input) -> pd.DataFrame:
    """Prokka veya Bakta GFF3 dosyasını okuyup gen tablosuna dönüştürür."""
    records = []
    if hasattr(file_input, "read"):
        raw = file_input.read()
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="ignore")
        lines = raw.splitlines()
    else:
        with open(file_input, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

    for line in lines:
        if line.startswith("##FASTA"):
            break
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        seqid, source, ftype, start, end, score, strand, phase, attributes = parts
        if ftype != "CDS":
            continue
        
        attr_dict = {}
        for item in attributes.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                attr_dict[k.strip()] = v.strip()
        
        records.append({
            "contig": seqid.strip(),
            "start": int(start),
            "end": int(end),
            "strand": strand.strip(),
            "gene_id": attr_dict.get("ID", ""),
            "gene_name": attr_dict.get("gene", attr_dict.get("Name", "hypothetical")),
            "product": attr_dict.get("product", "hypothetical protein")
        })
    return pd.DataFrame(records)

def parse_rgi(file_input) -> pd.DataFrame:
    """CARD / RGI (.txt / .tsv) çıktısını okuyup standartlaştırır."""
    df = pd.read_csv(file_input, sep="\t")
    def clean_contig(val):
        val = str(val).strip()
        return re.sub(r"^(contig_\d+)_\d+$", r"\1", val)

    df["contig_clean"] = df["Contig"].apply(clean_contig)
    return df

def parse_plasmidfinder(file_input) -> set:
    """Opsiyonel PlasmidFinder (.tsv/.txt) çıktısından plazmid pozitif contig isimlerini çıkarır."""
    if file_input is None:
        return set()
    df = pd.read_csv(file_input, sep="\t")
    for col in ["Contig", "contig", "Sequence", "SEQUENCE"]:
        if col in df.columns:
            return set(df[col].astype(str).str.strip().tolist())
    return set()
