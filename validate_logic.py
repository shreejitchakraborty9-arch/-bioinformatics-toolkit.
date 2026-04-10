import requests
import io
import sys
from Bio import SeqIO

def validate_cds(acc):
    try:
        # 1. Fetch Official NCBI Record
        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        query = "?db=nucleotide&id={0}&rettype=gb&retmode=text".format(acc)
        url = base + query
        
        print("Fetching {0} from NCBI...".format(acc))
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
        
        # In Python 2, use io.BytesIO or just pass the string if SeqIO supports it
        handle = io.BytesIO(resp.content)
        record = SeqIO.read(handle, "genbank")

        # 2. Ground Truth (Biopython Built-in Extraction)
        cds_feat = None
        for f in record.features:
            if f.type == "CDS":
                cds_feat = f
                break
        
        if not cds_feat:
            print("No CDS found for {0}".format(acc))
            return
            
        official_cds = str(cds_feat.extract(record.seq))

        # 3. Simulate App Logic (Coordinate-based splicing from our app.py and promoter.js)
        cds_regions = []
        # In Biopython locations might have parts
        if hasattr(cds_feat.location, 'parts'):
            for p in cds_feat.location.parts:
                cds_regions.append([int(p.start) + 1, int(p.end)])
        else:
            cds_regions.append([int(cds_feat.location.start) + 1, int(cds_feat.location.end)])
            
        strand = "+" if cds_feat.location.strand >= 0 else "-"
        
        # promoter.js extractCDS logic:
        if strand == "+":
            sorted_regions = sorted(cds_regions, key=lambda x: x[0])
        else:
            sorted_regions = sorted(cds_regions, key=lambda x: x[0], reverse=True)
        
        app_cds = ""
        for start, end in sorted_regions:
            # JavaScript seq.substring(start - 1, end)
            chunk = str(record.seq)[start-1:end]
            if strand == "-":
                from Bio.Seq import reverse_complement
                chunk = str(reverse_complement(chunk))
            app_cds += chunk

        # 4. Comparison
        length_match = (len(official_cds) == len(app_cds))
        mismatches = []
        min_len = min(len(official_cds), len(app_cds))
        for i in range(min_len):
            if official_cds[i] != app_cds[i]:
                mismatches.append((i, official_cds[i], app_cds[i]))
        
        identity = ((len(official_cds) - len(mismatches)) * 100.0 / len(official_cds)) if len(official_cds) > 0 else 0
        
        print("\n--- Validation Report for {0} ---".format(acc))
        print("Strand: {0}".format(strand))
        print("Official CDS Length: {0} bp".format(len(official_cds)))
        print("App Generated CDS Length: {0} bp".format(len(app_cds)))
        print("Identity: {0:.4f}%".format(identity))
        print("Mismatches: {0}".format(len(mismatches)))
        
        if mismatches:
            print("Example mismatches (first 5): {0}".format(mismatches[:5]))
        
    except Exception as e:
        print("Validation failed for {0}: {1}".format(acc, str(e)))

if __name__ == "__main__":
    validate_cds("NM_007294.4") # BRCA1 (+)
    validate_cds("NM_003106.4") # SOX2 (-)
