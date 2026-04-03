// lib/genomics/promoter-cds-analyzer.ts

type CDS = {
  start: number
  end: number
  length: number
  frame: number
  sequence: string
}

type Promoter = {
  start: number
  end: number
  sequence: string
}

type CisElement = {
  name: string
  motif: string
  position: number
}

const START_CODON = "ATG"
const STOP_CODONS = ["TAA", "TAG", "TGA"]

// ---------------------------------------------
// CDS / ORF Finder
// ---------------------------------------------

export function findCDS(sequence: string): CDS[] {

  const seq = sequence.toUpperCase()
  const results: CDS[] = []

  for (let frame = 0; frame < 3; frame++) {

    for (let i = frame; i < seq.length; i += 3) {

      const codon = seq.substring(i, i + 3)

      if (codon === START_CODON) {

        for (let j = i + 3; j < seq.length; j += 3) {

          const stop = seq.substring(j, j + 3)

          if (STOP_CODONS.includes(stop)) {

            const cdsSeq = seq.substring(i, j + 3)

            results.push({
              start: i + 1,
              end: j + 3,
              length: cdsSeq.length,
              frame: frame + 1,
              sequence: cdsSeq
            })

            break
          }

        }

      }

    }

  }

  return results
}

// ---------------------------------------------
// Promoter Region Finder
// default: 1000bp upstream
// ---------------------------------------------

export function extractPromoter(sequence: string, cdsStart: number, upstream = 1000): Promoter {

  const start = Math.max(0, cdsStart - upstream)

  return {
    start,
    end: cdsStart - 1,
    sequence: sequence.slice(start, cdsStart)
  }

}

// ---------------------------------------------
// Cis-acting element database
// (common plant motifs)
// ---------------------------------------------

const CIS_DATABASE = [

  { name: "TATA Box",    motif: "TATAAA" },
  { name: "CAAT Box",    motif: "CCAAT"  },
  { name: "GC Box",      motif: "GGGCGG" },
  { name: "ABRE",        motif: "ACGTG"  },
  { name: "G-Box",       motif: "CACGTG" },
  { name: "ARE",         motif: "TGGTTT" },
  { name: "MYB Binding", motif: "CNGTTR" }

]

// ---------------------------------------------
// Cis element scanner
// ---------------------------------------------

export function findCisElements(sequence: string): CisElement[] {

  const seq = sequence.toUpperCase()

  const hits: CisElement[] = []

  CIS_DATABASE.forEach(element => {

    const motif = element.motif.replace("N", "[ATGC]")

    const regex = new RegExp(motif, "g")

    let match

    while ((match = regex.exec(seq)) !== null) {

      hits.push({
        name: element.name,
        motif: element.motif,
        position: match.index + 1
      })

    }

  })

  return hits

}
