!pip install -q biopython py3Dmol ipywidgets requests

from IPython.display import display, HTML, clear_output
import ipywidgets as widgets
import requests, json, warnings
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
import pandas as pd
import py3Dmol
from io import StringIO
from Bio import BiopythonWarning

Entrez.email = "your_email@example.com"
warnings.filterwarnings('ignore', category=BiopythonWarning)

RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{}"
RCSB_FASTA_URL = "https://www.rcsb.org/fasta/entry/{}"
RCSB_CHEM_URL  = "https://data.rcsb.org/rest/v1/core/chemcomp/{}"

pdb_input = widgets.Text(value='4R4V', description='PDB ID:', placeholder='e.g., 4R4V')
run_button = widgets.Button(description='Run Query', button_style='primary')
out = widgets.Output(layout={'border': '1px solid black'})

def fetch_rcsb_entry(pdb_id):
    url = RCSB_ENTRY_URL.format(pdb_id)
    r = requests.get(url)
    r.raise_for_status()
    return r.json()

def fetch_fasta(pdb_id):
    url = RCSB_FASTA_URL.format(pdb_id)
    r = requests.get(url)
    r.raise_for_status()
    return r.text

def parse_fasta_text(fasta_text):
    records = []
    for block in fasta_text.strip().split('>'):
        if not block.strip():
            continue
        header, *seq_lines = block.splitlines()
        seq = ''.join(seq_lines).replace(' ', '').replace('\r','').strip()
        records.append({'id': header.split()[0], 'desc': header, 'seq': seq})
    return records

def try_get_refseq_from_entry(entry_json):
    refseqs = []
    polymer_entities = entry_json.get('rcsb_polymer_entity_container_identifiers', {})
    ref_ids = polymer_entities.get('reference_sequence_identifiers', [])
    for r in ref_ids:
        db = r.get('database')
        accession = r.get('identifier')
        if db and accession:
            refseqs.append((db, accession))
    return refseqs

def fetch_ncbi_nucleotide(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        fasta_text = handle.read()
        handle.close()
        fh = StringIO(fasta_text)
        record = SeqIO.read(fh, "fasta")
        return {'id': record.id, 'desc': record.description, 'seq': str(record.seq)}
    except Exception as e:
        return None

def fetch_chem_info(chem_id):
    try:
        r = requests.get(RCSB_CHEM_URL.format(chem_id))
        r.raise_for_status()
        data = r.json()
        name = data.get('chemicalName') or data.get('name') or chem_id
        formula = data.get('formula')
        synonyms = data.get('synonyms', [])
        return {'chem_id': chem_id, 'name': name, 'formula': formula, 'synonyms': synonyms}
    except Exception as e:
        return {'chem_id': chem_id, 'name': chem_id}

def show_gc_bar(gc_pct):
    fig, ax = plt.subplots(figsize=(4,3))
    ax.bar(["GC%","Remaining"], [gc_pct, 100-gc_pct])
    ax.set_ylim(0,100)
    ax.set_ylabel("Percentage (%)")
    ax.set_title("GC Content (against 100%)")
    for i, v in enumerate([gc_pct, 100-gc_pct]):
        ax.text(i, v+1, f"{v:.2f}%", ha='center')
    plt.show()

def show_length_histogram(lengths):
    fig, ax = plt.subplots(figsize=(5,3))
    ax.hist(lengths, bins=max(3, min(12, len(lengths))), edgecolor='black')
    ax.set_xlabel("Sequence length (aa)")
    ax.set_ylabel("Count")
    ax.set_title("Histogram of Chain Sequence Lengths")
    plt.show()

def create_py3dmol_viewer(pdb_id, ligands_info):
    view = py3Dmol.view(query=f"pdb:{pdb_id}")
    view.setBackgroundColor('0x000000')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.addStyle({'hetflag': True}, {'stick': {'radius':0.2}, 'sphere': {'scale':0.3}})
    js_callback = """
    function(atom, viewer) {
      if(!atom) return;
      var info = '';
      if (atom.resn) info += atom.resn + ' ';
      if (atom.chain) info += 'chain:' + atom.chain + ' ';
      if (atom.serial) info += 'serial:' + atom.serial + ' ';
      if (atom.elem) info += 'elem:' + atom.elem + ' ';
      viewer.removeAllLabels();
      try {
        viewer.addLabel(info, {position: atom, backgroundColor: 'black', fontColor:'white', inFront:true});
      } catch(err) {}
    }
    """
    view.setHoverable({}, True, js_callback)
    view.zoomTo()
    return view

def run_query(b):
    pdb_id = pdb_input.value.strip().upper()
    if not pdb_id:
        with out:
            clear_output()
            print("Please enter a PDB ID (e.g., 4R4V).")
        return

    with out:
        clear_output()
        print(f"Querying RCSB for PDB ID: {pdb_id} ...")
        try:
            entry = fetch_rcsb_entry(pdb_id)
        except Exception as e:
            print("Error fetching RCSB entry:", e)
            return

        try:
            fasta_text = fetch_fasta(pdb_id)
            fasta_records = parse_fasta_text(fasta_text)
        except Exception as e:
            fasta_records = []
            print("Failed to fetch FASTA from RCSB:", e)

        if fasta_records:
            print("\nFASTA Sequences (from RCSB) — showing header and first 80 aa of each:")
            for rec in fasta_records:
                print(">", rec['id'], rec['desc'])
                seq_snip = rec['seq'][:80] + ("..." if len(rec['seq'])>80 else "")
                print(seq_snip)
        else:
            print("No FASTA sequences found for this entry.")

        chain_lengths = [len(rec['seq']) for rec in fasta_records] if fasta_records else []
        if chain_lengths:
            print("\nChain lengths:", chain_lengths)
            show_length_histogram(chain_lengths)
        else:
            print("\nNo polymer chain sequences to plot histogram.")

        refseqs = try_get_refseq_from_entry(entry)
        nucleotide_seq_record = None
        if refseqs:
            candidate = None
            for db, acc in refseqs:
                if db.lower() in ('refseq', 'refseq_nucleotide', 'refseq_peptide') and acc.startswith('NM_'):
                    candidate = acc
                    break
            if not candidate:
                for db, acc in refseqs:
                    if acc.startswith(('NM_','XM_','NC_','LR_','LRG_','NZ_')):
                        candidate = acc
                        break
            if candidate:
                print(f"\nFound reference nucleotide accession in RCSB metadata: {candidate}")
                nucleotide_seq_record = fetch_ncbi_nucleotide(candidate)
                if nucleotide_seq_record:
                    print("Fetched nucleotide record:", nucleotide_seq_record['id'])
                else:
                    print("Failed to fetch nucleotide record for", candidate)
            else:
                print("\nRCSB metadata contains reference identifiers but no suitable RefSeq accession found.")
        else:
            print("\nNo reference sequence identifiers found in RCSB entry metadata.")

        if nucleotide_seq_record:
            nuc_seq = nucleotide_seq_record['seq']
            rna_seq = nuc_seq.replace('T','U').replace('t','u')
            gc_pct = gc_fraction(nuc_seq) * 100
            print(f"\nNucleotide (DNA) length: {len(nuc_seq)}")
            print("Showing first 200 bases of RNA sequence:")
            print(rna_seq[:200] + ("..." if len(rna_seq)>200 else ""))
            show_gc_bar(gc_pct)
            trim_len = len(nuc_seq) - (len(nuc_seq) % 3)
            coding = nuc_seq[:trim_len]
            protein_from_nuc = str(Seq(coding).translate(to_stop=True))
            print("\nTranslated protein (first 200 aa):")
            print(protein_from_nuc[:200] + ("..." if len(protein_from_nuc)>200 else ""))
        else:
            if fasta_records:
                prot_seq = fasta_records[0]['seq']
                print("\nNo nucleotide sequence available.")
            else:
                print("\nNo sequences available.")

        if fasta_records:
            df_prots = pd.DataFrame([{'id': r['id'], 'length': len(r['seq']), 'description': r['desc']} for r in fasta_records])
            display(df_prots)

        ligands = []
        chem_list = set()
        if 'nonpolymer_entities' in entry:
            for np in entry['nonpolymer_entities']:
                chem = np.get('chem_comp')
                if chem:
                    chem_list.add(chem)
        try:
            insts = entry.get('nonpolymer_entity_instances', [])
            for inst in insts:
                chem_id = inst.get('chem_comp')
                if chem_id:
                    chem_list.add(chem_id)
        except:
            pass

        ligands_info = []
        if chem_list:
            print("\nLigands found:")
            for chem in sorted(chem_list):
                info = fetch_chem_info(chem)
                ligands_info.append(info)
                print(f"- {chem}: {info.get('name')} ({info.get('formula')})")
            display(pd.DataFrame(ligands_info))
        else:
            print("\nNo ligand chemical components found.")

        print("\nRendering interactive 3D structure...")
        try:
            view = create_py3dmol_viewer(pdb_id, ligands_info)
            display(view.show())
        except Exception as e:
            print("Error creating viewer:", e)

        print("\nQuery complete.")

run_button.on_click(run_query)

display(widgets.HBox([pdb_input, run_button]))
display(out)
