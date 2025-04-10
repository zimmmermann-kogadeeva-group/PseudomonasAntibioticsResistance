
all: reference.fasta.gz annotation.tsv OtherGenomes

ILLUMNA_ADAPTERS_URL = 'https://cloud.stevenschmitt.de/s/axhzCB8h7CZKvui/download?path=%2FDocumentation%2FHelpful%20Files&files=Illumina_Adapters.fasta'
illumina_adapters.fasta:
	wget ${ILLUMNA_ADAPTERS_URL} -O $@

REF_URL = 'ftp://ftp.ensemblgenomes.org/pub/bacteria/release-52/fasta/bacteria_5_collection/pseudomonas_aeruginosa_pao1_gca_000006765/cdna/Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1.cdna.all.fa.gz'
reference.fasta.gz: 
	wget ${REF_URL} -O $@

ANNOT_URL = 'ftp://ftp.ensemblgenomes.org/pub/release-52/bacteria//tsv/bacteria_5_collection/pseudomonas_aeruginosa_pao1_gca_000006765/Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1.52.ena.tsv.gz'
annotation.tsv:
	wget ${ANNOT_URL} -O $@.gz && gunzip $@.gz

ANNOT_REF_URL = 'https://www.pseudomonas.com/downloads/pseudomonas/pgd_r_22_1/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107.gbk.gz'
PAO1.gbk:
	wget ${ANNOT_REF_URL} -O $@.gz && gunzip $@.gz

uniprot_mapping.csv:
	./get_uniprot_mapping.py $@


# Phylogenetic trees

OTHER_GENOMES_URL='https://www.pseudomonas.com/downloads/pseudomonas/pgd_r_22_1/Pseudomonas/complete/fna-complete.tar.gz'
OtherGenomes:
	wget ${OTHER_GENOMES_URL} && gunzip fna-complete.tar.gz && mv fna $@ && rm fna-complete.tar.gz

OTHER_GENOMES_META_URL='https://www.pseudomonas.com/downloads/pseudomonas/pgd_r_22_1/strain_summary.txt'
OtherGenomesMetadata/metadata.tsv:
	mkdir -p $(dir $@) && wget ${OTHER_GENOMES_META_URL} -O $@

# Public datasets

MEM_1_URL = 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE167137&format=file&file=GSE167137%5FP%5Faeruginosa%5Fexpression%5FPA%5F0M%5Fvs%5FPA%5F5M%5F%2Exlsx'
PublicDatasets/GSE167137/GSE167137_P_aeruginosa_expression_PA_0M_vs_PA_5M_.xlsx:
	mkdir -p $(dir $@) && wget ${MEM_1_URL} -O $@

MEM_3_URL = 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123544&format=file&file=GSE123544%5FProcessedDataMatrix%5Fclinical%5Fisolates%2Exlsx'
PublicDatasets/GSE123544/GSE123544_ProcessedDataMatrix_clinical_isolates.xlsx:
	mkdir -p $(dir $@) && wget ${MEM_3_URL} -O $@
