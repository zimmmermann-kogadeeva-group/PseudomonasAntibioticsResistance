
all: reference.fasta.gz annotation.tsv OtherGenomes

OTHER_GENOMES_URL='https://www.pseudomonas.com/downloads/pseudomonas/pgd_r_22_1/Pseudomonas/complete/fna-complete.tar.gz'
OtherGenomes:
	wget ${OTHER_GENOMES_URL} && gunzip fna-complete.tar.gz && mv fna $@ && rm fna-complete.tar.gz

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

