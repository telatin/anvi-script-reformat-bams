# Define shell
SHELL:=/bin/bash

# Define your directories
INPUT_DIR:=test/input
OUTPUT_DIR:=test

# Define files
REF:=$(INPUT_DIR)/ref.fa
REF_INDEX:=$(REF).bwt
RAW_BAM:=$(INPUT_DIR)/raw.bam
REFORMAT_FA:=$(OUTPUT_DIR)/REFORMAT.fa
LIST:=$(OUTPUT_DIR)/REPORT_FILE.txt

# Output
BAM:=$(OUTPUT_DIR)/output/new.bam

# Define dependencies
DEPS:=bwa samtools anvi-script-reformat-fasta

# Default rule
all: check_dependencies $(BAM)

# Rule to check for dependencies
check_dependencies:
	@echo "Checking for required dependencies..."
	@$(foreach dep,$(DEPS),\
		which $(dep) > /dev/null 2>&1 || (echo "ERROR: $(dep) is not installed or not in the PATH"; exit 1);)

# Rule for creating the BAM file using bwa mem, depends on the reference index
$(RAW_BAM): $(REF_INDEX)
	bwa mem $(REF) $(INPUT_DIR)/reads_R1.fq.gz $(INPUT_DIR)/reads_R2.fq.gz | samtools sort --write-index -o $(RAW_BAM) -

# Rule for creating the reference index using bwa index
$(REF_INDEX): $(REF)
	bwa index $(REF)

# Rule for making REFORMAT.fa from ref.fa
$(REFORMAT_FA) $(LIST): $(REF)
	anvi-script-reformat-fasta -r $(LIST) $(REF) -o $(REFORMAT_FA) --simplify

$(BAM): $(RAW_BAM) $(LIST)
	./bin/anvi-reformat-bam -l $(LIST) -i $(RAW_BAM) -o $(BAM)
	samtools index $(BAM)
	tree > docs/output.tree

# Phony targets for cleanliness
.PHONY: all clean

# Clean rule to remove generated files
clean:
	rm -f $(RAW_BAM)* $(REFORMAT_FA) $(REF).* $(LIST) $(BAM)*
	tree > docs/input.tree
