# Define shell
SHELL:=/bin/bash

# Define your directories

TEST_DIR:=test
INPUT_DIR:=test/input
OUTPUT_DIR:=test/output

# Define files
REF:=$(INPUT_DIR)/ref.fa
REF_INDEX:=$(REF).bwt
RAW_BAM:=$(INPUT_DIR)/raw.bam
REFORMAT_FA:=$(TEST_DIR)/REFORMAT.fa
LIST:=$(TEST_DIR)/REPORT_FILE.txt

# Output
BAM0:=$(OUTPUT_DIR)/stub.bam
BAM1:=$(OUTPUT_DIR)/final.bam
# Define dependencies
DEPS:=bwa samtools anvi-script-reformat-fasta

# Default rule
all: check_dependencies $(BAM1)

# Rule to check for dependencies
check_dependencies:
	@echo "Checking for required dependencies..."
	@$(foreach dep,$(DEPS),\
		which $(dep) > /dev/null 2>&1 || (echo "ERROR: $(dep) is not installed or not in the PATH"; exit 1);)
	@mkdir -p $(OUTPUT_DIR) 

# Rule for creating the BAM file using bwa mem, depends on the reference index
$(RAW_BAM): $(REF_INDEX)
	bwa mem $(REF) $(INPUT_DIR)/reads_R1.fq.gz $(INPUT_DIR)/reads_R2.fq.gz | samtools sort --write-index -o $(RAW_BAM) -

# Rule for creating the reference index using bwa index
$(REF_INDEX): $(REF)
	bwa index $(REF)

# Rule for making REFORMAT.fa from ref.fa
$(REFORMAT_FA) $(LIST): $(REF)
	anvi-script-reformat-fasta -r $(LIST) $(REF) -o $(REFORMAT_FA) --simplify-names 

$(BAM0) $(BAM1): $(RAW_BAM) $(LIST)
	./bin/reformat-bam-stub        -l $(LIST) -i $(RAW_BAM) -o $(BAM0)
	./bin/anvi-script-reformat-bam -l $(LIST) --force -o $(BAM1) $(RAW_BAM)
	samtools index --csi $(BAM1)
	samtools index --csi $(BAM0)
	tree > docs/output.tree

# Phony targets for cleanliness
.PHONY: all clean

# Clean rule to remove generated files
clean:
	rm -f $(RAW_BAM)* $(REFORMAT_FA) $(REF).* $(LIST) $(BAM0)* $(BAM1)*
	tree > docs/input.tree
