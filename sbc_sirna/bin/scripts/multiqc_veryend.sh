#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load python/3.6.8

cd ../../data/

multiqc --filename multiqc_mapping_quant_picard \
	aligned/pass2/ \
	featureCounts/ \
	picard/
