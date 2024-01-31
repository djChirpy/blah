#!/bin/bash                                                                    

source /broad/software/scripts/useuse
use Java-1.8
export SINGULARITY_CACHEDIR=/seq/vgb/swofford/singularity/cache
export SINGULARITY_TMPDIR=/seq/vgb/swofford/singularity/temp2
export SINGULARITY_LOCALCACHEDIR=/seq/vgb/swofford/singularity/temp2
export SINGULARITY_DOCKER_USERNAME=djchirpy
export SINGULARITY_DOCKER_PASSWORD=dockerIsAnnoying52

#$ -l h_vmem=30G
#$ -o /seq/vgb/swofford/logs/
#$ -e /seq/vgb/swofford/logs/
#$ -M swofford@broadinstitute.org
#$ -m ea
#$ -l h_rt=48:00:00
#$ -l os=RedHat7
#$ -V

cd /seq/vgb/swofford/TEs/gatkSVPipeline/ugerTest

#cd /seq/vgb/swofford/TEs/gatkSVPipeline/Local

java -Dconfig.file=/seq/vgb/swofford/TEs/gatkSVPipeline/cromwell.conf.uger -Xmx20G -jar /seq/vgb/software/cromwell-44.jar run /seq/vgb/swofford/TEs/gatkSVPipeline/gittingud/gatk-sv-v1-dog/module00a/Module00a.wdl -i /seq/vgb/swofford/TEs/gatkSVPipeline/Module00a.json

#java -Dconfig.file=/seq/vgb/swofford/TEs/gatkSVPipeline/cromwell.conf.uger -Xmx20G -jar /seq/vgb/software/cromwell-44.jar run /seq/vgb/swofford/TEs/gatkSVPipeline/activeEditsBasedOnJessica/gatk-sv-v1-master/module00a/Module00a.wdl -i /seq/vgb/swofford/TEs/gatkSVPipeline/Module00a.json

#java -Dconfig.file=/seq/vgb/swofford/TEs/gatkSVPipeline/cromwell.conf.nosubmit -Xmx20G -jar /seq/vgb/software/cromwell-44.jar run /seq/vgb/swofford/TEs/gatkSVPipeline/gittingud/gatk-sv-v1-dog/module00a/Module00a.wdl -i /seq/vgb/swofford/TEs/gatkSVPipeline/Module00a.json
