#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ZIP_FASTA {

    tag "Zipping up your ${TYPE}s..."
    label "genelab"


    input:
        val(TYPE)
        path(DIR)

    output:
        path("*.zip"), emit: zip_files, optional: true
        path("versions.txt"), emit: version

    script:
        """
        function zip_sample() {

            local SAMPLE=\$1
            local TYPE=\$2

            mkdir -p \${SAMPLE}-\${TYPE}s && \\
            cp -f \${SAMPLE}-\${TYPE}*.fasta \${SAMPLE}-\${TYPE}s && \\
            zip -r \${SAMPLE}-\${TYPE}s.zip \${SAMPLE}-\${TYPE}s

           }


        export -f zip_sample

        if [ ${TYPE} == 'bin' ]; then 
       
              WORKDIR=`pwd`
        else

              WORKDIR=${DIR}
        fi

        if [ `find -L \${WORKDIR} -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then


             if [ ${TYPE} == 'MAG' ]; then

                 find -L \${WORKDIR} -name '*.fasta' | xargs -I {} cp {} .

             fi

             SAMPLES=(`ls -1 *.fasta | sed -E 's/(.+)-${TYPE}.*.fasta/\\1/g'`)

             for SAMPLE in \${SAMPLES[*]};do
                  
                  zip_sample \${SAMPLE} ${TYPE}

             done

        fi
        GL-version | grep "GeneLab utils"| sed -E 's/^\\s+//' | sed 's/\x1B\[[0-9;]\{1,\}[A-Za-z]//g' > versions.txt
        """

}
