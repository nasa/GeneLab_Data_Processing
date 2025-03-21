#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//params.CAT_DL_LINK = "https://tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20210107.tar.gz"

/**************************************************************************************** 
***************************  Metagenomics databases set-up ******************************
****************************************************************************************/

// This process download CAT reference database.
process SETUP_CAT_DB {

    tag "Downloading and setting up contig annotation tool-s (CAT) database..."
    label "db_setup"

    input:
        val(CAT_DB_LINK)
    output:
        path("CAT_prepare_20210107/"), emit: cat_db
        path("CAT_prepare_20210107/CAT_DB_SETUP"), emit: completion_indicator
        path("versions.txt"), emit: version
    script:
        """
        printf "### Setting up CAT's reference database ###\\n\\n"
        printf "  Downloading reference db:\\n\\n"
        curl -L -C - -o CAT_prepare_20210107.tar.gz ${CAT_DB_LINK}

        printf "\\n\\n  Extracting reference db:\\n\\n"
        tar -xvzf CAT_prepare_20210107.tar.gz

        rm CAT_prepare_20210107.tar.gz CAT_prepare_20210107/2021-01-07_CAT_database/2021-01-07.nr.gz
        touch CAT_prepare_20210107/CAT_DB_SETUP
        curl --version  |head -n 1 | sed -E 's/(curl\\s.+)\\s\\(.+/\\1/' > versions.txt
        printf "### Set up completed successfully ###\\n\\n"
        """
}

// This process downloads KOFamScan db (minimally currently).
process SETUP_KOFAMSCAN_DB {

    tag "Downloading and setting up kofam scan-s database..."
    label "db_setup"

    output:
        path("kofamscan_db/"), emit: ko_db_dir
        path("kofamscan_db/KO_DB_SETUP"), emit: completion_indicator
        path("versions.txt"), emit: version
    script:
        """
        printf "### Setting up KOFamScan reference database ###\\n\\n" 
        
        # Using https instead of ftp for those whose systems don't have access to the ftp servers

        printf "\\n  Downloading ko_list file:\\n\\n"

        if ! curl -L -C - --connect-timeout 15 -o ko_list.gz ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz 
        then
            printf "\\n\\n  Downloading via http since ftp seemed to fail making the connection:\\n\\n"
            curl -L -C - -o ko_list.gz https://www.genome.jp/ftp/db/kofam/ko_list.gz
        fi

        printf "\\n\\n  Downloading profiles.tar.gz file:\\n\\n"


        if ! curl -L -C - --connect-timeout 15 -o profiles.tar.gz ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
        then
            printf "\\n\\n  Downloading via http since ftp seemed to fail making the connection:\\n\\n"
            curl -L -C - -o profiles.tar.gz https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
        fi

       [ -d kofamscan_db/ ] || mkdir kofamscan_db/
        printf "\\n\\n  Decompressing profiles.tar.gz file:\\n\\n" 
        tar -xzf profiles.tar.gz -C kofamscan_db/
        rm profiles.tar.gz

        gunzip ko_list.gz && \\
        mv ko_list kofamscan_db/ && \\
        touch kofamscan_db/KO_DB_SETUP
        curl --version  |head -n 1 | sed -E 's/(curl\\s.+)\\s\\(.+/\\1/' > versions.txt
        printf "### Set up completed successfully ###\\n\\n"
        """
}

// This process downloads the gtdb-tk db (minimally currently)
process SETUP_GTDBTK_DB {

    tag "Downloading and setting up genome taxonomy database toolkit-s (GTDBTK) database..."
    label "db_setup"

    input:
        val(GTDBTK_URL)
    output:
        path("GTDB-tk-ref-db/"), emit: gtdbtk_db_dir
        path("GTDB-tk-ref-db/SETUP_GTDBTK_DB_SETUP"), emit: completion_indicator
        path("versions.txt"), emit: version
    script:
        """
        [ -d GTDB-tk-ref-db/ ] || mkdir -p GTDB-tk-ref-db/

        # But still needs to be set for this particular session that is downloading and setting up the db
        export GTDBTK_DATA_PATH=GTDB-tk-ref-db/
        # Downloading
        download-GTDBTK-db.sh ${GTDBTK_URL} && touch GTDB-tk-ref-db/SETUP_GTDBTK_DB_SETUP
        gtdbtk -h |grep "GTDB-Tk" | sed -E 's/.+\\s+(GTDB-Tk v.+)\\s+.+/\\1/' > versions.txt
        printf "### Set up completed successfully ###\\n\\n"
        """
}

// The processes below download the databases required by humann3.
process SETUP_CHOCOPHLAN {

    tag "Downloading and setting up Humann-s chocoplan nucleotide database..."
    label "humann_setup"
    label "db_setup"

    output:
        path("humann3-db/chocophlan"), emit: chocophlan_dir
        path("humann3-db/CHOCOPHLAN_DB_SETUP"), emit: completion_indicator
        path("versions.txt"), emit: version
    script:
        """
        [ -d humann3-db/ ] || mkdir -p humann3-db/
        printf "### Setting up humann3 reference databases ###\\n\\n"

        if [ ! -f humann3-db/CHOCOPHLAN_DB_SETUP ]; then
            printf "  Downloading full chocophlan db:\\n\\n" 
            # No need to update locations since I pass them as arguaments to the script
            humann3_databases --update-config no --download chocophlan full humann3-db/ && \\
            touch humann3-db/CHOCOPHLAN_DB_SETUP
            humann3 --version  > versions.txt
            printf "### Set up completed successfully ###\\n\\n"
        fi
        
        """
}


process SETUP_UNIREF {

    tag "Downloading and setting up Humann-s uniref protein database..."
    label "humann_setup"
    label "db_setup"

    output:
        path("humann3-db/uniref/"), emit: uniref_dir
        path("humann3-db/UNIREF_DB_SETUP"), emit: completion_indicator
        path("versions.txt"), emit: version
    script:
        """
        [ -d humann3-db/ ] || mkdir -p humann3-db/
        printf "### Setting up humann3's uniref database ###\\n\\n"
        if [ ! -f humann3-db/UNIREF_DB_SETUP ];then
            printf "\\n\\n  Downloading uniref90_ec_filtered_diamond db:\\n\\n"
            # No need to update locations since I pass them as arguaments to the script
            humann3_databases  --update-config no --download uniref uniref90_ec_filtered_diamond humann3-db/ && \\
            touch humann3-db/UNIREF_DB_SETUP
            humann3 --version  > versions.txt
            printf "### Set up completed successfully ###\\n\\n"
        fi
        """
}

process SETUP_UTILITY_MAPPING {

    tag "Downloading and setting up Humann-s utilities mapping database..."
    label "humann_setup"
    label "db_setup"

    output:
        path("humann3-db/utility_mapping/"), emit: utilities_dir
        path("humann3-db/UTILITY_MAPPING_SETUP"), emit: completion_indicator
        path("versions.txt"), emit: version
    script:
        """
        [ -d humann3-db/ ] || mkdir -p humann3-db/
        printf "### Setting up humann3's utilities mapping database ###\\n\\n"
        if [ ! -f humann3-db/UTILITY_MAPPING_SETUP ];then
            printf "\\n\\n  Downloading full utility_mapping db:\\n\\n"

            # Containers are read only but conda environments can be modified
            if [ ${params.use_conda} == 'true'];then

               humann3_databases --update-config yes --download utility_mapping full humann3-db/ && \\
               touch humann3-db/UTILITY_MAPPING_SETUP

             else 
           
               humann3_databases --update-config no --download utility_mapping full humann3-db/ && \\
               touch humann3-db/UTILITY_MAPPING_SETUP

            fi
            humann3 --version  > versions.txt
            printf "### Set up completed successfully ###\\n\\n"
        fi
        """
}


process SETUP_METAPHLAN {
    tag "Downloading and setting up Humann-s utilities mapping database..."
    label "humann_setup"
    label "db_setup"

    output:
        path("metaphlan4-db/"), emit: metaphlan_db_dir
        path("metaphlan4-db/METAPHLAN4_DB_SETUP"), emit: completion_indicator
        path("versions.txt"), emit: version
    script:
        """
        [ -d metaphlan4-db/ ] || mkdir -p metaphlan4-db/
        printf "### Setting up metaphlan's reference database ###\\n\\n"

        if [ ! -f metaphlan4-db/METAPHLAN4_DB_SETUP ];then
            printf "\\n\\n  Downloading metaphlan db:\\n\\n"
            metaphlan --install --bowtie2db metaphlan4-db/ && \\
            touch metaphlan4-db/METAPHLAN4_DB_SETUP
            metaphlan --version > versions.txt
            printf "### Set up completed successfully ###\\n\\n"
        fi
        """
}


workflow make_humann_db {

    main:
        SETUP_CHOCOPHLAN()
        SETUP_UNIREF()
        SETUP_UTILITY_MAPPING()
        SETUP_METAPHLAN()

        software_versions_ch = Channel.empty()
        SETUP_CHOCOPHLAN.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SETUP_UNIREF.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SETUP_METAPHLAN.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SETUP_UTILITY_MAPPING.out.version | mix(software_versions_ch) | set{software_versions_ch}

    emit:
       chocophlan_dir  = SETUP_CHOCOPHLAN.out.chocophlan_dir
       uniref_dir = SETUP_UNIREF.out.uniref_dir
       metaphlan_db_dir = SETUP_METAPHLAN.out.metaphlan_db_dir
       utilities_dir = SETUP_UTILITY_MAPPING.out.utilities_dir
       versions = software_versions_ch

}

workflow make_databases {

    take:
        CAT_DB_LINK
        GTDBTK_URL

    main:
        SETUP_CAT_DB(CAT_DB_LINK)
        SETUP_KOFAMSCAN_DB()
        SETUP_GTDBTK_DB(GTDBTK_URL)
        make_humann_db()

        software_versions_ch = Channel.empty()
        SETUP_CAT_DB.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SETUP_KOFAMSCAN_DB.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SETUP_GTDBTK_DB.out.version | mix(software_versions_ch) | set{software_versions_ch}
        make_humann_db.out.versions | mix(software_versions_ch) | set{software_versions_ch}

        
   
    emit:
        cat_db = SETUP_CAT_DB.out.cat_db
        kofam_db = SETUP_KOFAMSCAN_DB.out.ko_db_dir
        gtdbtk_db_dir = SETUP_GTDBTK_DB.out.gtdbtk_db_dir
        chocophlan_dir  = make_humann_db.out.chocophlan_dir
        uniref_dir = make_humann_db.out.uniref_dir
        metaphlan_db_dir = make_humann_db.out.metaphlan_db_dir
        utilities_dir = make_humann_db.out.utilities_dir
        versions = software_versions_ch

    }



workflow {
     make_databases(Channel.of(params.CAT_DB_LINK), Channel.of(params.GTDBTK_LINK))
}
