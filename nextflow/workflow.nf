#!/usr/bin/env nextflow
 
/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Check the hostnames against configured profiles

def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}


def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'SCAPE-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'SCAPE Workflow Summary'
    section_href: '"https://confluence.phibred.com/display/OMICS/Expression+Data+Resource+White+Paper'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

import java.text.SimpleDateFormat
Date now = new Date()
SimpleDateFormat timestamp_fmt = new SimpleDateFormat("yyyyMMdd_HHmmss")
def timestamp = timestamp_fmt.format(now)

def scape_username = System.getProperty("user.name")
//String hostname = java.net.InetAddress.getLocalHost().getHostName();

my_LIMS_url_URL = 'https://url.research.phibred.com/v4/api/requests/' + params.LIMS_ID


// FOR DEBUGGING
// println timestamp
// println scape_username
// println hostname
// println "${my_LIMS_url_URL}"

my_s3_input_fastq_timestamp = params.s3_input_fastq_timestamp

my_trimQ = params.trimQ
my_minLength = params.minLength

my_S3_Out_URL = params.outdir + "/"
my_S3_Out_URL = my_S3_Out_URL.replaceFirst(/s3:\/\//, "https://s3.console.aws.amazon.com/s3/buckets/")
my_S3_adfs = "https://sso.agcompany.net/adfs/ls/IdpInitiatedSignOn.aspx?loginToRp=urn:amazon:webservices"


if (!"${my_trimQ}".isInteger()) {
   my_trimQ = 20
}

if (!"${my_minLength}".isInteger()) {
   my_minLength = 35
}

// Start the Process execution

// Step 0 - Parse software version numbers

process get_software_versions {
    echo true
    tag "Get_Software_Versions"
    validExitStatus 0,1,127
    publishDir "${params.outdir}/0_software_versions/", mode: 'copy', pattern: '*.xt'

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file '*.txt' into software_versions_text

    script:
    """
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    bbversion.sh --version > v_bbduk.txt
    fastp --version >&  v_fastp.txt
    salmon --version > v_salmon.txt
    multiqc --version > v_multiqc.txt

    for X in `ls *.txt`; do
        cat \$X >> all_versions.txt;
    done
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

// Step 1.1 Sample_ids_from_url

if (params.include_samples_file) {
   opt_file = file(params.include_samples_file)

   // WE GET Subset of Samples
   
   process get_SUBSET_samples_from_url {

    publishDir "${params.outdir}/1_samples_url", mode:'copy', pattern:'*.txt'
    tag "Get_Subset_Sample_IDs"

    input:
    file opt from opt_file

    output:
    path "url.txt" into records
    file ("*.txt")

    script:
    
    """
    cat ${opt} | perl -pe 's#\015##g' > temp_list.txt
    curl --silent  -X GET ${my_LIMS_url_URL} -H "accept:application/json" | jq '.request.samples[].sample_id' | perl -pe 's#\"##g' > list_all_samples.txt

    cat list_all_samples.txt temp_list.txt | sort | uniq -d > wanted_list.txt

    if [ -s wanted_list.txt ]; then    
        cat wanted_list.txt | perl -pe 's#^#https://gsd.research.phibred.com/api/submit/#g'  > url.txt
    fi

    """
    } // end of process get_subset_samples_from_url

} else {

  // WE GET ALL Samples from url for a LIMS_ID

   process get_ALL_samples_from_s3 {

    publishDir "${params.outdir}/1_samples", mode:'copy', pattern:'*.txt'
    tag "Get_ALL_Sample_IDs"

    output:
    path "url.txt" into records

    script:
    """
    curl --silent  -X GET ${URL} -H "accept:application/json" | jq '.request.samples[].sample_id' | perl -pe 's#\"##g' | perl -pe 's#^#https://gsd.research.phibred.com/api/submit/#g'  > url.txt
    """
    } // end process get_ALL_samples_from_url
}

// Step - 1.2 Extract 

process extract {
    publishDir "${params.outdir}/2_CURL_url", mode:'copy'
    tag "CURL_SAMPLES_url"

    input:
    path "url.txt" from records

    output:
    path "s3_urls.txt" into curl_records_1
    path "job_urls.txt" into curl_records_2
    path "*.txt"

    script:
    """
    cp url.txt paths.txt
    for path in \$(cat paths.txt); do
	    curl --connect-timeout 18000  --retry 5  --silent -X POST \${path} 2>/dev/null | jq '.results_path' | perl -pe 's#\"\$#/#g' | perl -pe 's#\"##g' 
    done | tee s3_urls.txt

    cat s3_urls.txt | perl -pe 's#.*/(\\S+)/#https://url/here' > job_urls.txt
    """
}


// Process - 1.3 Check Data download from url job status


process check_jobs {
    publishDir "${params.outdir}/3_CURL_url_jobs", mode:'copy'
    tag "SLEEPING_FOR_url_DUMP"

    errorStrategy 'retry'
    maxRetries 5

    input:
    path "s3_urls.txt" from curl_records_1
    path "job_urls.txt" from curl_records_2

    output:
    path "*.txt"
    path pass_s3_urls into get_S3_aws_url

    script:
    """
    aws_url_job_stats.sh job_urls.txt > stdout_job_urls.txt
    cp  s3_urls.txt  pass_s3_urls
    """
}

// Step - 1.4 Get the S3 URLs

process get_S3_urls {
    publishDir "${params.outdir}/4_get_s3_results", mode:'copy'
    tag "Get_S3_URLS_url"

    input:
    path pass_s3_urls from get_S3_aws_url

    output:
    
    path ("*.txt") into download_url
    path processed_urls

    script:
    """
    create_URL_files_per_line.pl pass_s3_urls
    cp  pass_s3_urls  processed_urls
    """
}

// Step - 1.5 Dump data from url bucket to Scape bucket

process dump_S3 {
    label 'low_cm'
    publishDir "${params.outdir}/5_get_s3_dump", mode:'copy', pattern: "*S3_Copy_stats.txt"

    tag "S3_COPY_from_url"
    errorStrategy 'retry'
    maxRetries 5

    input:
    path s3_url from download_url.flatten()

    output:
    path "*S3_Copy_stats.txt" into concat_S3_Paths

    script:
    if (params.singleEnd) {
    
    	"""
	for path in \$(cat ${s3_url}); do
	aws s3 cp \${path} s3://input_data/${scape_username}/${my_s3_input_fastq_timestamp}/ --exclude "*" --include "*single*" --recursive | tail -n 1  |  perl -pe 's#.*s3.*/(\\S+)-.*gz.* to (s3\\S+)#\$1\t\$2#g' >& ${s3_url.simpleName}_S3_Copy_stats.txt
	done
	head -n 1 ${s3_url} >&  ${s3_url.simpleName}_S3_head_stats.txt
        """
    } else { // we have paired end now
 
    	"""
  	# head -n 1 ${s3_url} | xargs -I % aws s3 cp % s3://pd-nfs/scape/WORKSPACE/input_data/${scape_username}/${my_s3_input_fastq_timestamp}/ --exclude "*" --include "*pairs*" --recursive | tail -n 1  |  perl -pe 's#.*s3.*/(\\S+)-.*gz.* to (s3.*/).*gz#\$1\t\$2\$1-pairsR1.fq.gz\t\$2\$1-pairsR2.fq.gz#g' >& ${s3_url.simpleName}_S3_Copy_stats.txt
	
	for path in \$(cat ${s3_url}); do
	aws s3 cp \${path} s3://pd-nfs/scape/WORKSPACE/input_data/${scape_username}/${my_s3_input_fastq_timestamp}/  --exclude "*" --include "*pairs*" --recursive | tail -n 1  |  perl -pe 's#.*s3.*/(\\S+)-.*gz.* to (s3.*/).*gz#\$1\t\$2\$1-pairsR1.fq.gz\t\$2\$1-pairsR2.fq.gz#g' >& ${s3_url.simpleName}_S3_Copy_stats.txt
	done
	head -n 1 ${s3_url} >&  ${s3_url.simpleName}_S3_head_stats.txt
        """
    }
}

if (params.singleEnd) {
concat_S3_Paths
	.splitCsv(sep:'\t')
	.map{ row-> tuple(row[0], row[1]) }
	.into { raw_reads_fastqc; raw_reads_remove_contamination }
} else {
concat_S3_Paths
	.splitCsv(sep:'\t')
	.map{ row-> tuple(row[0], [ row[1], row[2] ] ) }
	.into { raw_reads_fastqc; raw_reads_remove_contamination }
}

// Step - 2 Fastqc_Before

process fastqc_before {
	  label 'mid_cm'

	  tag "${name}_FASTQC_BEFORE"
    
	  publishDir "${params.outdir}/FASTQC_Before", mode: 'copy',
               saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

	  input:
	  tuple val(name), path(sample) from raw_reads_fastqc

	  output:
	  path "*_fastqc.{zip,html}" into  preqc_results

          script:
	  """
	  fastqc -t 16 -q $sample
	  """
}

// STEP 3 - Remove Contaminants

process remove_contaminants {
    label 'low_cm'
    validExitStatus 0,1

    tag "${name}_BBMAP"
    cpus 16
    time '16h'
    memory '30 GB'
    publishDir "${params.outdir}/qc/remove_contaminants/bbmap/stats/",  mode: 'copy', pattern: "*.txt"

    input:
    tuple val(name), path(reads) from raw_reads_remove_contamination
    path 'contaminant.fa' from params.contaminant_fasta

    output:
    tuple val(name), path ("*.clean.fastq.gz") into clean_reads
    path "*.txt" into contaminant_stats

    script:
    if (params.singleEnd){
       if(params.contaminant_software == 'bbmap') {

       """
                bbmap.sh  \\
                  in=${name}-singles.fq.gz \\
                  outu=${name}.clean.fastq.gz \\
                  ref=contaminant.fa \\
                  maxindel=1 minid=0.95 nodisk \\
                  statsfile=${name}_bbmap.trimstats.txt >& ${name}_bbmap.stdout_stats.txt
       """

       }
    }else { // we have paired end now
       if(params.contaminant_software == 'bbmap') {

       """
        bbmap.sh  \\
                  in1=${name}-pairsR1.fq.gz  in2=${name}-pairsR2.fq.gz \\
                  outu1=${name}_1.clean.fastq.gz  outu2=${name}_2.clean.fastq.gz \\
                  ref=contaminant.fa \\
                  maxindel=1 minid=0.95 nodisk \\
                  statsfile=${name}_bbmap.trimstats.txt >& ${name}_bbmap.stdout_stats.txt
       """
       }
    }
} // end of Step 2 - Process remove_contaminants - channel(s) output clean_reads


// STEP 4 - Remove Adapters

if(params.adapter_software == 'bbduk') {
     process remove_adapters_bbduk {
         label 'low_cm'
         validExitStatus 0,1

         tag "${name}_BBDUK"
         cpus 16
         time '16h'
         memory '30 GB'
	
        publishDir "${params.outdir}/qc/adapter_quality_trimming/bbduk", mode: 'copy', pattern: "*stats.txt"

        input:
        tuple val(name), path(clean_read) from clean_reads
	path 'adapter.fa' from params.adapter_fasta

        output:
        tuple val(name), path ("*_trimmed_clean.fastq.gz") into analysis_ready_reads, post_trimming_qc
        path "*stats.txt" into adapter_removal_stats

        script:
        if (params.singleEnd){
           """
           bbduk.sh in=${name}.clean.fastq.gz  out=${name}_trimmed_clean.fastq.gz \\
                  ref=adapter.fa \\
                  minlength=${my_minLength}  trimq=${my_trimQ}   ktrim=r qtrim=rl  k=23   mink=11 hdist=1 tbo \\
		  stats=${name}_bbduk.trimstats.txt   refstats=${name}_bbduk.refstats.txt >& ${name}_bbduk.stdout_stats.txt	
           """
        }else { // we have paired end now
           """
           bbduk.sh in1=${name}_1.clean.fastq.gz  in2=${name}_2.clean.fastq.gz  \\
                  out1=${name}_1_trimmed_clean.fastq.gz  out2=${name}_2_trimmed_clean.fastq.gz \\
                  ref=adapter.fa \\
                  minlength=${my_minLength}  trimq=${my_trimQ} ktrim=r qtrim=rl  k=23 mink=11 hdist=1 tpe tbo \\
                  stats=${name}_bbduk.trimstats.txt   refstats=${name}_bbduk.refstats.txt >& ${name}_bbduk.stdout_stats.txt
           """
           } // end of else section for if SE or PE
    } // end of bbduk block
}

if(params.adapter_software == 'fastp') {
     process remove_adapters_fastp {
         label 'low_cm'
         validExitStatus 0,1
         tag "${name}_FASTP"
         cpus 16
         time '16h'
         memory '30 GB'

        publishDir "${params.outdir}/qc/adapter_quality_trimming/fastp", mode: 'copy', pattern: "*.{json,html}"

        input:
        tuple val(name), path(clean_read) from clean_reads

        output:
        tuple val(name), path ("*_trimmed_clean.fastq.gz") into analysis_ready_reads, post_trimming_qc
        path "*.{json,html}" into adapter_removal_stats

        script:
        if (params.singleEnd){
             """
             fastp --in1=${name}.clean.fastq.gz \\
             --out1=${name}_trimmed_clean.fastq.gz \\
             --qualified_quality_phred=${my_trimQ} --length_required=${my_minLength}

	     mv fastp.html ${name}_fastp.html
	     mv fastp.json ${name}_fastp.json
             """
         } else { // we have paired end now
            """
            fastp --in1=${name}_1.clean.fastq.gz  --in2=${name}_2.clean.fastq.gz  \\
              --out1=${name}_1_trimmed_clean.fastq.gz  --out2=${name}_2_trimmed_clean.fastq.gz \\
              --qualified_quality_phred=${my_trimQ} --length_required=${my_minLength}

	    mv fastp.html ${name}_fastp.html
	    mv fastp.json ${name}_fastp.json
            """
       } //  end of else section for if SE or PE
    } // end of fastp block
}// end of Step 3 - Process remove_adapters - channel(s) output analysis_ready_reads, post_trimming_qc


// STEP 5a - FastQC After

process fastqc_after {
    label 'mid_cm'

    tag "${name}_FASTQC_AFTER"

    publishDir "${params.outdir}/FASTQC_After", mode: 'copy',
              saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    tuple val(name), path(reads) from post_trimming_qc

    output:
    path "*_fastqc.{zip,html}" into fastqc_after_results

    script:
    """
    fastqc -t 16 -q $reads
    """
}

// STEP 5b - Salmon quant

my_salmon_libtype = params.salmon_libtype
my_quant_method = params.quant_method

process quant {
    label 'mid_cm'

    tag "${pair_id}_SALMON"
    time '20h'
    publishDir "${params.outdir}/quant", mode:'copy'

    input:
    path 'my_salmon_index' from params.salmon_index_path
    tuple pair_id, path(reads) from analysis_ready_reads

    output:
    path "salmon_${pair_id}" into quant_ch, quant_ch_summary

    script:
    if (params.singleEnd){
       if(my_quant_method == 'salmon') {

       """
                salmon --no-version-check  quant --threads 16 --seqBias --validateMappings --numBootstraps 100 -l $my_salmon_libtype  \\
               --writeUnmappedNames -i my_salmon_index   \\
               -r ${pair_id}_trimmed_clean.fastq.gz                  \\
               -o salmon_${pair_id}

       """

       }
    }else { // we have paired end now

     if(my_quant_method == 'salmon') {

       """
        salmon --no-version-check quant --threads 16 --seqBias --validateMappings --numBootstraps 100 -l $my_salmon_libtype \\
               --writeUnmappedNames -i my_salmon_index   \\
               -1 ${pair_id}_1_trimmed_clean.fastq.gz  -2 ${pair_id}_2_trimmed_clean.fastq.gz   \\
               -o salmon_${pair_id}

       """
       }
      }
}

/*
 * STEP 6 - MultiQC
*/

Channel.fromPath(params.multiqc_config, checkIfExists: true).set { ch_config_for_multiqc }

process multiQC {
    label 'high_cm'
    tag "MultiQC"

    validExitStatus 0,1,143
    errorStrategy 'ignore'
    publishDir "${params.outdir}/multiqc/", mode: 'copy'

    input:
    file multiqc_config from ch_config_for_multiqc
    file ('pre_fastqc/*') from  preqc_results.collect().ifEmpty([])
    file ('bbmap_stats/*') from contaminant_stats.collect().ifEmpty([])
    file ('adapter_removal/*') from adapter_removal_stats.collect().ifEmpty([])
    file ('fastqc_after/*') from fastqc_after_results.collect().ifEmpty([])
    file ('quant/*') from quant_ch.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml

    output:
    file 'multiqc_report.html' into multiqc_report


    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + "multiqc_report" : ''

    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config --interactive
    """
}

// Step 7 - TxImport

if(params.tx2gene_file) {

     process txImport_summary_with_tx2gene {
        label 'mid_cm'

	tag "TxImport"

        publishDir "${params.outdir}/summary", mode: 'copy'

        input:
        path 'quant/*' from quant_ch_summary.collect()
	path 'my_tx2gene' from params.tx2gene_file

        output:
        path "*.csv"
	val true into done_ch

        script:
        """
        forAWS_getSalmonQuants_txImport.r  quant my_tx2gene
        """
        }
}else {
        process txImport_summary_no_tx2gene {
        label 'mid_cm'
	tag "TxImport_No_tx2gene"

        publishDir "${params.outdir}/summary", mode: 'copy'

        input:
        path 'quant/*' from quant_ch_summary.collect()

        output:
        path "*.csv"
	val true into done_ch	

        script:
        """
        forAWS_getSalmonQuants_txImport.r  quant
        """
    }
}

// Step 8 - Deleting S3 Input FQ files

process delete_S3_fq {

    tag "Delete_Input_FQ"
    //label 'low_cm'

    publishDir "${params.outdir}/deleted_S3_input_fastq", mode: 'copy', pattern: "*delete_S3_objects.txt"

    input:
    val flag from done_ch

    output:
    path "*delete_S3_objects.txt"

    script:
    """
    aws s3 rm s3://input_data/${scape_username}/${my_s3_input_fastq_timestamp}/  --recursive >& ${timestamp}_delete_S3_objects.txt
    """
}


// Completion e-mail notification

// admin Email list
params.admin="pranay.pashine@corteva.com,anand.venkatraman@corteva.com"

//Send email Notification
workflow.onComplete {

    def msg = """
--------------------------------------------------------------------------------------------------------------------------

Please login first using the Corteva AWS ADFS Sign On - ${my_S3_adfs}

Check your PIPELINE Output at ${my_S3_Out_URL} 

To download from S3 to On-Prem, do: aws s3 cp ${params.outdir} <your_onprem_location> --recursive

--------------------------------------------------------------------------------------------------------------------------

The command used to launch the workflow was as follows:

${workflow.commandLine}

====================
   Execution Summary
====================

Launch time		${workflow.start}

Completed at		${workflow.complete}

Total CPU-Hours	${workflow.stats.computeTimeFmt ?: '-'}

Tasks stats		Success ${workflow.stats.succeedCountFmt}; Cached ${workflow.stats.cachedCountFmt}; Ignored ${workflow.stats.ignoredCountFmt}; Failed ${workflow.stats.failedCountFmt}

Duration		${workflow.duration}

projectDir		${workflow.projectDir}

Script File Path	${workflow.scriptFile}

Script Name		${workflow.scriptName}

Workflow Profile	${workflow.profile}

Container Image	${workflow.container}

Script ID		${workflow.scriptId}

runName		${workflow.runName}

Nextflow Version	${workflow.nextflow.version}

outDir			${params.outdir}

Nextflow Compile 	${workflow.nextflow.timestamp}

SCAPE Version		${params.scapeVersion}

SCAPE BuildDate	${params.scapeBuildDate}

v_FastQC		${params.v_FastQC}

v_BBTools		${params.v_BBTools}

v_FASTP		${params.v_FASTP}

v_MultiQC		${params.v_MultiQC}

salmon_Index_Version	${params.salmon_Index_Version}

salmon_Quant_Version    ${params.Salmon_Quant_Version}

--------------------------------------------------------------------------------------------------------------------------
        """
    def status = "NA"
    def local_report = file('multiqc_report.html')
    multiqc_report.getVal().copyTo( local_report )
    
    if(workflow.success) {
        status = "SUCCESS"
        sendMail {
            to "${params.email}"
            bcc "${params.admin}"
            subject "${status} - ${params.LIMS_ID} - SCAPE Scheduler Prod - AWS url Workflow Completion"
            attach   local_report
            body
            """
Execution Completed Successfully.

            ${msg}
            """
            .stripIndent()
                }
        }
        else {
             status = "FAILED"

        sendMail {
            to "${params.email}"
            bcc "${params.admin}"
            subject "${status} - ${params.LIMS_ID} - SCAPE from AWS url - Workflow Completion"
            attach "${workflow.launchDir}/.nextflow.log"
            body
            """
Execution Completed Unsuccessfully, check the attached log file for full information on error  -> ${workflow.errorReport}

            ${msg}
            """
            .stripIndent()
                }
        }
}

