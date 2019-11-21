#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
========================================================================================
                        ICGC-ARGO DNA SEQ ALIGNMENT
========================================================================================
#### Homepage / Documentation
https://github.com/icgc-argo/nextflow-dna-seq-alignment
#### Authors
Alexandru Lepsa @lepsalex <alepsa@oicr.on.ca>
----------------------------------------------------------------------------------------

#### WORK IN PROGRESS
TODO: Replace A1/A2 with actual analysis types
TODO: Docker container version to latest?

Required Parameters (no default):
--study_id                              song study ID
--analysis_id                           song A1 analysis ID
--song_url                              song server URL
--score_url                             score server URL
--api_token                             song/score API Token
--aligned_basename                      final aligned filename

General Parameters (with defaults):
--reference_dir                         reference genome directory
--aligned_lane_prefix                   prefix for alignment (defaults to "grch38-aligned")
--cpus                                  cpus given to all process containers (default 1)
--memory                                memory (MB) given to all process containers (default 1024)

Download Parameters (object):
--download
{
    song_container_version              song docker container version, defaults set below
    score_container_version             score docker container version, defaults set below
    song_url                            song url for download process (defaults to main song_url param)
    score_url                           score url for download process (defaults to main score_url param)
    api_token                           song/score API token for download process (defaults to main api_token param)
    cpus                                cpus for process container, defaults to cpus parameter
    memory                              memory (MB) for process container, defaults to memory parameter
}

Preprocess Parameters (object):
---preprocess
{
    container_version                   docker container version, defaults set below
    reads_max_discard_fraction          preprocess reads max discard function
    cpus                                cpus for preprocess container, defaults to cpus parameter
    memory                              memory (MB) for preprocess container, defaults to memory parameter
}

Align Parameters (object):
--align
{
    container_version                   docker container version, defaults set below
    cpus                                cpus for align container, defaults to cpus parameter
    memory                              memory (MB) for align container, defaults to memory parameter
}

Merge Parameters (object):
--merge
{
    container_version                   docker container version, defaults set below
    output_format                       options are ['cram', 'bam']
    markdup                             TODO: write description
    lossy                               TODO: write description
    cpus                                cpus for merge container, defaults to cpus parameter
    memory                              memory (MB) for merge container, defaults to memory parameter
}

Upload Parameters (object):
--upload
{
    song_container_version              song docker container version, defaults set below
    score_container_version             score docker container version, defaults set below
    song_url                            song url for upload process (defaults to main song_url param)
    score_url                           score url for upload process (defaults to main score_url param)
    api_token                           song/score API token for upload process (defaults to main api_token param)
    cpus                                cpus for process container, defaults to cpus parameter
    memory                              memory (MB) for process container, defaults to memory parameter
}

*/

params.reference_dir = 'reference'
params.aligned_lane_prefix = 'grch38-aligned'
params.cpus = 1
params.memory = 1024

params.download = [
    'song_container_version': 'latest',
    'score_container_version': 'latest',
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    'cpus': params.cpus,
    'mem': params.memory,
    *:(params.download ?: [:])
]

params.preprocess = [
    'container_version': '0.1.5.0',
    'reads_max_discard_fraction': 0.05,
    'cpus': params.cpus,
    'mem': params.memory,
    *:(params.preprocess ?: [:])
]

params.align = [
    'container_version': '0.1.2',
    'cpus': params.cpus,
    'mem': params.memory,
    *:(params.align ?: [:])
]

params.merge = [
    'container_version': '0.1.4',
    'output_format': ['cram'],
    'markdup': 'OPTIONAL_INPUT',
    'lossy': 'OPTIONAL_INPUT',
    'cpus': params.cpus,
    'mem': params.memory,
    *:(params.align ?: [:])
]

// params.upload.container_version = 'latest'
// params.upload.song_url = parms.song_url
// params.upload.score_url = parms.score_url
// params.upload.api_token = params.api_token
// params.upload.cpus = params.cpu
// params.upload.memory = params.memory

// Include all modules and pass params
include song_score_download as download from './data-processing/modules/song_score_download' params(params = params.download)                                                                             
include preprocess from './dna-seq-processing/modules/seq_data_to_lane_bam' params(params.preprocess)
include bwaMemAligner as align from './dna-seq-processing/modules/bwa_mem_aligner.nf' params(params.align)
include bamMergeSortMarkdup as merge from './dna-seq-processing/modules/bam_merge_sort_markdup.nf' params(params.merge)
// include songScoreUpload as upload from './data-processing/modules/song_score_upload' params(params.upload)

ref_gnome = Channel.fromPath("${params.reference_dir}/*").collect()

workflow {
    // download files and metadata from song/score (A1)
    download(params.study_id, params.analysis_id)

    // run files through preprocess step (split to lanes)
    preprocess(download.out.analysis_json_and_files)

    // align each lane independently
    align(preprocess.out.lane_bams, ref_gnome, params.aligned_lane_prefix)

    // // collect aligned lanes for merge and markdup
    merge(align.out.aligned_file.collect(), ref_gnome, params.aligned_basename)

    // // upload aligned file and metadata to song/score (A2)
    // download.out.analysis_json.view()
    // upload(merge.out)
}
