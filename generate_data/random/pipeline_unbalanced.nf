params.bacteria    = 'ECOLI_2';
params.outDir      = "$baseDir/datasets/${params.bacteria1}/UP/random"
// params.indexTrain  = "$baseDir/train.py"
// params.region      = '35+s+10'

//Variables
// indexBedFile       = file("$baseDir/bacteria/${params.bacteria}/${params.bacteria}.bed")
// indexFastaFile     = file("$baseDir/bacteria/${params.bacteria}/${params.bacteria}.fasta")
// indexBedFile       = file("$baseDir/bacteria/${params.bacteria}/${params.bacteria}.bed")
scanFastaFile      = file("$baseDir/datasets/${params.bacteria1}/UP/positive.fasta")
promoterBedFile    = file("$baseDir/datasets/${params.bacteria1}/positive.bed")
indexFastaFile     = file("$baseDir/datasets/${params.bacteria1}/complete.fasta")
lengthFile         = file("${params.bacteria2}.genome")
// trainFile          = file(params.indexTrain)


// negativeFile       = file("${params.outDir}/negative_40.fasta")

// noidFile           = file("$baseDir/scan_data_new_1/${params.bacteria}/positive_no_id.fasta")
n = (Math.random()*(32-27+1))+27
// n = len_list

Channel
   .from(indexFastaFile)
   .splitFasta( record: [id: true, seqString: true ])
   .collectFile(name: lengthFile) { record -> record.id + "\t" + record.seqString.length() + "\n"}
   .set{indexGenomeFile}

// Channel
//    .from(indexFastaFile)
//    .splitFasta( record: [id: true, seqString: true ])
//    .collectFile(name: lengthFile) { record -> record.id + "\t" + record.seqString.length() + "\n"}
//    .set{indexGenomeFile2}

println """\
PROMOTER PREDICTION = ${params.bacteria}
===============================

AUTHOR: RUBEN CHEVEZ
DATE:   OCT 15, 2018

"""
.stripIndent()

println """\
${scanFastaFile}
===============================


"""
.stripIndent()

process getNegativePromotersBED {
    // container 'genomicpariscentre/bedtools'
    input:
        file promoterBedFile
        file scanFastaFile 
        file "${params.bacteria2}.genome" from indexGenomeFile
        publishDir params.outDir, mode: 'copy' 
    output:
        file "negative.bed" into negative_promoter_bed
    script: 
    """
    bedtools random -l 100 -n `grep -i -v ">" ${scanFastaFile} | echo \$(( \$(wc -l) * 10 )) ` -g ${params.bacteria2}.genome > negative.bed
    
    """
    // bedtools subtract -f 0.13 -s -A -a tmp_negative_promoter.bed -b ${promoterBedFile} > negative.bed//"tmp_negative_promoter.bed"
}

process getNegativePromoterFASTA {
    // container 'genomicpariscentre/bedtools'
    input:
        file negative_promoter_bed
        file indexFastaFile   
        publishDir params.outDir, mode: 'copy' 
    output:
        file "negative.fasta" into negative_promoter_sequence
    script:
    """
    bedtools getfasta -s  -bed ${negative_promoter_bed} -fi ${indexFastaFile}  -fo  "negative.fasta" 
    """
}