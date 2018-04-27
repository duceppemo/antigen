#!/bin/bash

baseDir='/home/bioinfo/Desktop/tags3'


cpu=$(nproc) #total number of cores
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
memJava="-Xmx"$mem"g"
prog=""${HOME}"/prog"

# Download reference genomes
# NC_002945

# Get fasta file
esearch -db nuccore -query "NC_002945 [ACCN]" \
    | efetch -format fasta \
    > "${baseDir}"/NC_002945.fasta

#Get full genbank file
esearch -db nuccore -query "NC_002945 [ACCN]" \
    | efetch -format gbwithparts \
    > "${baseDir}"/NC_002945.gbk


# Extract each locustags in a separate fasta file
python << END

import os
from Bio import SeqIO

# genbank_file = '$baseDir/NC_002945.gbk'
genbank_file = '/home/bioinfo/Desktop/nadia/ORFs/annotation/MBWGS010.gbk'
out_folder_faa = '$baseDir/faa'
out_folder_ffn = '$baseDir/ffn'

# Create output folders for faa and fnn files
if not os.path.exists(out_folder_faa):
    os.makedirs(out_folder_faa)

if not os.path.exists(out_folder_ffn):
    os.makedirs(out_folder_ffn)

'''
locus_tag_list = ['BQ2027_MB0296', 'BQ2027_MB1662', 'BQ2027_MB3645C', 'BQ2027_MB3247C',
                  'BQ2027_MB2898', 'BQ2027_MB1918C', 'BQ2027_MB1230', 'BQ2027_MB2900',
                  'BQ2027_MB3834C', 'BQ2027_MB0134C', 'BQ2027_MB2057C', 'BQ2027_MB0293',
                  'BQ2027_MB1227', 'BQ2027_MB2903C', 'BQ2027_MB0130', 'BQ2027_MB0939C',
                  'BQ2027_MB0295', 'BQ2027_MB3766C', 'BQ2027_MB1207C', 'BQ2027_MB1967',
                  'BQ2027_MB3904', 'BQ2027_MB3905', 'BQ2027_MB0041C', 'BQ2027_MB0064',
                  'BQ2027_MB0129', 'BQ2027_MB0184', 'BQ2027_MB0285C', 'BQ2027_MB0286C',
                  'BQ2027_MB0287C', 'BQ2027_MB0305', 'BQ2027_MB0368', 'BQ2027_MB0426',
                  'BQ2027_MB0545', 'BQ2027_MB0593C', 'BQ2027_MB0690', 'BQ2027_MB0768',
                  'BQ2027_MB0777C', 'BQ2027_MB0845C', 'BQ2027_MB0856', 'BQ2027_MB0858',
                  'BQ2027_MB0896C', 'BQ2027_MB1002', 'BQ2027_MB1006C', 'BQ2027_MB1018C',
                  'BQ2027_MB1044C', 'BQ2027_MB1050', 'BQ2027_MB1096C', 'BQ2027_MB1097C',
                  'BQ2027_MB1116', 'BQ2027_MB1118', 'BQ2027_MB1121', 'BQ2027_MB1262C',
                  'BQ2027_MB1275C', 'BQ2027_MB1360C', 'BQ2027_MB1476C', 'BQ2027_MB1485C',
                  'BQ2027_MB1487C', 'BQ2027_MB1575C', 'BQ2027_MB1679C', 'BQ2027_MB1783C',
                  'BQ2027_MB1784C', 'BQ2027_MB1788', 'BQ2027_MB1831C', 'BQ2027_MB2006C',
                  'BQ2027_MB2108', 'BQ2027_MB2125C', 'BQ2027_MB2323', 'BQ2027_MB2410C',
                  'BQ2027_MB2506C', 'BQ2027_MB2517C', 'BQ2027_MB2606C', 'BQ2027_MB2667C',
                  'BQ2027_MB2734', 'BQ2027_MB2740C', 'BQ2027_MB2878', 'BQ2027_MB3183C',
                  'BQ2027_MB3402', 'BQ2027_MB3453C', 'BQ2027_MB3469C', 'BQ2027_MB3481',
                  'BQ2027_MB3537', 'BQ2027_MB3538', 'BQ2027_MB3541', 'BQ2027_MB3543',
                  'BQ2027_MB3562', 'BQ2027_MB3751', 'BQ2027_MB3761C', 'BQ2027_MB3792',
                  'BQ2027_MB3833C', 'BQ2027_MB3841', 'BQ2027_MB3909C', 'BQ2027_MB3921C',
                  'BQ2027_MB0419C', 'BQ2027_MB0975C', 'BQ2027_MB1397C', 'BQ2027_MB1530',
                  'BQ2027_MB1724', 'BQ2027_MB2377C', 'BQ2027_MB2554C']
'''
'''
locus_tag_list = ['BQ2027_MB0248C', 'BQ2027_MB0249', 'BQ2027_MB0254C', 'BQ2027_MB0276',
                  'BQ2027_MB0340', 'BQ2027_MB0439', 'BQ2027_MB0471', 'BQ2027_MB0476',
                  'BQ2027_MB0478', 'BQ2027_MB0485', 'BQ2027_MB0515C', 'BQ2027_MB0584',
                  'BQ2027_MB0592', 'BQ2027_MB0598C', 'BQ2027_MB0649C', 'BQ2027_MB0658',
                  'BQ2027_MB0723', 'BQ2027_MB0727', 'BQ2027_MB0740', 'BQ2027_MB0754',
                  'BQ2027_MB0784C', 'BQ2027_MB0824', 'BQ2027_MB0854C', 'BQ2027_MB0908C',
                  'BQ2027_MB0920', 'BQ2027_MB0929', 'BQ2027_MB0956C', 'BQ2027_MB1010',
                  'BQ2027_MB1099C', 'BQ2027_MB1109C', 'BQ2027_MB1123', 'BQ2027_MB1164C',
                  'BQ2027_MB1213', 'BQ2027_MB1301C', 'BQ2027_MB1387', 'BQ2027_MB1427',
                  'BQ2027_MB1656', 'BQ2027_MB1857', 'BQ2027_MB1951C', 'BQ2027_MB1950',
                  'BQ2027_MB2164C', 'BQ2027_MB2238', 'BQ2027_MB2244', 'BQ2027_MB2265',
                  'BQ2027_MB2270', 'BQ2027_MB2397C', 'BQ2027_MB2587', 'BQ2027_MB2624C',
                  'BQ2027_MB2656', 'BQ2027_MB2659C', 'BQ2027_MB2801C', 'BQ2027_MB2810',
                  'BQ2027_MB2913C', 'BQ2027_MB2965C', 'BQ2027_MB2970C', 'BQ2027_MB3026C',
                  'BQ2027_MB3071', 'BQ2027_MB3074C', 'BQ2027_MB3219C', 'BQ2027_MB3389',
                  'BQ2027_MB3424C', 'BQ2027_MB3451C', 'BQ2027_MB3473C', 'BQ2027_MB3486C',
                  'BQ2027_MB3487C', 'BQ2027_MB3687C', 'BQ2027_MB3743C', 'BQ2027_MB3749C',
                  'BQ2027_MB3826', 'BQ2027_MB3871', 'BQ2027_MB3876']
'''
'''
locus_tag_list = ['BQ2027_MB1228', 'BQ2027_MB1396C', 'BQ2027_MB1782C', 'BQ2027_MB1789C',
                  'BQ2027_MB1836', 'BQ2027_MB3626C', 'BQ2027_MB0014C', 'BQ2027_MB0007',
                  'BQ2027_MB0312C', 'BQ2027_MB1951C', 'BQ2027_MB2970C', 'BQ2027_MB3505',
                  'BQ2027_MB2648C', 'BQ2027_MB3184C', 'BQ2027_MB3707', 'BQ2027_MB3595']
'''
locus_tag_list = ['GNBNAFKM_00129', 'GNBNAFKM_00208', 'GNBNAFKM_00213', 'GNBNAFKM_00324',
                  'GNBNAFKM_00531', 'GNBNAFKM_00600', 'GNBNAFKM_00604', 'GNBNAFKM_01019',
                  'GNBNAFKM_01070', 'GNBNAFKM_01176', 'GNBNAFKM_01633', 'GNBNAFKM_01643',
                  'GNBNAFKM_01724', 'GNBNAFKM_01799', 'GNBNAFKM_01968', 'GNBNAFKM_01973',
                  'GNBNAFKM_02111', 'GNBNAFKM_02273', 'GNBNAFKM_02274', 'GNBNAFKM_02282',
                  'GNBNAFKM_02518', 'GNBNAFKM_02865', 'GNBNAFKM_03022', 'GNBNAFKM_03269',
                  'GNBNAFKM_03372', 'GNBNAFKM_03375', 'GNBNAFKM_03588', 'GNBNAFKM_03590',
                  'GNBNAFKM_03872', 'GNBNAFKM_03910', 'GNBNAFKM_03989', 'GNBNAFKM_04068',
                  'GNBNAFKM_01235', 'GNBNAFKM_01247', 'GNBNAFKM_00325', 'GNBNAFKM_03804',
                  'GNBNAFKM_02403', 'GNBNAFKM_04029', 'GNBNAFKM_02281', 'GNBNAFKM_02401',
                  'GNBNAFKM_03047', 'GNBNAFKM_02669', 'GNBNAFKM_01516', 'GNBNAFKM_00239']


suffix_faa = '.faa'
suffix_ffn = '.ffn'

out_all_faa = '$baseDir/NC_002945.faa'
out_all_file_faa = open(out_all_faa, 'w')

out_all_ffn = '$baseDir/NC_002945.ffn'
out_all_file_ffn = open(out_all_ffn, 'w')

# Get protein sequences of CDS
for record in SeqIO.parse(genbank_file, "genbank"):
    for f in record.features:
        if f.type == "CDS":
            # All protein
            out_all_file_faa.write('>' + f.qualifiers["locus_tag"][0] + "\n")
            out_all_file_faa.write(f.qualifiers["translation"][0] + "\n")

            # All coding sequences
            seq = str(f.extract(record.seq))
            out_all_file_ffn.write('>' + f.qualifiers["locus_tag"][0] + "\n")
            out_all_file_ffn.write(seq + "\n")

            if f.qualifiers["locus_tag"][0] in locus_tag_list:
                # Protein
                out_faa = os.path.join(out_folder_faa, f.qualifiers["locus_tag"][0] + suffix_faa)
                out_file_faa = open(out_faa, 'w')
                out_file_faa.write('>' + f.qualifiers["locus_tag"][0] + "\n")
                out_file_faa.write(f.qualifiers["translation"][0])
                out_file_faa.close()

                # DNA
                out_ffn = os.path.join(out_folder_ffn, f.qualifiers["locus_tag"][0] + suffix_ffn)
                out_file_ffn = open(out_ffn, 'w')
                seq = str(f.extract(record.seq))
                out_file_ffn.write('>' + f.qualifiers["locus_tag"][0] + "\n")
                out_file_ffn.write(seq)
                out_file_ffn.close()

out_all_file_faa.close()
out_all_file_ffn.close()

END

# Fix fasta to make sequence lines max of 80 characters
find "${baseDir}/faa" -type f -name "*.faa" \
    | parallel 'perl /home/bioinfo/scripts/formatFasta.pl -i {} -o {}.tmp -w 80; mv {}.tmp {}'
find "${baseDir}/ffn" -type f -name "*.ffn" \
    | parallel 'perl /home/bioinfo/scripts/formatFasta.pl -i {} -o {}.tmp -w 80; mv {}.tmp {}'

# Merge all locus_tags in one file
cat ${baseDir}/faa/*.faa > ${baseDir}/faa/all_locus_tags.faa
cat ${baseDir}/ffn/*.ffn > ${baseDir}/ffn/all_locus_tags.ffn


########## Trimming ##########


#trimm reads
bbduk.sh \
    "$memJava" \
    threads="$cpu" \
    in1='/media/6tb_raid10/data/Mycobaterium_bovis/canada/MBWGS010_R1.fastq.gz' \
    in2='/media/6tb_raid10/data/Mycobaterium_bovis/canada/MBWGS010_R2.fastq.gz' \
    ref="${prog}"/bbmap/resources/nextera.fa.gz \
    ktrim=r k=23 mink=11 hdist=1 tbo tpe \
    qtrim=lr trimq=10 \
    minlen=64 \
    out1='$baseDir/fastq/MBWGS010_Trimmed_1P.fastq.gz' \
    out2='$baseDir/fastq/MBWGS010_Trimmed_2P.fastq.gz' \
    pigz=t unpigz=t


########## Maping ##########


# Map reads of MBWG010 (2011/0690) onto this reference
genome="$baseDir/ffn/all_locus_tags.ffn"
r1='$baseDir/fastq/MBWGS010_Trimmed_1P.fastq.gz'
r2='$baseDir/fastq/MBWGS010_Trimmed_2P.fastq.gz'

bwa index "$genome"

#map reads
bwa mem -t "$cpu" -r 1 -a -M "$genome" "$r1" "$r2" | \
    samtools view -@ "$cpu" -b -h -F 4 - | \
    samtools sort -@ "$cpu" -o "$baseDir/fastq/MBWGS010.bam" -

#remove duplicates
samtools rmdup \
        "$baseDir/fastq/MBWGS010.bam" \
        "$baseDir/fastq/MBWGS010_nodup.bam"

#index bam file
samtools index "$baseDir/fastq/MBWGS010_nodup.bam"

# Inspect in IGV
igv


########## Assembly ##########

spades.py \
    --pe1-1 '$baseDir/fastq/MBWGS010_Trimmed_1P.fastq.gz' \
    --pe1-2 '$baseDir/fastq/MBWGS010_Trimmed_2P.fastq.gz' \
    --careful \
    -t "$cpu" \
    -o "$baseDir/assembly" \
    -k "21,33,55,77,99"


#Check contigs quickly with blast on nr to see if contamination
blastn \
    -query "$baseDir/assembly/scaffolds.fasta" \
    -db "nt" \
    -out "$baseDir/scaffolds.blastn" \
    -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
    -num_threads "$cpu" \
    -max_target_seqs 1



########## Polishing ##########


#Pilon
genome="$baseDir/assembly/scaffolds.fasta" 
r1='$baseDir/fastq/MBWGS010_Trimmed_1P.fastq.gz'
r2='$baseDir/fastq/MBWGS010_Trimmed_2P.fastq.gz'

#index reference genome
bwa index "$genome"

[ -d "$baseDir/pilon" ] || mkdir -p "$baseDir/pilon"

#map reads
bwa mem -t "$cpu" -r 1 -a -M "$genome" "$r1" "$r2" | \
    samtools view -@ "$cpu" -b -h -F 4 - | \
    samtools sort -@ "$cpu" -o "$baseDir/pilon/MBWGS010.bam" -

#remove duplicates
samtools rmdup \
        "$baseDir/pilon/MBWGS010.bam" \
        "$baseDir/pilon/MBWGS010_nodup.bam"

#index bam file
samtools index "$baseDir/pilon/MBWGS010_nodup.bam"

#Correct contigs using pilon based on the Illumina reads
java "$memJava" -jar "${prog}"/pilon/pilon-dev.jar \
    --threads "$cpu" \r
    --genome "$genome" \
    --bam "$baseDir/pilon/MBWGS010_nodup.bam" \
    --outdir "$baseDir/pilon" \
    --output "MBWGS010_pilon" \
    --changes \
    --verbose \
    | tee "$baseDir/pilon/pilon.log"


########## Annotation ##########


# Create output annotation folder if does not already exist
[ -d "$baseDir/annotation/" ] || mkdir -p "$baseDir/annotation/"

#Find and annotate genes with prokka
prokka \
    --force \
    --cpus "$cpu" \
    --rfam \
    --kingdom "Bacteria" \
    --genus "Mycobacteria" \
    --species "bovis" \
    --strain "MBWGS010" \
    --gram "pos" \
    --prefix "MBWGS010" \
    --compliant \
    --outdir "$baseDir/annotation/" \
    "$baseDir/pilon/MBWGS010_pilon.fasta"


########## Blasting ##########


# # Blast locus_tag sequences from NC_002945 on an annotated genome assembly of MBWGS010

# #create blast database with concatenated CDS from 
# makeblastdb \
#     -in "$baseDir/annotation/MBWGS010.ffn" \
#     -input_type "fasta" \
#     -dbtype "nucl" \
#     -parse_seqids -hash_index

# # Run blast to find reference genes in assembly. Keep only the best hits
# blastn \
#     -query ""${baseDir}"/ffn/all_locus_tags.ffn" \
#     -db ""${baseDir}"/annotation/MBWGS010.ffn" \
#     -out ""${baseDir}"/orfs.blastn" \
#     -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
#     -num_threads "$cpu" \
#     -max_target_seqs 1

# # only keep the first line of blast result is many for each query
# cat "${baseDir}"/orfs.blastn \
#     | sort -uk1,1 \
#     > "${baseDir}"/orfs_uniq.blastn

# # add header
# echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > "${baseDir}"/orfs.blastn.tmp
# cat "${baseDir}"/orfs.blastn >> "${baseDir}"/orfs.blastn.tmp
# mv "${baseDir}"/orfs.blastn.tmp "${baseDir}"/orfs.blastn

# echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > "${baseDir}"/orfs_uniq.blastn.tmp
# cat "${baseDir}"/orfs_uniq.blastn >> "${baseDir}"/orfs_uniq.blastn.tmp
# mv "${baseDir}"/orfs_uniq.blastn.tmp "${baseDir}"/orfs_uniq.blastn

# # make ffn file "2-line" style
# cat "${baseDir}"/annotation/MBWGS010.ffn \
#     | awk '{if(substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' \
#     > "${baseDir}"/annotation/MBWGS010_2lines.ffn

# # terms to search
# cat "${baseDir}"/orfs_uniq.blastn \
#     | sed 1d \
#     | cut -f 2 \
#     > "${baseDir}"/list.txt

# # get subject header and sequences from list
# cat "${baseDir}"/annotation/MBWGS010_2lines.ffn \
#     | LC_ALL=C grep --no-group-separator -w -A 1 -F -f ""${baseDir}"/list.txt" \
#     > "${baseDir}"/hits.ffn

# # Add query to hit header
# for i in $(cat ""${baseDir}"/list.txt"); do
#     query=$(cat "${baseDir}"/orfs_uniq.blastn | grep -w -F "$i" | cut -f 1)
#     subject=$(cat "${baseDir}"/annotation/MBWGS010_2lines.ffn \
#                 | LC_ALL=C grep --no-group-separator -w -A 1 -F "$i")
#     header=$(echo "$subject" | tr '\n' '@' | cut -d '@' -f 1)
#     seq=$(echo "$subject" | tr '\n' '@' | cut -d '@' -f 2)
#     echo -e ""$header" ("$query")\n"$seq"" >> "${baseDir}"/hits_final.ffn
# done

# # Get the protein sequences
# cat "${baseDir}"/hits_final.ffn \
#     | grep -E "^>" \
#     | cut -d " " -f 1 \
#     | tr -d ">" \
#     > "${baseDir}"/toFetch.list

# # make faa file "2-line" style
# cat "${baseDir}"/annotation/MBWGS010.faa \
#     | awk '{if(substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' \
#     > "${baseDir}"/annotation/MBWGS010_2lines.faa

# # get proteing sequences in the list
# cat "${baseDir}"/annotation/MBWGS010_2lines.faa \
#     | LC_ALL=C grep --no-group-separator -w -A 1 -F -f ""${baseDir}"/toFetch.list" \
#     > "${baseDir}"/hits.faa

# # Add query to hit header
# for i in $(cat ""${baseDir}"/toFetch.list"); do
#     query=$(cat "${baseDir}"/orfs_uniq.blastn | grep -w -F "$i" | cut -f 1)
#     subject=$(cat "${baseDir}"/annotation/MBWGS010_2lines.faa \
#                 | LC_ALL=C grep --no-group-separator -w -A 1 -F "$i")
#     header=$(echo "$subject" | tr '\n' '@' | cut -d '@' -f 1)
#     seq=$(echo "$subject" | tr '\n' '@' | cut -d '@' -f 2)
#     echo -e ""$header" ("$query")\n"$seq"" >> "${baseDir}"/hits_final.faa
# done


#manually remove one hit for a tRNA, which had no protein sequence

###################################################################################################################################
# An attempt to parallelize the blast code - Work in progress


bash ~/scripts/refseqDownloader.sh \
    -t "bacteria" \
    -q "Mycobacterium bovis" \
    -o /home/bioinfo/Desktop/nadia/ECSVM_positive/refseq_mbovis \
    -u -n 48

function blast ()
{
    ref="$1"
    query="$2"

    db=$(basename "$ref")
    # db_name="${db%.ffn}"
    db_name="${db%.faa}"

    makeblastdb \
        -in "$ref" \
        -input_type "fasta" \
        -dbtype "prot" \
        -parse_seqids -hash_index

    out="/home/bioinfo/Desktop/nadia/ECSVM_positive/blast/"${db_name}"_ECSVM_hits.blastn.tsv"
    outfmt="6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    blastp \
        -query "$query" \
        -db "$ref" \
        -out "$out" \
        -outfmt "$outfmt" \
        -num_threads 4 \
        -max_target_seqs 1

    # only keep the first line of blast result is many for each query
    cat "$out" \
        | sort -uk1,1 \
        > "${out%.blastn.tsv}"_uniq.blastn.tsv

    # add header
    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > "${out}".tmp
    cat "$out" >> "${out}".tmp
    mv "${out}".tmp "$out"

    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > "${out%.blastn.tsv}"_uniq.blastn.tsv.tmp
    cat "${out%.blastn.tsv}"_uniq.blastn.tsv >> "${out%.blastn.tsv}"_uniq.blastn.tsv.tmp
    mv "${out%.blastn.tsv}"_uniq.blastn.tsv.tmp "${out%.blastn.tsv}"_uniq.blastn.tsv
}

export -f blast

find /home/bioinfo/Desktop/nadia/ECSVM_positive/refseq_mbovis/faa -type f -name "*.faa" |
parallel    --bar \
            --env blast \
            --jobs 12 \
            'blast {} /home/bioinfo/Desktop/nadia/ECSVM_positive/faa/all_locus_tags.faa'


function crossref_ffn2faa ()
{
    blast="$1"
    name=$(sed 's/_ECSVM.*//' <<< $(basename "$blast"))
    ref_ffn="/home/bioinfo/Desktop/nadia/ECSVM_positive/refseq_mbovis/ffn/"${name}".ffn"
    ref_faa="/home/bioinfo/Desktop/nadia/ECSVM_positive/refseq_mbovis/faa/"${name}".faa"

    # make ffn file "2-line" style
    cat "$ref_ffn" \
        | awk '{if(substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' \
        > "${ref_ffn%.ffn}"_2lines.ffn

    # terms to search
    cat "$blast" \
        | sed 1d \
        | cut -f 2 \
        | sed 's/^.*cds_//' \
        | grep -F '.' \
        | sed 's/\.1.*/\.1/' \
        > "${blast%.blastn.tsv}".list

    # get subject header and sequences from list
    cat "${ref_ffn%.ffn}"_2lines.ffn \
        | LC_ALL=C grep --no-group-separator -w -A 1 -F -f "${blast%.blastn.tsv}".list \
        > "${blast%.blastn.tsv}".hits.ffn

    # Add query to hit header
    [ -e "${blast%.blastn.tsv}".hits.final.ffn ] && rm -rf "${blast%.blastn.tsv}".hits.final.ffn
    for i in $(cat "${blast%.blastn.tsv}".list); do
        query=$(cat "$blast" | grep -w -F "$i" | cut -f 1)
        subject=$(cat "${ref_ffn%.ffn}"_2lines.ffn \
                    | LC_ALL=C grep --no-group-separator -w -A 1 -F "$i")
        header=$(echo "$subject" | tr '\n' '@' | cut -d '@' -f 1)
        seq=$(echo "$subject" | tr '\n' '@' | cut -d '@' -f 2)
        echo -e ""$header" ("$query")\n"$seq"" >> "${blast%.blastn.tsv}".hits.final.ffn
    done

    # Get the protein sequences

    #get accession numbers
    cat "${blast%.blastn.tsv}".hits.final.ffn \
        | grep -E "^>" \
        | cut -d " " -f 1 \
        | tr -d ">" \
        | sed 's/^.*cds_//' \
        | grep -F '.' \
        | sed 's/\.1.*/\.1/' \
        > "${blast%.blastn.tsv}".tofetch

    # make faa file "2-line" style
    cat "$ref_faa" \
        | awk '{if(substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' \
        > "${ref_faa%.faa}"_2lines.faa

    # get protein sequences in the list
    cat "${ref_faa%.faa}"_2lines.faa \
        | LC_ALL=C grep --no-group-separator -w -A 1 -F -f "${blast%.blastn.tsv}".tofetch \
        >  "${blast%.blastn.tsv}".hits.faa

    # Add query to hit header
    [ -e "${blast%.blastn.tsv}".hits.final.faa ] && rm -rf "${blast%.blastn.tsv}".hits.final.faa
    for i in $(cat "${blast%.blastn.tsv}".tofetch); do
        query=$(cat "$blast" | grep -w -F "$i" | cut -f 1)
        subject=$(cat "${ref_faa%.faa}"_2lines.faa \
                    | LC_ALL=C grep --no-group-separator -w -A 1 -F "$i")
        header=$(echo "$subject" | tr '\n' '@' | cut -d '@' -f 1)
        seq=$(echo "$subject" | tr '\n' '@' | cut -d '@' -f 2)
        echo -e ""$header" ("$query")\n"$seq"" >> "${blast%.blastn.tsv}".hits.final.faa
    done
}

export -f crossref_ffn2faa

find /home/bioinfo/Desktop/nadia/ECSVM_positive/blast -type f -name "*_uniq.blastn.tsv" | \
parallel    --bar \
            --env crossref_ffn2faa \
            'crossref_ffn2faa {}'


function crossref_faa2faa ()
{
    blast="$1"
    name=$(sed 's/_ECSVM.*//' <<< $(basename "$blast"))
    ref_faa="/home/bioinfo/Desktop/nadia/ECSVM_positive/refseq_mbovis/faa/"${name}".faa"

    #get accession numbers
    cat "$blast" \
        | sed 1d \
        | cut -f 2 \
        | tr -d "|" \
        | sed 's/ref//' \
        > "${blast%.blastn.tsv}".list

    # make faa file "2-line" style
    cat "$ref_faa" \
        | awk '{if(substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' \
        > "${ref_faa%.faa}"_2lines.faa

    # get protein sequences in the list
    cat "${ref_faa%.faa}"_2lines.faa \
        | LC_ALL=C grep --no-group-separator -w -A 1 -F -f "${blast%.blastn.tsv}".list \
        >  "${blast%.blastn.tsv}".hits.faa

    # Add query to hit header
    [ -e "${blast%.blastn.tsv}".hits.final.faa ] && rm -rf "${blast%.blastn.tsv}".hits.final.faa
    for i in $(cat "${blast%.blastn.tsv}".list); do
        query=$(cat "$blast" | grep -w -F "$i" | cut -f 1)
        subject=$(cat "${ref_faa%.faa}"_2lines.faa \
                    | LC_ALL=C grep --no-group-separator -w -A 1 -F "$i")
        header=$(echo "$subject" | tr '\n' '@' | cut -d '@' -f 1)
        seq=$(echo "$subject" | tr '\n' '@' | cut -d '@' -f 2)
        echo -e ""$header" ("$query")\n"$seq"" >> "${blast%.blastn.tsv}".hits.final.faa
    done
}

export -f crossref_faa2faa

find /home/bioinfo/Desktop/nadia/ECSVM_positive/blast -type f -name "*_uniq.blastn.tsv" | \
parallel    --bar \
            --env crossref_faa2faa \
            'crossref_faa2faa {}'



###################################################################################################################################


########## Localization Prediction ##########


#run tmhmm
decodeanhmm.Linux_x86_64 \
    -f /home/bioinfo/prog/tmhmm-2.0c/lib/TMHMM2.0.options \
    /home/bioinfo/prog/tmhmm-2.0c/lib/TMHMM2.0.model \
    "${baseDir}"/hits_final.faa \
    > "${baseDir}"/hits_final.faa.tmhmm.txt

#run SignalP
signalp \
    -t gram+ \
    "${baseDir}"/hits_final.faa \
    > "${baseDir}"/hits_final.faa.signalp.txt

#run Phobius
perl /home/bioinfo/prog/interproscan-5.24-63.0/bin/phobius/1.01/phobius.pl \
    "${baseDir}"/hits_final.faa \
    > "${baseDir}"/hits_final.faa.phobius.txt

perl /home/bioinfo/prog/phobius/jphobius -r \
    "${baseDir}"/hits_final.faa \
    > "${baseDir}"/hits_final.faa.jphobius.txt

#ran PSORTb online manualy
# "${baseDir}"/hits_final.faa.psortb.txt
# reformat output
cat hits_final.faa.psortb.txt \
    | awk '/SeqID/ {print}; /Final/ {getline; print}; /Secondary/ {getline; print}' \
    | sed 's/SeqID: /@/g' \
    | tr -d "\r" \
    | tr -s " " \
    | awk '{if (!/^@/) {gsub("^ ", ""); gsub(" ", "\t"); print} else {print}}' \
    | tr "\n" "\t" \
    | tr "@" "\n" \
    | sed 's/[[:space:]]*$//' \
    | sed 's/ /\t/' \
    > hits_final.faa.psortb_short.txt
