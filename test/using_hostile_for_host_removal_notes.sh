hostile clean --fastq1 ./ERR4686055_1.fastq --fastq2 ERR4686055_2.fastq \
--aligner bowtie2 \
--index /home/ac.jsheriff/dev_container/projects/red_fox/ncbi_dataset/data/GCF_048418805.1/GCF_048418805.1_VulVul3_index*


datasets download genome accession GCF_000001405.40 --include cds,protein,gff3

for f in *.zip; do 
    unzip -o "$f" && rm "$f"; 
done

nohup /home/ac.jsheriff/dev_container/projects/host_removal/scripts/create_custom_bowtie2_indices.sh \
    > /home/ac.jsheriff/dev_container/projects/host_removal/files/create_custom_bowtie2_indices.log 2>&1 &

ps aux | grep ./scripts/create_custom_bowtie2_indices.sh 


# Run Test scripts

App-FastqUtils xx /home/ac.jsheriff/dev_container/modules/p3_fqutils/app_specs/FastqUtils.json \
    /home/ac.jsheriff/dev_container/projects/host_removal/test/Fastq_Utils_Test/bowtie2/ERR12389203_zebrafish_bowtie2/ERR12389203_zebrafish_bowtie2.json


App-FastqUtils xx /home/ac.jsheriff/dev_container/modules/p3_fqutils/app_specs/FastqUtils.json \
    /home/ac.jsheriff/dev_container/projects/host_removal/test/Fastq_Utils_Test/minimap2/SRR39975832_cat_minimap2/input_param.json




# copy results to test with Kraken2

p3-cp /home/ac.jsheriff/dev_container/projects/host_removal/test/Fastq_Utils_Test/bowtie2/SRR40478378_Pig_bowtie2/work/SRR40478378* ws:"/jsheriff@bvbrc/BVBRC_UIC_Test"

p3-mkdir "/jsheriff@bvbrc/BVBRC_UIC_Test"

die "Failure uploading $local_file to shock: " . $sres->status_line;


For test - create synthetic FASTQ to test on alpha - look at IDs after filtering to see if any reads still exist

js_run_dev_tests \
      -s p3_fqutils \
      -t /home/ac.jsheriff/dev_container/modules/p3_fqutils/test/hostile_tests \
      -r /home/ac.jsheriff/dev_container/modules \
      -o /jsheriff@bvbrc/Jamal_Dev/Dev_Tests/FastqUtils_Hostile_DevTest_30Oct2026 \
      -d /home/ac.jsheriff/dev_container/dev_tests/FastqUtils_Hostile_DevTest_30Oct2026