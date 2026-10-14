mkdir -p ~/lsf_logs/compress
bsub -q long \
     -n 8 \
     -M 8000 \
     -R "select[mem>=8000] rusage[mem=8000] span[hosts=1]" \
     -G team222 \
     -J "compress_results" \
     -o ~/lsf_logs/compress/compress_%J.out \
     -e ~/lsf_logs/compress/compress_%J.err \
     ~/git_repos/mosquito_virome_yucatan_LEVE/scripts/individual_analyses/compress_results.sh
