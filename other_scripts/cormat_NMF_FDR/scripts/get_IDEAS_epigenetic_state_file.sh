cd /working_directory
# download the epigenetic state file from this link: https://usevision.org/data/hg38/IDEASstates/ideasJointMay2021/


# convert space serperated file to tab serperated file and remove the first column and last column and 40th and 41st column (remove NEU)
tail -n+2 S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state | cut -d ' ' -f 2- | rev | cut -d ' ' -f 2- | rev \
| sed 's/ /\t/g' | cut -f 1-38,41- > S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state.bed


