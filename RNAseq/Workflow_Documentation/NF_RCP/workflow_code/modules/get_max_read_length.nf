process GET_MAX_READ_LENGTH {
  // tag("Dataset-wide")
  
  input:
    path(fastqc_datazips)

  output:
    env(MAX_LENGTH), emit: length

  script:
    """
    MAX_LENGTH=0
    for zip in $fastqc_datazips
    do
      unzip -p \$zip '*/fastqc_data.txt' | awk '/Sequence length/ {split(\$0,a,"\t"); if(a[2]~/-/){split(a[2],b,"-"); len=b[2]}else{len=a[2]}; if(len>max){max=len}} END{print max}' >> lengths.txt
    done
    MAX_LENGTH=\$(sort -rn lengths.txt | head -n 1)
    """
}