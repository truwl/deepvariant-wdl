import "htsget_DeepVariant_GLnexus.wdl" as main

workflow DVGLx_test {  
  call main.htsget_DeepVariant_GLnexus
  call check_outputs {
    input:
      pvcf_gz = htsget_DeepVariant_GLnexus.pvcf_gz
  }
}

task check_outputs {
  File pvcf_gz

  command {
    set -ex -o pipefail
    [ -s "${pvcf_gz}" ] || exit 1
  }
}

