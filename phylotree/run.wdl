version 1.0

task idseq_main {
  runtime {
    docker: docker_image_id
  }
  input {
    String docker_image_id
    String s3_wd_uri
    String dag_branch
    String message
    File fastqs_0
    File? fastqs_1
  }
  command<<<
  echo ~{message}
  echo Phylotree workflow goes here!
  >>>
  output {
    Array[String]+ messages = read_lines(stdout())
  }
}
