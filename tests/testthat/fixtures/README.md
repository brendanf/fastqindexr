Fixture regeneration commands (run from package root):

```sh
Rscript -e "con <- gzfile('tests/testthat/fixtures/cli_fixture.fastq.gz','wt'); writeLines(c('@r1','ACGT','+','!!!!','@r2','TTAA','+','####','@r3','GGGG','+','$$$$','@r4','CCCC','+','%%%%'), con); close(con)"
"../lifeplan_COI/bin/fastqindex_0.9.0b" index -f="tests/testthat/fixtures/cli_fixture.fastq.gz" -i="tests/testthat/fixtures/cli_fixture.fastq.gz.fqi"

Rscript -e "con <- gzfile('tests/testthat/fixtures/cli_fixture_part2.fastq.gz','wt'); writeLines(c('@r5','AAAA','+','++++','@r6','TTTT','+','----'), con); close(con)"
"../lifeplan_COI/bin/fastqindex_0.9.0b" index -f="tests/testthat/fixtures/cli_fixture_part2.fastq.gz" -i="tests/testthat/fixtures/cli_fixture_part2.fastq.gz.fqi"
```
