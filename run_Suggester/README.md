This is a python script to run Suggester.

It takes three pasers from command as,
```
python run_Suggester/run_Suggester.py {directory of suggester}/ncbi-sra-run-suggester/ncbi-sra-run-suggester.jar {directory of input data} {directory to save output file}
```

This script also remove any information other than the SRA IDs, and automatically name the output file in a style of
yyyymmdd_yyyymmdd.txt. It also outputs log files for any stdrr.