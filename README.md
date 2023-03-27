# HarvestVariants_AutoPipe

Inside the root directory, run ```hv_pipleline.py``` in command:
```
python hv_pipeline.py --email {email associated with your NCBI account} --api_ky {api key of your NCBI account} --pdat {yyyy/mm/dd:yyyy/mm/dd} --suggester_jar {path/to/ncbi-sra-run-suggester.jar}
```

The ```pdat``` is a required parser for the start and end dates of publication date. The ```suggester_jar``` requires the path to where you store the ```run-suggester.jar``` in your local or server directory.

The output files including intermediate output files will be instored in the follwing file tree stucture where the directories will be automatically created if not exist,
```bash
├── Data
│   ├── hv
│   ├── pango
│   ├── suggester
│   └── xml
```
where the xml files (including both the original files and pretty print files converted after xml_pp tool will be stored in the `xml` folder, and the SRA ID list txt files will be stored in the suggester folder, etc.)

The email and api_key parameters are optional. You don't need an email address of registered NCBI account or the API Key associated with the account to run the tool.  
However, it's highly recommended to do so as NCBI Entrez gives permission to access up to 10 requests per second to E-utilities  
with associated email and api key but only 3 requests without them.

Please also be considered when you need to access large amounts of requests. NCBI recommends that users do so on the weekends  
or between 9:00 PM and 5:00 AM Eastern time during weekdays. For more info, please visit: https://www.ncbi.nlm.nih.gov/books/NBK25497/
