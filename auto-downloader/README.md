This is an auto downloading tool using Biopython Entrez module. For more information please visit:  
https://biopython.org/docs/1.75/api/Bio.Entrez.html

You don't need an email address of registered NCBI account or the API Key associated with the account to run the tool.  
However, it's highly recommended to do so as NCBI Entrez gives permission to access up to 10 requests per second to E-utilities  
with associated email and api key but only 3 requests without them.

Please also be considered when you need to access large amounts of requests. NCBI recommends that users do so on the weekends  
or between 9:00 PM and 5:00 AM Eastern time during weekdays. For more info, please visit: https://www.ncbi.nlm.nih.gov/books/NBK25497/

This tool needs to be run under python 3 or above and biopython enviornments. To run the tool,
```
python auto_downloader/downloader.py --email {your NCBI email address} --api_key {your NCBI api key} --pdat yy/mm/dd:yy/mm/dd --output_dir {directory to save file}
```
The pdat is a required parser for the start and end dates of publication date. Email and api_key are optional. The output location is also optional 
and if you do not specific a location, it will be saved in your current directory where you run the tool. Default filename is in yyyymmdd_yyyymmdd_original.xml 
style which indicates the start and end pdat.