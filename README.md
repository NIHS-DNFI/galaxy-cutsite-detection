# Galaxy for Cut Site Detection
Galaxy for Cut Site Detection is a Galaxy (https://usegalaxy.org/) based data analysis pipeline for off-target sequence such as SITE-Seq. 
Galaxy for Cut Site Detection identifies DSB sites with termination positions of aligned reads. 
Optionally, user can focus the cut sites only within the annotated region (e.g. exon).


## Usage
### Requirements
If you are using Windows machine, you must first enable the “Windows Subsystem for Linux” optional feature before installing any Linux distributions on Windows. Please refer to the following page for WSL installation (https://learn.microsoft.com/en-us/windows/wsl/install).
Galaxy for Cut Site Detection requires Docker. Download from Docker Desktop official site (https://www.docker.com/products/docker-desktop).

<br>

### Running Galaxy for Cut Site Detection container
First, pull a Galaxy for Cut Site Detection docker image from our repository (https://hub.docker.com/repository/docker/nihsdnfi/galaxy-cutsite-detection).
```vb
$ docker pull nihsdnfi/galaxy-cutsite-detection:0.XX
```
<br>

### Run Galaxy for Cut Site Detection  
```vb
$ docker run -d --name galaxy-csd -p 8080:80 nihsdnfi/galaxy-cutsite-detection:0.XX
```
After a few minutes of executing above command, you will be able to access Galaxy for Cut Site Detection via ```localhost:8080``` from your web browser.
You will need to launch your browser manually and connect to ```localhost:8080```.

<br>

## Documents
More detailed operation manual can be found in manual.pdf.


## License
Galaxy for Cut Site Detection is licensed under MIT License. 


## Citation
Jumpei Narushima, Shinya Kimata, Yuh Shiwa, Takahiro Gondo, Satoru Akimoto, Keisuke Soga, Satoko Yoshiba, Kosuke Nakamura, Norihito Shibata and Kazunari Kondo. (2021). bioRxiv. Strategy for detecting off-target sites in genome-edited rice. doi: https://doi.org/10.1101/2021.05.28.446070
