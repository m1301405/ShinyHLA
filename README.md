# ShinyHLA: an R/Shiny-based web application for HLA typing of NGS data

[![Docker](https://img.shields.io/badge/Docker-b0909132%2Fshinyhla-blue)](https://hub.docker.com/r/b0909132/shinyhla)
![License](https://img.shields.io/badge/license-MIT-green)
[![Documentation](https://img.shields.io/badge/docs-Notion-blue)](https://animated-reminder-9ff.notion.site/ebd//7b6f9fcf9d034442b2b64c89f8662a6c)

ShinyHLA is an R/Shiny-based web application designed to simplify HLA typing analysis from next-generation sequencing (NGS) data.  
The platform integrates multiple widely used HLA typing algorithms and supports multiple versions of the IPD-IMGT/HLA reference database within a unified graphical interface.

ShinyHLA enables users to perform HLA typing, compare results across tools and database versions, and explore pharmacogenomic annotations through an interactive and reproducible workflow.

---

# Key Features

### Multi-tool HLA typing integration

ShinyHLA integrates several widely used HLA typing algorithms within a single graphical interface:

- [OptiType v1.3.5](https://github.com/FRED-2/OptiType)
- [arcasHLA v0.6.0](https://github.com/RabadanLab/arcasHLA)
- [HLA-HD v1.7.0](https://w3.genome.med.kyoto-u.ac.jp/HLA-HD/)
- [T1K v1.0.8](https://github.com/mourisl/T1K)

These tools are bundled inside the ShinyHLA Docker environment and **do not require manual installation**.

---

### Multi-version IPD-IMGT/HLA database support

ShinyHLA supports multiple versions of the IPD-IMGT/HLA database, allowing users to evaluate potential differences in allele assignments caused by database updates.

---

### User-friendly graphical interface

ShinyHLA provides a graphical user interface implemented in **R Shiny**, enabling users to perform HLA typing analyses without command-line experience.

---

### Cross-analysis results integration

Typing results generated from different tools, NGS data types, and database versions can be consolidated using an interactive **pivot table**.

This allows flexible comparison across:

- HLA typing tools  
- NGS data types  
- MHC classes  
- IPD-IMGT/HLA versions  

Results can be exported as `.xlsx` files for further analysis.

---

### Integrated PharmGKB pharmacogenomic annotations

ShinyHLA integrates pharmacogenomic knowledge from the **PharmGKB** database, enabling users to explore clinically relevant associations between HLA alleles and drug responses.

Three annotation modules are available:

- Summary annotations  
- Variant annotations  
- Drug label annotations  

---

# Online Web Application

ShinyHLA is freely available as an online web server for testing and demonstration:

https://shiny.cgu.edu.tw/app/shinyHLA

The web server allows users to perform HLA typing directly through a graphical interface without installing any software.

Users can:

- Upload WES or RNA-seq data  
- Perform HLA typing using multiple tools  
- Compare results across IPD-IMGT/HLA database versions  
- Explore PharmGKB pharmacogenomic annotations  
- Export typing results and integrated comparison tables  

---

# Local Deployment (Docker)

For users who prefer to run ShinyHLA locally, the application is available as a Docker container.

### 1. Pull the Docker image

```bash
docker pull b0909132/shinyhla:latest
```

### 2. Run the container

```bash
docker run -p 3838:3838 --name shinyhla b0909132/shinyhla:latest
```

### 3. Access the ShinyHLA interface

Open a browser and navigate to:

```
http://localhost:3838
```

Running ShinyHLA locally allows users to analyze large datasets and maintain full control over their computing environment.

---

# Optional Preprocessing

To improve computational efficiency, ShinyHLA provides an optional preprocessing workflow that extracts reads mapped to the HLA genomic region before HLA typing.

```
chr6:29602228–33410226
```

This preprocessing step reduces the input data size while retaining reads relevant to HLA loci, thereby accelerating downstream HLA typing analyses.

The preprocessing module is provided as a separate Docker-based web application.

### 1. Pull the preprocessing image

```bash
docker pull b0909132/shinyhla_preprocessing:latest
```

### 2. Launch the preprocessing container

After downloading the image:

- Open **Docker Desktop**
- Navigate to the **Images** tab
- Locate `b0909132/shinyhla_preprocessing`
- Click **Run**
- Assign a local port (e.g., `3839`)
- Start the container

### 3. Access the preprocessing web interface

Once the container is running, open:

```
http://localhost:3839
```

### 4. Perform preprocessing through the GUI

Within the preprocessing interface:

- Select the sequencing data format (**BAM** or **CRAM**)
- Upload the sequencing file and index file
- Click **Upload**
- Click **Extraction**
- Click **Download** to obtain the processed file

The processed output file can then be uploaded to ShinyHLA for downstream HLA typing.

Detailed GUI screenshots and instructions are available in the online documentation.

---

# Documentation

Full documentation for ShinyHLA is available online:

https://animated-reminder-9ff.notion.site/ebd//7b6f9fcf9d034442b2b64c89f8662a6c

The documentation includes:

- Implementation details of the ShinyHLA framework  
- Step-by-step user guides  
- Demonstration using HG00096 datasets  
- PharmGKB annotation exploration  
- Docker deployment instructions  
- Preprocessing workflow  

---

# Built-in Demonstration Dataset

ShinyHLA includes built-in demonstration datasets for interactive exploration.

The demo data are based on the **HG00096 sample**:

- Whole-exome sequencing (WES) from the **1000 Genomes Project**  
- RNA-seq data from the **Geuvadis consortium**

Precomputed typing results are included for multiple IPD-IMGT/HLA versions (**v3.14.0 and v3.49–v3.60**) to enable instant demonstration without running the full analysis pipeline.

These datasets demonstrate:

- HLA typing using multiple tools  
- Cross-tool comparison  
- Database version comparison  
- PharmGKB pharmacogenomic annotation exploration  

---

# Contact

For questions, feedback, or bug reports, please contact:

m1301405@cgu.edu.tw

---

# License

This project is released under the MIT License.
