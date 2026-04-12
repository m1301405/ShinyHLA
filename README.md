# ShinyHLA: an R/Shiny-based web application for HLA typing of NGS data

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

This integration allows systematic comparison of HLA typing results across different algorithms.  
These tools are bundled inside the ShinyHLA Docker environment and do not require manual installation.

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

# Documentation

Full documentation for ShinyHLA is available online:

**Documentation**  
https://animated-reminder-9ff.notion.site/ebd//7b6f9fcf9d034442b2b64c89f8662a6c

The documentation includes:

- Framework implementation details  
- Step-by-step user guides  
- Demonstration using HG00096 datasets  
- PharmGKB annotation exploration  
- Docker deployment instructions  
- Preprocessing workflow with GUI screenshots  

---

# Quick Start (Docker)

ShinyHLA is distributed as a Docker container to ensure reproducibility and simplify installation.

### Pull the Docker image

```bash
docker pull b0909132/shinyhla:latest
```

### Run the container

```bash
docker run -p 3838:3838 b0909132/shinyhla:latest
```

### Access the ShinyHLA interface

Open a browser and navigate to:

```text
http://localhost:3838
```

---

# Optional Preprocessing

To improve computational efficiency, ShinyHLA provides an optional preprocessing workflow that extracts reads mapped to the HLA genomic region before HLA typing.  
This step reduces the input data size while retaining reads relevant to HLA loci, thereby accelerating downstream analyses.

The preprocessing module is provided as a separate Docker-based web application.

### 1. Pull the preprocessing image

```bash
docker pull b0909132/shinyhla_preprocessing
```

### 2. Start the preprocessing container in Docker Desktop

After the image is downloaded:

- Open **Docker Desktop**
- Go to the **Images** tab
- Locate `b0909132/shinyhla_preprocessing`
- Click **Run**
- Assign a local port (for example, `3838`) that does not conflict with other running services
- Start the container

### 3. Open the preprocessing web application

Once the container status shows that the service is running, open the corresponding local link in your browser, for example:

```text
http://localhost:3838
```

### 4. Perform preprocessing in the GUI

Within the preprocessing web interface:

- Select the input data format (**BAM** or **CRAM**)
- Upload the sequencing file together with its index file
- Click **Upload**
- Click **Extraction** to process the file
- Click **Download** to obtain the extracted result

The processed output file can then be uploaded to ShinyHLA for downstream HLA typing.

Detailed GUI screenshots and step-by-step instructions are available in the online documentation.

---

# Supported Input Data

ShinyHLA supports the following input formats:

| NGS Type | Supported Formats |
|--------|------------------|
| Whole Exome Sequencing (WES) | FASTQ, BAM |
| RNA Sequencing (RNA-seq) | FASTQ, BAM |

BAM files should be aligned to the **GRCh38 reference genome** for compatibility with the Quality Control module.

---

# Built-in Demonstration Dataset

ShinyHLA includes built-in demonstration datasets for interactive exploration.

The demo data are based on the **HG00096 sample**:

- Whole-exome sequencing (WES) from the **1000 Genomes Project**  
- RNA-seq data from the **Geuvadis consortium**

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
