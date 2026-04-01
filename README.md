# ShinyHLA: an R/Shiny-based web application for HLA typing of NGS data
ShinyHLA is an R/Shiny-based tool for HLA typing that supports WES and RNA-seq data, integrates leading HLA typing tools and reference databases, and provides a user-friendly graphical interface for result exploration and analysis.

- *Data compatibility: Supports HLA typing for both WES and RNA-seq data.*
- *Integration of tools and databases: Incorporates OptiType, arcasHLA, HLA-HD, and SpecHLA, together with the IPD-IMGT/HLA reference database.*
- *User-friendly interface: Provides a graphical user interface (GUI) that eliminates the need for programming experience.*
- *Enhanced result interpretation: Integrates pivot-table-based comparison and PharmGKB annotations to facilitate interpretation of HLA typing results in pharmacogenomic contexts.*

## Website
The official web documentation for ShinyHLA is available online.

- Documentation: [ShinyHLA Documentation](https://animated-reminder-9ff.notion.site/ebd//7b6f9fcf9d034442b2b64c89f8662a6c)

It provides:

1. **Step-by-step guides:** Detailed instructions for tasks such as file upload, preprocessing, quality control, and HLA typing.
2. **Tutorials:** Clear walkthroughs for both beginners and advanced users, demonstrating core features such as HLA typing, result exploration, PharmGKB-based pharmacogenomic annotation, and pivot table usage for cross-tool result comparison.
3. **Docker deployment instructions:** Simple and clear guidelines for setting up Docker environments to optimize preprocessing and streamline the analysis workflow.

## Requirements
ShinyHLA uses the following software and reference resources:

### 1. HLA typing tools
- [OptiType v1.3.5](https://github.com/FRED-2/OptiType)
- [arcasHLA v0.6.0](https://github.com/RabadanLab/arcasHLA)
- [HLA-HD v1.7.0](https://w3.genome.med.kyoto-u.ac.jp/HLA-HD/)
- [T1K v1.0.8](https://github.com/mourisl/T1K)

### 2. R packages
- Please refer to [ShinyHLA Documentation](https://animated-reminder-9ff.notion.site/ebd//7b6f9fcf9d034442b2b64c89f8662a6c), section 1.7.

### 3. Sequence Processing Tools
- [seqtk v1.3-r106](https://github.com/lh3/seqtk)
- [PanDepth v2.25](https://github.com/HuiyangYu/PanDepth)

### 4. Human genome reference
- [Gencode](https://www.gencodegenes.org/human/release_38.html): gencode.v38.annotation.gff3

## Installation via Docker

1. Install Docker on your computer and ensure that it is properly configured.
2. Run `docker pull b0909132/shinyhla` to download the Docker image.
3. Launch the container using:  
   `docker run -d -p 3838:3838 --name shinyhla b0909132/shinyhla:latest`
4. Open your web browser and navigate to:  
   `http://localhost:3838`.

For more detailed installation and usage instructions, please refer to [ShinyHLA Documentation](https://animated-reminder-9ff.notion.site/ebd//7b6f9fcf9d034442b2b64c89f8662a6c), section 3.

## Getting help
If you have any questions, please feel free to contact us. We will respond as soon as possible.

([m1301405@cgu.edu.tw](mailto:m1301405@cgu.edu.tw))
